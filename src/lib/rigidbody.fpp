!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rigidbody_mod
!> @brief   setup and numerical primitives for rigid-body dynamics
!! @authors Genesis Developers
!
!  (c) Copyright 2026 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module rigidbody_mod

  use rigidbody_str_mod
  use fileio_rigidbody_mod
  use molecules_str_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! subroutines
  public :: setup_rigidbody
  public :: rigidbody_diameter
  public :: quat_to_rotmatrix
  public :: quat_normalize
  public :: quat_derivative
  public :: rigidbody_angvel
  public :: fit_rigidbody_quat
  public :: quat_multiply
  public :: quat_rotate_by_body_omega

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rigidbody
  !> @brief        build rigid-body information (membership, mass, inertia,
  !!               principal-axis reference geometry) from a global,
  !!               non-decomposed molecule
  !! @param[in]    molecule  : molecule information
  !! @param[in]    rb_index  : rigid-body atom-index groups
  !! @param[in]    rb_ref    : rigid-body body-fixed reference coordinates
  !! @param[out]   rigidbody : rigid-body information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rigidbody(molecule, rb_index, rb_ref, rigidbody)

    ! formal arguments
    type(s_molecule),          intent(in)    :: molecule
    type(s_rigidbody_index),   intent(in)    :: rb_index
    type(s_rigidbody_ref),     intent(in)    :: rb_ref
    type(s_rigidbody),         intent(inout) :: rigidbody

    ! local variables
    integer                    :: ib, i, n, iatm, natom, nbody
    real(wp)                   :: total_mass, m
    real(wp)                   :: com_ref(3)
    real(wp)                   :: r(3)
    real(wp)                   :: inertia_mat(3,3)
    real(wp)                   :: eval(3), evec(3,3), work(24)
    integer                    :: ierr


    nbody = rb_index%num_groups
    natom = rb_ref%num_atoms

    do ib = 1, nbody
      if (size(rb_index%group(ib)%idx) /= natom) &
        call error_msg('Setup_Rigidbody> group size does not match '// &
                        'the reference coordinate file')
    end do

    call alloc_rigidbody(rigidbody, RigidBodyList, nbody, natom)
    call alloc_rigidbody(rigidbody, RigidBodyPhys, nbody, natom)

    rigidbody%is_used    = .true.
    rigidbody%num_bodies = nbody
    rigidbody%max_natom  = natom

    do ib = 1, nbody

      n = natom
      rigidbody%natom(ib) = n

      do i = 1, n
        iatm = rb_index%group(ib)%idx(i)
        if (iatm < 1 .or. iatm > molecule%num_atoms) &
          call error_msg('Setup_Rigidbody> atom index out of range in '// &
                          'rigid-body index file')
        rigidbody%atomlist(i,ib)    = iatm
        ! captured here (rather than looked up at runtime) because
        ! s_molecule is deallocated after spdyn setup completes
        rigidbody%atom_mass(i,ib)   = molecule%mass(iatm)
        rigidbody%atom_charge(i,ib) = molecule%charge(iatm)
        rigidbody%atom_cls_no(i,ib) = molecule%atom_cls_no(iatm)
      end do

      ! total mass and center of mass in the template (body-fixed) frame
      !
      total_mass  = 0.0_wp
      com_ref(1:3) = 0.0_wp

      do i = 1, n
        iatm = rigidbody%atomlist(i,ib)
        m = molecule%mass(iatm)
        total_mass   = total_mass   + m
        com_ref(1:3) = com_ref(1:3) + m * rb_ref%coord(1:3,i)
      end do

      if (total_mass <= 0.0_wp) &
        call error_msg('Setup_Rigidbody> total mass of a rigid body is zero')

      com_ref(1:3) = com_ref(1:3) / total_mass

      ! inertia tensor about the center of mass, in the template frame
      !
      inertia_mat(1:3,1:3) = 0.0_wp

      do i = 1, n
        iatm = rigidbody%atomlist(i,ib)
        m = molecule%mass(iatm)
        r(1:3) = rb_ref%coord(1:3,i) - com_ref(1:3)

        inertia_mat(1,1) = inertia_mat(1,1) + m*(r(2)*r(2) + r(3)*r(3))
        inertia_mat(2,2) = inertia_mat(2,2) + m*(r(1)*r(1) + r(3)*r(3))
        inertia_mat(3,3) = inertia_mat(3,3) + m*(r(1)*r(1) + r(2)*r(2))
        inertia_mat(1,2) = inertia_mat(1,2) - m*r(1)*r(2)
        inertia_mat(1,3) = inertia_mat(1,3) - m*r(1)*r(3)
        inertia_mat(2,3) = inertia_mat(2,3) - m*r(2)*r(3)
      end do

      inertia_mat(2,1) = inertia_mat(1,2)
      inertia_mat(3,1) = inertia_mat(1,3)
      inertia_mat(3,2) = inertia_mat(2,3)

      ! diagonalize to get principal moments (eval) and principal axes (evec)
      !
      evec(1:3,1:3) = inertia_mat(1:3,1:3)

#ifdef LAPACK
      call dsyev('V', 'U', 3, evec, 3, eval, work, 24, ierr)
      if (ierr /= 0) &
        call error_msg('Setup_Rigidbody> failed to diagonalize the '// &
                        'inertia tensor')
#else
      call error_msg('Setup_Rigidbody> ERROR: this feature needs LAPACK.')
#endif

      rigidbody%mass(ib)     = total_mass
      rigidbody%inv_mass(ib) = 1.0_wp / total_mass

      do i = 1, 3
        rigidbody%inertia(i,ib) = eval(i)
        if (eval(i) > EPS) then
          rigidbody%inv_inertia(i,ib) = 1.0_wp / eval(i)
        else
          rigidbody%inv_inertia(i,ib) = 0.0_wp
        end if
      end do

      ! express the reference geometry in the principal-axis frame:
      ! ref_coord_principal = evec^T * (template_coord - com_ref)
      !
      do i = 1, n
        r(1:3) = rb_ref%coord(1:3,i) - com_ref(1:3)
        rigidbody%ref_coord(1,i,ib) = evec(1,1)*r(1)+evec(2,1)*r(2)+evec(3,1)*r(3)
        rigidbody%ref_coord(2,i,ib) = evec(1,2)*r(1)+evec(2,2)*r(2)+evec(3,2)*r(3)
        rigidbody%ref_coord(3,i,ib) = evec(1,3)*r(1)+evec(2,3)*r(2)+evec(3,3)*r(3)
      end do

    end do

    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Rigidbody> Summary of rigid-body information'
      write(MsgOut,'(A20,I10)') '  num_bodies      = ', rigidbody%num_bodies
      write(MsgOut,'(A20,I10)') '  natom_per_body  = ', rigidbody%max_natom
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_rigidbody

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      rigidbody_diameter
  !> @brief        largest pairwise distance among a rigid body's reference
  !!               atoms (used by spdyn to check compatibility with the cell
  !!               size of the spatial domain decomposition)
  !! @param[in]    rigidbody : rigid-body information
  !! @param[in]    ibody     : body index
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function rigidbody_diameter(rigidbody, ibody) result(diameter)

    ! formal arguments
    type(s_rigidbody),        intent(in) :: rigidbody
    integer,                  intent(in) :: ibody

    ! return value
    real(wp)                  :: diameter

    ! local variables
    integer                    :: i, j, n
    real(wp)                   :: d(3), d2, maxd2


    n = rigidbody%natom(ibody)
    maxd2 = 0.0_wp

    do i = 1, n-1
      do j = i+1, n
        d(1:3) = rigidbody%ref_coord(1:3,i,ibody) - rigidbody%ref_coord(1:3,j,ibody)
        d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
        maxd2 = max(maxd2, d2)
      end do
    end do

    diameter = sqrt(maxd2)

    return

  end function rigidbody_diameter

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quat_to_rotmatrix
  !> @brief        convert a unit quaternion (w,x,y,z) to a rotation matrix
  !!               that maps body(principal)-frame vectors to the space frame
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine quat_to_rotmatrix(q, R)

    ! formal arguments
    real(wp),                 intent(in)  :: q(4)
    real(wp),                 intent(out) :: R(3,3)

    ! local variables
    real(wp)                  :: w, x, y, z


    w = q(1); x = q(2); y = q(3); z = q(4)

    R(1,1) = 1.0_wp - 2.0_wp*(y*y + z*z)
    R(1,2) = 2.0_wp*(x*y - w*z)
    R(1,3) = 2.0_wp*(x*z + w*y)
    R(2,1) = 2.0_wp*(x*y + w*z)
    R(2,2) = 1.0_wp - 2.0_wp*(x*x + z*z)
    R(2,3) = 2.0_wp*(y*z - w*x)
    R(3,1) = 2.0_wp*(x*z - w*y)
    R(3,2) = 2.0_wp*(y*z + w*x)
    R(3,3) = 1.0_wp - 2.0_wp*(x*x + y*y)

    return

  end subroutine quat_to_rotmatrix

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      quat_normalize
  !> @brief        normalize a quaternion (guards against zero-norm input)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function quat_normalize(q) result(qn)

    ! formal arguments
    real(wp),                 intent(in) :: q(4)

    ! return value
    real(wp)                  :: qn(4)

    ! local variables
    real(wp)                  :: norm


    norm = sqrt(q(1)*q(1) + q(2)*q(2) + q(3)*q(3) + q(4)*q(4))

    if (norm > EPS) then
      qn(1:4) = q(1:4) / norm
    else
      qn(1:4) = (/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
    end if

    return

  end function quat_normalize

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    quat_derivative
  !> @brief        time derivative of the orientation quaternion given the
  !!               space-frame angular velocity: dq/dt = 1/2 * omega_space ⊗ q
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine quat_derivative(q, omega_space, qdot)

    ! formal arguments
    real(wp),                 intent(in)  :: q(4)
    real(wp),                 intent(in)  :: omega_space(3)
    real(wp),                 intent(out) :: qdot(4)

    ! local variables
    real(wp)                  :: wx, wy, wz


    wx = omega_space(1); wy = omega_space(2); wz = omega_space(3)

    qdot(1) = -0.5_wp*(wx*q(2) + wy*q(3) + wz*q(4))
    qdot(2) =  0.5_wp*(wx*q(1) + wy*q(4) - wz*q(3))
    qdot(3) =  0.5_wp*(wy*q(1) + wz*q(2) - wx*q(4))
    qdot(4) =  0.5_wp*(wz*q(1) + wx*q(3) - wy*q(2))

    return

  end subroutine quat_derivative

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rigidbody_angvel
  !> @brief        compute body-frame and space-frame angular velocity from
  !!               the orientation quaternion and space-frame angular
  !!               momentum, following the standard rigid-body scheme
  !!               (e.g. Miller et al., J. Chem. Phys. 116, 8649 (2002)):
  !!               rotate L into the body (principal-axis) frame, apply the
  !!               diagonal inverse inertia there, rotate the result back.
  !! @param[in]    q             : orientation quaternion (w,x,y,z)
  !! @param[in]    angmom_space  : angular momentum, space frame
  !! @param[in]    inv_inertia   : diagonal inverse inertia, body frame (3)
  !! @param[out]   omega_space   : angular velocity, space frame
  !! @param[out]   R             : rotation matrix (body -> space) used
  !! @param[out]   omega_body    : angular velocity, body (principal) frame
  !!                               (optional; used for the exponential-map
  !!               orientation update, see quat_rotate_by_body_omega)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rigidbody_angvel(q, angmom_space, inv_inertia, omega_space, R, &
                              omega_body)

    ! formal arguments
    real(wp),                 intent(in)  :: q(4)
    real(wp),                 intent(in)  :: angmom_space(3)
    real(wp),                 intent(in)  :: inv_inertia(3)
    real(wp),                 intent(out) :: omega_space(3)
    real(wp),                 intent(out) :: R(3,3)
    real(wp), optional,       intent(out) :: omega_body(3)

    ! local variables
    real(wp)                  :: angmom_body(3), wbody(3)


    call quat_to_rotmatrix(q, R)

    ! angmom_body = R^T * angmom_space
    angmom_body(1) = R(1,1)*angmom_space(1) + R(2,1)*angmom_space(2) &
                   + R(3,1)*angmom_space(3)
    angmom_body(2) = R(1,2)*angmom_space(1) + R(2,2)*angmom_space(2) &
                   + R(3,2)*angmom_space(3)
    angmom_body(3) = R(1,3)*angmom_space(1) + R(2,3)*angmom_space(2) &
                   + R(3,3)*angmom_space(3)

    wbody(1:3) = inv_inertia(1:3) * angmom_body(1:3)
    if (present(omega_body)) omega_body(1:3) = wbody(1:3)

    ! omega_space = R * omega_body
    omega_space(1) = R(1,1)*wbody(1) + R(1,2)*wbody(2) + R(1,3)*wbody(3)
    omega_space(2) = R(2,1)*wbody(1) + R(2,2)*wbody(2) + R(2,3)*wbody(3)
    omega_space(3) = R(3,1)*wbody(1) + R(3,2)*wbody(2) + R(3,3)*wbody(3)

    return

  end subroutine rigidbody_angvel

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      quat_multiply
  !> @brief        Hamilton product q1 * q2 (w,x,y,z convention); composing
  !!               a space-frame orientation q1 with a body-frame rotation
  !!               increment q2 on the right gives the new orientation
  !!               (used by the exponential-map rotation update)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function quat_multiply(q1, q2) result(q)

    ! formal arguments
    real(wp),                 intent(in) :: q1(4), q2(4)

    ! return value
    real(wp)                  :: q(4)


    q(1) = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4)
    q(2) = q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3)
    q(3) = q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2)
    q(4) = q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(1)

    return

  end function quat_multiply

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      quat_rotate_by_body_omega
  !> @brief        advance an orientation quaternion by exactly rotating it
  !!               through angle |omega_body|*dt about the axis
  !!               omega_body/|omega_body| (the body-frame exponential map),
  !!               composed on the right of q -- a geometric integrator for
  !!               SO(3) with much better long-term energy behavior than a
  !!               linearized (forward-Euler) quaternion derivative step,
  !!               and, unlike that step, never needs renormalization
  !!               (kept anyway as a numerical-precision safety net)
  !! @param[in]    q          : current orientation quaternion
  !! @param[in]    omega_body : angular velocity, body (principal) frame
  !! @param[in]    dt         : time step
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function quat_rotate_by_body_omega(q, omega_body, dt) result(qn)

    ! formal arguments
    real(wp),                 intent(in) :: q(4)
    real(wp),                 intent(in) :: omega_body(3)
    real(wp),                 intent(in) :: dt

    ! return value
    real(wp)                  :: qn(4)

    ! local variables
    real(wp)                  :: wnorm, half_angle, dq(4)


    wnorm = sqrt(sum(omega_body(1:3)*omega_body(1:3)))

    if (wnorm > EPS) then
      half_angle = 0.5_wp * wnorm * dt
      dq(1)   = cos(half_angle)
      dq(2:4) = (sin(half_angle)/wnorm) * omega_body(1:3)
    else
      dq = (/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
    end if

    qn = quat_normalize(quat_multiply(q, dq))

    return

  end function quat_rotate_by_body_omega

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fit_rigidbody_quat
  !> @brief        find the orientation quaternion q such that
  !!               quat_to_rotmatrix(q) best superimposes refcoord
  !!               (body-fixed principal frame) onto movcoord (current
  !!               space-frame positions), by the mass-weighted quaternion
  !!               RMSD method (Kearsley 1989 / Kabsch), used once at
  !!               startup to initialize a rigid body's orientation from
  !!               its actual starting atom positions. Both coordinate
  !!               arrays are recentered on their own mass-weighted center
  !!               internally, so callers need not pre-center them.
  !! @param[in]    n        : number of atoms
  !! @param[in]    refcoord : body-fixed principal-frame coordinates (3,n)
  !! @param[in]    movcoord : current space-frame coordinates (3,n)
  !! @param[in]    mass     : per-atom mass (n)
  !! @param[out]   quat     : orientation quaternion (w,x,y,z), normalized
  !! @param[out]   ierr     : 0 on success
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fit_rigidbody_quat(n, refcoord, movcoord, mass, quat, ierr)

    ! parameters
    integer,                  parameter    :: lwork = max(1,4*3-1)

    ! formal arguments
    integer,                  intent(in)   :: n
    real(wp),                 intent(in)   :: refcoord(:,:), movcoord(:,:)
    real(wp),                 intent(in)   :: mass(:)
    real(wp),                 intent(out)  :: quat(4)
    integer,                  intent(out)  :: ierr

    ! local variables
    real(wp)                  :: total_mass, com_ref(3), com_mov(3)
    real(wp)                  :: dref(3), dmov(3), dadd(3), dsub(3)
    real(wp)                  :: sym_matrix(4,4), eval(4), work(lwork)
    integer                   :: i


    if (n == 0) then
      ierr = -1
      return
    end if

    total_mass  = 0.0_wp
    com_ref(1:3) = 0.0_wp
    com_mov(1:3) = 0.0_wp

    do i = 1, n
      total_mass   = total_mass   + mass(i)
      com_ref(1:3) = com_ref(1:3) + mass(i) * refcoord(1:3,i)
      com_mov(1:3) = com_mov(1:3) + mass(i) * movcoord(1:3,i)
    end do

    if (total_mass <= 0.0_wp) then
      ierr = -2
      return
    end if

    com_ref(1:3) = com_ref(1:3) / total_mass
    com_mov(1:3) = com_mov(1:3) / total_mass

    sym_matrix(1:4,1:4) = 0.0_wp

    do i = 1, n
      dref(1:3) = refcoord(1:3,i) - com_ref(1:3)
      dmov(1:3) = movcoord(1:3,i) - com_mov(1:3)
      dsub(1:3) = mass(i) * (dref(1:3) - dmov(1:3))
      dadd(1:3) = mass(i) * (dref(1:3) + dmov(1:3))

      sym_matrix(1,1) = sym_matrix(1,1) + dsub(1)*dsub(1) &
                                        + dsub(2)*dsub(2) &
                                        + dsub(3)*dsub(3)
      sym_matrix(1,2) = sym_matrix(1,2) + dadd(2)*dsub(3) - dsub(2)*dadd(3)
      sym_matrix(1,3) = sym_matrix(1,3) + dsub(1)*dadd(3) - dadd(1)*dsub(3)
      sym_matrix(1,4) = sym_matrix(1,4) + dadd(1)*dsub(2) - dsub(1)*dadd(2)
      sym_matrix(2,2) = sym_matrix(2,2) + dsub(1)*dsub(1) &
                                        + dadd(2)*dadd(2) &
                                        + dadd(3)*dadd(3)
      sym_matrix(2,3) = sym_matrix(2,3) + dsub(1)*dsub(2) - dadd(1)*dadd(2)
      sym_matrix(2,4) = sym_matrix(2,4) + dsub(1)*dsub(3) - dadd(1)*dadd(3)
      sym_matrix(3,3) = sym_matrix(3,3) + dadd(1)*dadd(1) &
                                        + dsub(2)*dsub(2) &
                                        + dadd(3)*dadd(3)
      sym_matrix(3,4) = sym_matrix(3,4) + dsub(2)*dsub(3) - dadd(2)*dadd(3)
      sym_matrix(4,4) = sym_matrix(4,4) + dadd(1)*dadd(1) &
                                        + dadd(2)*dadd(2) &
                                        + dsub(3)*dsub(3)
    end do

    sym_matrix(2,1) = sym_matrix(1,2)
    sym_matrix(3,1) = sym_matrix(1,3)
    sym_matrix(3,2) = sym_matrix(2,3)
    sym_matrix(4,1) = sym_matrix(1,4)
    sym_matrix(4,2) = sym_matrix(2,4)
    sym_matrix(4,3) = sym_matrix(3,4)

#ifdef LAPACK
    call dsyev('V', 'U', 4, sym_matrix, 4, eval, work, lwork, ierr)
#else
    call error_msg('Fit_Rigidbody_Quat> ERROR: this feature needs LAPACK.')
#endif

    if (ierr /= 0) then
      ierr = -3
      return
    end if

    ! smallest eigenvalue's eigenvector is the best-fit quaternion; with
    ! refcoord = body-fixed principal frame and movcoord = space frame,
    ! quat_to_rotmatrix(quat) maps body-fixed vectors to the space frame
    ! (see the derivation note above fit_rigidbody_quat's declaration)
    quat(1:4) = quat_normalize(sym_matrix(1:4,1))

    ierr = 0

    return

  end subroutine fit_rigidbody_quat

end module rigidbody_mod
