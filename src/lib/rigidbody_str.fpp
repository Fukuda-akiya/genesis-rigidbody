!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rigidbody_str_mod
!> @brief   structure of rigid-body information
!! @authors Genesis Developers
!
!  (c) Copyright 2026 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module rigidbody_str_mod

  use messages_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_rigidbody

    logical                       :: is_used    = .false.
    integer                       :: num_bodies = 0
    integer                       :: max_natom  = 0

    ! membership: global atom index of each member atom, 0-padded up to max_natom
    integer,          allocatable :: natom(:)          ! (num_bodies)
    integer,          allocatable :: atomlist(:,:)     ! (max_natom, num_bodies)

    ! per-atom static properties, indexed the same way as atomlist; kept
    ! here (rather than looked up from s_molecule at runtime) because
    ! s_molecule is deallocated after spdyn setup completes, while this
    ! structure is replicated on every rank and kept for the whole run
    real(wp),         allocatable :: atom_mass(:,:)    ! (max_natom, num_bodies)
    real(wp),         allocatable :: atom_charge(:,:)  ! (max_natom, num_bodies)
    integer,          allocatable :: atom_cls_no(:,:)  ! (max_natom, num_bodies)

    ! body-fixed reference coordinates, expressed in each body's own principal
    ! axis frame and centered on that body's own center of mass
    real(wp),         allocatable :: ref_coord(:,:,:)  ! (3, max_natom, num_bodies)

    ! mass / inertia (principal moments), computed once at setup
    real(wp),         allocatable :: mass(:)           ! (num_bodies)
    real(wp),         allocatable :: inv_mass(:)       ! (num_bodies)
    real(wp),         allocatable :: inertia(:,:)      ! (3, num_bodies)
    real(wp),         allocatable :: inv_inertia(:,:)  ! (3, num_bodies)

  end type s_rigidbody

  ! parameters for allocatable variables
  integer,      public, parameter :: RigidBodyList  = 1
  integer,      public, parameter :: RigidBodyPhys  = 2

  ! subroutines
  public :: init_rigidbody
  public :: alloc_rigidbody
  public :: dealloc_rigidbody
  public :: dealloc_rigidbody_all

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_rigidbody
  !> @brief        initialize rigid-body information
  !! @param[out]   rigidbody : rigid-body information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_rigidbody(rigidbody)

    ! formal arguments
    type(s_rigidbody),        intent(inout) :: rigidbody


    rigidbody%is_used    = .false.
    rigidbody%num_bodies = 0
    rigidbody%max_natom  = 0

    return

  end subroutine init_rigidbody

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_rigidbody
  !> @brief        allocate rigid-body information
  !! @param[inout] rigidbody : rigid-body information
  !! @param[in]    variable  : selected variable
  !! @param[in]    var_size  : size of the selected variable (num_bodies)
  !! @param[in]    var_size1 : 2nd size of the selected variable (max_natom)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_rigidbody(rigidbody, variable, var_size, var_size1)

    ! formal arguments
    type(s_rigidbody),        intent(inout) :: rigidbody
    integer,                  intent(in)    :: variable
    integer,                  intent(in)    :: var_size
    integer, optional,        intent(in)    :: var_size1

    ! local variables
    integer                   :: alloc_stat
    integer                   :: dealloc_stat


    alloc_stat   = 0
    dealloc_stat = 0

    select case (variable)

    case (RigidBodyList)

      if (allocated(rigidbody%atomlist)) then
        if (size(rigidbody%atomlist(:,1)) /= var_size1 .or. &
            size(rigidbody%natom)         /= var_size) then
          deallocate(rigidbody%natom,       &
                     rigidbody%atomlist,    &
                     rigidbody%atom_mass,   &
                     rigidbody%atom_charge, &
                     rigidbody%atom_cls_no, &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(rigidbody%atomlist)) &
        allocate(rigidbody%natom(var_size),               &
                 rigidbody%atomlist(var_size1, var_size),  &
                 rigidbody%atom_mass(var_size1, var_size),   &
                 rigidbody%atom_charge(var_size1, var_size), &
                 rigidbody%atom_cls_no(var_size1, var_size), &
                 stat = alloc_stat)

      rigidbody%natom      (1:var_size)              = 0
      rigidbody%atomlist   (1:var_size1, 1:var_size) = 0
      rigidbody%atom_mass  (1:var_size1, 1:var_size) = 0.0_wp
      rigidbody%atom_charge(1:var_size1, 1:var_size) = 0.0_wp
      rigidbody%atom_cls_no(1:var_size1, 1:var_size) = 0

    case (RigidBodyPhys)

      if (allocated(rigidbody%mass)) then
        if (size(rigidbody%mass) /= var_size .or. &
            size(rigidbody%ref_coord(1,:,1)) /= var_size1) then
          deallocate(rigidbody%ref_coord,   &
                     rigidbody%mass,        &
                     rigidbody%inv_mass,    &
                     rigidbody%inertia,     &
                     rigidbody%inv_inertia, &
                     stat = dealloc_stat)
        end if
      end if

      if (.not. allocated(rigidbody%mass)) &
        allocate(rigidbody%ref_coord  (3, var_size1, var_size), &
                 rigidbody%mass       (var_size),               &
                 rigidbody%inv_mass   (var_size),                &
                 rigidbody%inertia    (3, var_size),             &
                 rigidbody%inv_inertia(3, var_size),             &
                 stat = alloc_stat)

      rigidbody%ref_coord  (1:3, 1:var_size1, 1:var_size) = 0.0_wp
      rigidbody%mass       (1:var_size)                   = 0.0_wp
      rigidbody%inv_mass   (1:var_size)                   = 0.0_wp
      rigidbody%inertia    (1:3, 1:var_size)              = 0.0_wp
      rigidbody%inv_inertia(1:3, 1:var_size)              = 0.0_wp

    end select

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_rigidbody

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rigidbody
  !> @brief        deallocate rigid-body information
  !! @param[inout] rigidbody : rigid-body information
  !! @param[in]    variable  : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rigidbody(rigidbody, variable)

    ! formal arguments
    type(s_rigidbody),        intent(inout) :: rigidbody
    integer,                  intent(in)    :: variable

    ! local variables
    integer                   :: dealloc_stat


    dealloc_stat = 0

    select case (variable)

    case (RigidBodyList)

      if (allocated(rigidbody%atomlist)) then
        deallocate(rigidbody%natom,       &
                   rigidbody%atomlist,    &
                   rigidbody%atom_mass,   &
                   rigidbody%atom_charge, &
                   rigidbody%atom_cls_no, &
                   stat = dealloc_stat)
      end if

    case (RigidBodyPhys)

      if (allocated(rigidbody%mass)) then
        deallocate(rigidbody%ref_coord,   &
                   rigidbody%mass,        &
                   rigidbody%inv_mass,    &
                   rigidbody%inertia,     &
                   rigidbody%inv_inertia, &
                   stat = dealloc_stat)
      end if

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_rigidbody

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rigidbody_all
  !> @brief        deallocate all rigid-body information
  !! @param[inout] rigidbody : rigid-body information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rigidbody_all(rigidbody)

    ! formal arguments
    type(s_rigidbody),        intent(inout) :: rigidbody


    call dealloc_rigidbody(rigidbody, RigidBodyList)
    call dealloc_rigidbody(rigidbody, RigidBodyPhys)

    return

  end subroutine dealloc_rigidbody_all

end module rigidbody_str_mod
