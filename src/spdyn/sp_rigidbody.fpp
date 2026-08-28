!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_rigidbody_mod
!> @brief   [RIGIDBODY] control section and spdyn-side rigid-body setup
!! @authors Genesis Developers
!
!  (c) Copyright 2026 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_rigidbody_mod

  use rigidbody_str_mod
  use rigidbody_mod
  use fileio_rigidbody_mod
  use sp_boundary_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_rgbd_info
    character(MaxFilename)   :: indexfile = ''
    character(MaxFilename)   :: reffile   = ''
  end type s_rgbd_info

  ! subroutines
  public :: show_ctrl_rigidbody
  public :: read_ctrl_rigidbody
  public :: setup_rigidbody_spdyn

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_rigidbody
  !> @brief        show RIGIDBODY section usage
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine show_ctrl_rigidbody(show_all)

    ! formal arguments
    logical,                 intent(in) :: show_all


    if (show_all) then
      write(MsgOut,'(A)') '[RIGIDBODY]'
      write(MsgOut,'(A)') '# indexfile = index.dat  # rigid-body atom-index groups'
      write(MsgOut,'(A)') '# reffile   = ref.dat    # body-fixed reference coordinates'
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine show_ctrl_rigidbody

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_rigidbody
  !> @brief        read RIGIDBODY section in the control file
  !! @param[in]    handle    : handle for control file
  !! @param[out]   rgbd_info : RIGIDBODY section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_rigidbody(handle, rgbd_info)

    ! parameters
    character(*),            parameter     :: Section = 'Rigidbody'

    ! formal arguments
    integer,                  intent(in)    :: handle
    type(s_rgbd_info),        intent(inout) :: rgbd_info


    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_string(handle, Section, 'indexfile', rgbd_info%indexfile)
    call read_ctrlfile_string(handle, Section, 'reffile',   rgbd_info%reffile)

    call end_ctrlfile_section(handle)

    if (main_rank .and. rgbd_info%indexfile /= '') then
      write(MsgOut,'(A)') 'Read_Ctrl_Rigidbody> Parameters of Rigid-Body Dynamics'
      write(MsgOut,'(A20,A)') '  indexfile       = ', trim(rgbd_info%indexfile)
      write(MsgOut,'(A20,A)') '  reffile         = ', trim(rgbd_info%reffile)
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_rigidbody

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_rigidbody_spdyn
  !> @brief        build rigid-body information for spdyn and check that it
  !!               is compatible with the spatial domain decomposition (the
  !!               [RIGIDBODY] feature is a no-op when indexfile/reffile are
  !!               not given)
  !! @param[in]    rgbd_info : RIGIDBODY section control parameters information
  !! @param[in]    molecule  : molecule information (global, pre-decomposition)
  !! @param[in]    boundary  : boundary condition information
  !! @param[out]   rigidbody : rigid-body information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_rigidbody_spdyn(rgbd_info, molecule, boundary, rigidbody)

    ! formal arguments
    type(s_rgbd_info),        intent(in)    :: rgbd_info
    type(s_molecule),         intent(in)    :: molecule
    type(s_boundary),         intent(in)    :: boundary
    type(s_rigidbody),        intent(inout) :: rigidbody

    ! local variables
    type(s_rigidbody_index)   :: rb_index
    type(s_rigidbody_ref)     :: rb_ref
    real(wp)                  :: diameter, cell_edge
    integer                   :: ib


    call init_rigidbody(rigidbody)

    if (rgbd_info%indexfile == '' .or. rgbd_info%reffile == '') return

    call input_rigidbody_index   (rgbd_info%indexfile, rb_index)
    call input_rigidbody_refcoord(rgbd_info%reffile,   rb_ref)

    call setup_rigidbody(molecule, rb_index, rb_ref, rigidbody)

    ! spdyn's spatial domain decomposition places every atom of a rigid
    ! body into a single cell, chosen from one representative atom (see
    ! setup_atom_by_rigidbody in sp_domain.fpp), and the nonbonded pairlist
    ! only searches neighboring cells; a rigid body whose diameter reaches
    ! or exceeds the cell edge could then have real, nearby interactions
    ! fall outside that search and be silently missed. Rather than widen
    ! the pairlist search (which would also touch SETTLE/SHAKE code paths),
    ! this is rejected at setup time -- see tests/rigid-body/SPEC.md section 2.
    !
    cell_edge = real(min(boundary%cell_size_x, boundary%cell_size_y, &
                         boundary%cell_size_z), wp)

    do ib = 1, rigidbody%num_bodies
      diameter = rigidbody_diameter(rigidbody, ib)
      if (diameter >= cell_edge) then
        write(MsgOut,'(A,I0,A,F10.4,A,F10.4)')                          &
          'Setup_Rigidbody_Spdyn> body ', ib, ' diameter = ', diameter, &
          ' >= domain cell edge = ', cell_edge
        call error_msg('Setup_Rigidbody_Spdyn> rigid body is too large '// &
                        'for the current spatial domain decomposition')
      end if
    end do

    call dealloc_rigidbody_index(rb_index)
    call dealloc_rigidbody_ref(rb_ref)

    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Rigidbody_Spdyn> Summary of rigid-body setup'
      write(MsgOut,'(A20,I10)')  '  num_bodies      = ', rigidbody%num_bodies
      write(MsgOut,'(A20,F10.4)')'  min cell edge   = ', cell_edge
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine setup_rigidbody_spdyn

end module sp_rigidbody_mod
