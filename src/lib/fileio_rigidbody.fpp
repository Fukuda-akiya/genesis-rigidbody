!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_rigidbody_mod
!> @brief   read atom-index groups and body-fixed reference coordinates
!!          for rigid-body dynamics
!! @authors Genesis Developers
!
!  (c) Copyright 2026 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_rigidbody_mod

  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  !
  ! s_rbidx_group : atom indices (global, 1-based, referring to the same
  !                 numbering as PDB/PSF/CRD) belonging to one rigid-body copy
  type, public :: s_rbidx_group
    integer,        allocatable :: idx(:)
  end type s_rbidx_group

  ! s_rigidbody_index : one entry per rigid-body copy, as read from an
  !                     index file (repeated blocks of: count, blank line,
  !                     that many atom-index lines, blank line, ... until EOF)
  type, public :: s_rigidbody_index
    integer                                :: num_groups = 0
    type(s_rbidx_group),      allocatable  :: group(:)
  end type s_rigidbody_index

  ! s_rigidbody_ref : body-fixed frame reference coordinates of a single
  !                   representative rigid body, as read from a reference
  !                   coordinate file. atom_no(:) is informational only
  !                   (the i-th coordinate maps positionally onto the i-th
  !                   index of every group in s_rigidbody_index)
  type, public :: s_rigidbody_ref
    integer                                :: num_atoms = 0
    integer,                  allocatable  :: atom_no(:)
    real(wp),                 allocatable  :: coord(:,:)
  end type s_rigidbody_ref

  ! subroutines
  public  :: input_rigidbody_index
  public  :: input_rigidbody_refcoord
  public  :: dealloc_rigidbody_index
  public  :: dealloc_rigidbody_ref

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rigidbody_index
  !> @brief        open, read, and close a rigid-body atom-index file
  !! @param[in]    filename  : filename of the rigid-body index file
  !! @param[out]   rb_index  : rigid-body index information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rigidbody_index(filename, rb_index)

    ! formal arguments
    character(*),               intent(in)    :: filename
    type(s_rigidbody_index),    intent(inout) :: rb_index

    ! local variables
    integer                     :: unit_no


    call open_file(unit_no, filename, IOFileInput)
    call read_rigidbody_index(unit_no, rb_index)
    call close_file(unit_no)

    if (main_rank) then
      write(MsgOut,'(A)') 'Input_Rigidbody_Index> Summary of rigidbody index file'
      write(MsgOut,'(A20,I10)') '  num_groups      = ', rb_index%num_groups
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine input_rigidbody_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rigidbody_index
  !> @brief        read data from a rigid-body index file
  !! @param[in]    unit_no  : unit number of the rigid-body index file
  !! @param[out]   rb_index : rigid-body index information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rigidbody_index(unit_no, rb_index)

    ! formal arguments
    integer,                    intent(in)    :: unit_no
    type(s_rigidbody_index),    intent(inout) :: rb_index

    ! local variables
    integer                     :: i, ig, n, ios, ngroup, dummy


    ! first pass: count the number of blocks (rigid-body copies)
    !
    ! Note: each index value is read into a real variable (not skipped with
    ! a bare "read(unit,*)") so that list-directed input's blank-record
    ! skipping applies consistently with the second pass below -- the file
    ! format has a blank line both after the count and after the last index
    ! of each block, and a bare read only advances exactly one record.
    !
    ngroup = 0
    do
      read(unit_no, *, iostat=ios) n
      if (ios /= 0) exit
      if (n <= 0) call error_msg('Read_Rigidbody_Index> invalid group size')
      ngroup = ngroup + 1
      do i = 1, n
        read(unit_no, *, iostat=ios) dummy
        if (ios /= 0) &
          call error_msg('Read_Rigidbody_Index> unexpected end of file')
      end do
    end do

    if (ngroup == 0) &
      call error_msg('Read_Rigidbody_Index> no rigid-body group found')

    rewind(unit_no)

    rb_index%num_groups = ngroup
    allocate(rb_index%group(ngroup))

    ! second pass: read atom indices of each block
    !
    do ig = 1, ngroup
      read(unit_no, *) n
      allocate(rb_index%group(ig)%idx(n))
      do i = 1, n
        read(unit_no, *) rb_index%group(ig)%idx(i)
      end do
    end do

    return

  end subroutine read_rigidbody_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rigidbody_refcoord
  !> @brief        open, read, and close a rigid-body reference coordinate file
  !! @param[in]    filename  : filename of the rigid-body reference coordinate file
  !! @param[out]   rb_ref    : rigid-body reference coordinate information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rigidbody_refcoord(filename, rb_ref)

    ! formal arguments
    character(*),               intent(in)    :: filename
    type(s_rigidbody_ref),      intent(inout) :: rb_ref

    ! local variables
    integer                     :: unit_no


    call open_file(unit_no, filename, IOFileInput)
    call read_rigidbody_refcoord(unit_no, rb_ref)
    call close_file(unit_no)

    if (main_rank) then
      write(MsgOut,'(A)') 'Input_Rigidbody_Refcoord> Summary of rigidbody reference coordinate file'
      write(MsgOut,'(A20,I10)') '  num_atoms       = ', rb_ref%num_atoms
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine input_rigidbody_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_rigidbody_refcoord
  !> @brief        read data from a rigid-body reference coordinate file
  !! @param[in]    unit_no : unit number of the rigid-body reference coordinate file
  !! @param[out]   rb_ref  : rigid-body reference coordinate information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_rigidbody_refcoord(unit_no, rb_ref)

    ! formal arguments
    integer,                    intent(in)    :: unit_no
    type(s_rigidbody_ref),      intent(inout) :: rb_ref

    ! local variables
    integer                     :: i, n


    read(unit_no, *)          ! 'The number of mass points <n>'
    read(unit_no, *) n        ! number of atoms
    read(unit_no, *)          ! 'Coordinates in the body-fixed frame <n>'

    if (n <= 0) &
      call error_msg('Read_Rigidbody_Refcoord> invalid number of atoms')

    rb_ref%num_atoms = n
    allocate(rb_ref%atom_no(n), rb_ref%coord(3,n))

    do i = 1, n
      read(unit_no, *) rb_ref%atom_no(i), rb_ref%coord(1:3,i)
    end do

    return

  end subroutine read_rigidbody_refcoord

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rigidbody_index
  !> @brief        deallocate rigid-body index information
  !! @param[inout] rb_index : rigid-body index information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rigidbody_index(rb_index)

    ! formal arguments
    type(s_rigidbody_index),    intent(inout) :: rb_index

    ! local variables
    integer                     :: ig


    if (allocated(rb_index%group)) then
      do ig = 1, size(rb_index%group)
        if (allocated(rb_index%group(ig)%idx)) deallocate(rb_index%group(ig)%idx)
      end do
      deallocate(rb_index%group)
    end if

    return

  end subroutine dealloc_rigidbody_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_rigidbody_ref
  !> @brief        deallocate rigid-body reference coordinate information
  !! @param[inout] rb_ref : rigid-body reference coordinate information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_rigidbody_ref(rb_ref)

    ! formal arguments
    type(s_rigidbody_ref),      intent(inout) :: rb_ref


    if (allocated(rb_ref%atom_no)) deallocate(rb_ref%atom_no)
    if (allocated(rb_ref%coord))   deallocate(rb_ref%coord)

    return

  end subroutine dealloc_rigidbody_ref

end module fileio_rigidbody_mod
