!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Ecotype Simulation models the sequence diversity within a bacterial clade
!    as the evolutionary result of net ecotype formation and periodic
!    selection, yielding a certain number of ecotypes.
!
!    Copyright (C) 2009-2013  Fred Cohan, Wesleyan University
!                             Danny Krizanc, Wesleyan University
!                             Jason M. Wood, Montana State University
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> The divergencematrix program generates a divergence matrix of the provided
!> sequence data.
!>
!> @pre  Requires that the correctpcr.out file exists and that it contains the
!>         sequence data necessary to generate the divergence matrix.
!> @post Generates the divergencematrix.dat file that contains the divergence
!>         matrix.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program divergencematrix
  use methods
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 20)             :: output_format
  character(len = 256)            :: correctpcr_out
  character(len = 256)            :: divergencematrix_dat
  character(len = 1), allocatable :: sequence(:,:)
  integer                         :: allocate_status
  integer                         :: numstrain
  integer                         :: length
  integer                         :: i
  integer                         :: istrain
  integer                         :: jstrain
  integer                         :: knuc
  integer                         :: kstrain
  integer, parameter              :: correctpcr_unit = 1
  integer, parameter              :: output_unit = 2
  logical                         :: file_exists
  real, allocatable               :: divergematrix(:,:)
  real, allocatable               :: identitymatrix(:,:)
  ! Provide default file names to use.
  correctpcr_out = 'correctpcr.out'
  divergencematrix_dat = 'divergencematrix.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, correctpcr_out)
      case (2)
        call getArgument (2, divergencematrix_dat)
      case (3)
        call getArgument (3, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: correctpcr.out divergencematrix.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input file exist.
  inquire (file = trim (correctpcr_out), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The correctpcr.out file was not found at: ", &
      trim (correctpcr_out)
    ! Error, exit the program.
    stop
  end if
  ! Open Input file
  open (unit = correctpcr_unit, file = trim (correctpcr_out), &
    access = 'sequential', form = 'formatted')
  ! Open Output file
  open (unit = output_unit, file = trim (divergencematrix_dat), &
    access = 'sequential', form = 'formatted')
  ! Read in correctpcr_out
  read (unit = correctpcr_unit, fmt = *) numstrain, length
  ! numstrain is the number of strains in the sample;
  ! length is the number of nucleotide sites
  ! Allocate variables.
  allocate (sequence(numstrain, length), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  allocate (divergematrix(numstrain, numstrain), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: divergematrix!"
    ! Error, exit the program.
    stop
  end if
  allocate (identitymatrix(numstrain, numstrain), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: identitymatrix!"
    ! Error, exit the program.
    stop
  end if
  ! Create the format for printing sequence data.
  write (unit = sequence_format, fmt = "(a1,i0,a3)") "(", length, "a1)"
  write (unit = output_format, fmt = "(a1,i0,a10)") "(", numstrain, &
    "(1x,f7.4))"
  ! Read in the sequence data
  do jstrain = 1, numstrain
    read (unit = correctpcr_unit, fmt = trim (sequence_format)) &
      (sequence(jstrain, knuc), knuc = 1, length)
  end do
  call diverge (numstrain, length, sequence, divergematrix)
  do jstrain = 1, numstrain
    do istrain = 1, numstrain
      identitymatrix(jstrain, istrain) = 1.0 - divergematrix(jstrain, istrain)
    end do
  end do
  ! Output to divergencematrix_dat
  write (unit = output_unit, fmt = *) numstrain, length
  do jstrain = 1, numstrain
    write (unit = output_unit, fmt = trim (output_format)) &
      (identitymatrix(jstrain, kstrain), kstrain = 1, numstrain)
  end do
  ! Deallocate memory.
  deallocate (sequence, stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  deallocate (divergematrix, stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) &
      "Failed to deallocate memory for: divergematrix!"
    ! Error, exit the program.
    stop
  end if
  deallocate (identitymatrix, stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) &
      "Failed to deallocate memory for: identitymatrix!"
    ! Error, exit the program.
    stop
  end if
  ! Close data files.
  close (unit = correctpcr_unit)
  close (unit = output_unit)
  ! Successful termination of program.
  stop
end program divergencematrix
