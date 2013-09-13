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
!> The readsynec program extracts sequence data and the name of each strain
!> into separate files, making sure that the sequence data is in lower case.
!>
!> @pre  Requires that the fasta.dat file exists and that it contains
!>         sequence data.
!>
!> @post Generates the population.dat and nameofstrains.dat files that contain
!>         just the sequence data, and the name of each strain respectively.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program readsynec
  use methods
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 20)             :: strainname_format
  character(len = 256)            :: fasta_dat
  character(len = 256)            :: population_dat
  character(len = 256)            :: namesofstrains_dat
  character(len = 1), allocatable :: sequence(:,:)
  character(len = 1), allocatable :: strainname(:,:)
  integer                         :: allocate_status
  integer                         :: i
  integer                         :: jallele
  integer                         :: idigit
  integer                         :: knuc
  integer                         :: numstrain
  integer                         :: length
  integer, parameter              :: max_strainname = 20
  integer, parameter              :: fasta_unit = 1
  integer, parameter              :: population_unit = 2
  integer, parameter              :: nameofstrains_unit = 3
  logical                         :: file_exists
  ! Provide default file names to use.
  fasta_dat = 'removegaps.dat'
  population_dat = 'population.dat'
  namesofstrains_dat = 'namesofstrains.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, fasta_dat)
      case (2)
        call getArgument (2, population_dat)
      case (3)
        call getArgument (3, namesofstrains_dat)
      case (4)
        call getArgument (4, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: fasta.dat population.dat namesofstrains.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input file exist.
  inquire (file = trim (fasta_dat), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The fasta.dat file was not found at: ", &
      trim (fasta_dat)
    ! Error, exit the program.
    stop
  end if
  ! Open Input file
  open (unit = fasta_unit, file = trim (fasta_dat), access = 'sequential', &
    form = 'formatted')
  ! Open Output files
  open (unit = population_unit, file = trim (population_dat), &
    access = 'sequential', form = 'formatted')
  open (unit = nameofstrains_unit, file = trim (namesofstrains_dat), &
    access = 'sequential', form = 'formatted')
  ! Read in fasta_dat and output to population_dat and nameofstrains_dat
  read (unit = fasta_unit, fmt = *) numstrain, length
  write (unit = population_unit, fmt = *) numstrain, length
  write (unit = nameofstrains_unit, fmt = *) numstrain, length
  ! Allocate Variables.
  allocate (sequence(numstrain, length), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  allocate (strainname(numstrain, max_strainname), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: strainname!"
    ! Error, exit the program.
    stop
  end if
  ! Create the format for printing sequence data.
  write (unit = strainname_format, fmt = "(a1,i0,a3)") "(", max_strainname, &
    "a1)"
  write (unit = sequence_format, fmt = "(a1,i0,a3)") "(", length, "a1)"
  do jallele = 1, numstrain
    read (unit = fasta_unit, fmt = trim (strainname_format)) &
      (strainname(jallele, idigit), idigit = 1, max_strainname)
    write (unit = nameofstrains_unit, fmt = trim (strainname_format)) &
      (strainname(jallele, idigit), idigit = 2, max_strainname)
    read (unit = fasta_unit, fmt = trim (sequence_format)) &
      (sequence(jallele, knuc), knuc = 1, length)
    do knuc = 1, length
      if (sequence(jallele, knuc) .eq. 'A') sequence(jallele, knuc) = 'a'
      if (sequence(jallele, knuc) .eq. 'G') sequence(jallele, knuc) = 'g'
      if (sequence(jallele, knuc) .eq. 'T') sequence(jallele, knuc) = 't'
      if (sequence(jallele, knuc) .eq. 'C') sequence(jallele, knuc) = 'c'
    end do
    write (unit = population_unit, fmt = trim (sequence_format)) &
      (sequence(jallele, knuc), knuc = 1, length)
  end do
  ! Deallocate memory.
  deallocate (sequence, stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: sequence!"
    ! Error, exit the program.
    stop
  end if
  deallocate (strainname, stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: strainname!"
    ! Error, exit the program.
    stop
  end if
  ! Close data files.
  close (unit = fasta_unit)
  close (unit = population_unit)
  close (unit = nameofstrains_unit)
  ! Successful termination of program.
  stop
end program readsynec
