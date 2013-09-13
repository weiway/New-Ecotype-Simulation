!    Ecotype Simulation models the sequence diversity within a bacterial clade as
!    the evolutionary result of net ecotype formation, periodic selection,
!    and drift, yielding a certain number of ecotypes.
! 
!    Copyright (C) 2009  Fred Cohan, Wesleyan University
!                        Carlo Francisco, Wesleyan University
!                        Danny Krizanc, Wesleyan University
!                        Andrew Warner, Wesleyan University
!                        Jason Wood, Montana State University
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program readsynec
  use Methods
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 20)             :: strainname_format
  character(len = 256)            :: fasta_dat
  character(len = 256)            :: population_dat
  character(len = 256)            :: namesofstrains_dat
  character(len = 1), allocatable :: sequence(:,:)
  character(len = 1), allocatable :: strainname(:,:)
  integer                         :: allocate_status
  integer                         :: allele
  integer                         :: alleledigit
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
  ! Read command line arguments.
  if (iargc() .ge. 3) then
    ! Expect filenames for fasta.dat, population.dat, and namesofstrains.dat.
    call getArgument (1, fasta_dat)
    call getArgument (2, population_dat)
    call getArgument (3, namesofstrains_dat)
  else if (iargc() .eq. 0) then
    ! No filenames supplied, so use defaults.
    write (unit = *, fmt = *) "Using default file names."
    fasta_dat = 'fasta.dat'
    population_dat = 'population.dat'
    namesofstrains_dat = 'namesofstrains.dat'
  else
    ! An unexpected number of arguments was supplied.
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: fasta.dat population.dat namesofstrains.dat"
    stop
  end if
  ! Check for optional debug command line argument.
  if (iargc() .eq. 4) then
    call getArgument (4, debug)
  else
    debug = .false.
  end if
  ! Verify that input file exist.
  inquire (file = trim (fasta_dat), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The fasta.dat file was not found at: ", trim (fasta_dat)
    stop
  end if
  ! Open Input file
  open (unit = fasta_unit, file = trim (fasta_dat), access = 'sequential', form = 'formatted')
  ! Open Output files
  open (unit = population_unit, file = trim (population_dat), access = 'sequential', form = 'formatted')
  open (unit = nameofstrains_unit, file = trim (namesofstrains_dat), access = 'sequential', form = 'formatted')
  ! Read in fasta_dat and output to population_dat and nameofstrains_dat
  read (unit = fasta_unit, fmt = *) numstrain, length
  write (unit = population_unit, fmt = *) numstrain, length
  write (unit = nameofstrains_unit, fmt = *) numstrain, length
  ! Allocate Variables.
  allocate (sequence(numstrain, length), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (strainname(numstrain, max_strainname), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  ! Create the format for printing sequence data.
  write (unit = strainname_format, fmt = "(a1,i0,a3)") "(", max_strainname, "a1)"
  write (unit = sequence_format, fmt = "(a1,i0,a3)") "(", length, "a1)"
  do jallele = 1, numstrain
    read (unit = fasta_unit, fmt = trim(strainname_format)) (strainname(jallele, idigit), idigit = 1, max_strainname)
    write (unit = nameofstrains_unit, fmt = trim (strainname_format)) (strainname(jallele, idigit), idigit = 2, max_strainname)
    read (unit = fasta_unit, fmt = trim(sequence_format)) (sequence(jallele, knuc), knuc = 1, length)
    do knuc = 1, length
      if (sequence(jallele, knuc) .eq. 'A') sequence(jallele, knuc) = 'a'
      if (sequence(jallele, knuc) .eq. 'G') sequence(jallele, knuc) = 'g'
      if (sequence(jallele, knuc) .eq. 'T') sequence(jallele, knuc) = 't'
      if (sequence(jallele, knuc) .eq. 'C') sequence(jallele, knuc) = 'c'
    end do
    write (unit = population_unit, fmt = trim(sequence_format)) (sequence(jallele, knuc), knuc = 1, length)
  end do
  close (unit = fasta_unit)
  close (unit = population_unit)
  close (unit = nameofstrains_unit)
  stop
end program readsynec
