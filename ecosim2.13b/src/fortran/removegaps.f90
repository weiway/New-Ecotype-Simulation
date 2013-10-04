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
!> The removegaps program removes all nucleotides in an aligned set of
!> sequences that contain any character other than A,C,G, or T.
!>
!> @pre  Requires that the sequencesfasta.txt and numbers.dat files exist and
!>         that the contain the sequence data and other variables necessary
!>         to perform the removal of gaps.
!>
!> @post Generates the fasta.dat file that contains the sequences without
!>         gaps.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program removegaps
  use methods
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 20)             :: strainname_format
  character(len = 256)            :: sequencesfasta_txt
  character(len = 256)            :: numbers_dat
  character(len = 256)            :: fasta_dat
  character(len = 1), allocatable :: strainname(:,:)
  character(len = 1), allocatable :: sequence(:,:)
  logical                         :: badnuc
  integer                         :: allocate_status
  integer                         :: numstrain
  integer                         :: length
  integer                         :: i
  integer                         :: jnuc
  integer                         :: knuc
  integer                         :: mnuc
  integer                         :: jseq
  integer                         :: idigit
  integer, parameter              :: max_strainname = 20
  integer, parameter              :: sequences_unit = 1
  integer, parameter              :: numbers_unit = 2
  integer, parameter              :: output_unit = 3
  logical                         :: file_exists
  ! Provide default file names to use.
  sequencesfasta_txt = 'sequences.dat'
  numbers_dat = 'numbers.dat'
  fasta_dat = 'removegaps.dat'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, sequencesfasta_txt)
      case (2)
        call getArgument (2, numbers_dat)
      case (3)
        call getArgument (3, fasta_dat)
      case (4)
        call getArgument (4, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: sequences.dat numbers.dat removegaps.dat"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input files exist.
  inquire (file = trim (sequencesfasta_txt), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) &
      "The sequences.dat file was not found at: ", &
      trim (sequencesfasta_txt)
    ! Error, exit the program.
    stop
  end if
  inquire (file = trim (numbers_dat), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The numbers.dat file was not found at: ", &
      trim (numbers_dat)
    ! Error, exit the program.
    stop
  end if
  ! Open Input files
  open (unit = sequences_unit, file = trim (sequencesfasta_txt), &
    access = 'sequential', form = 'formatted')
  open (unit = numbers_unit, file = trim (numbers_dat), &
    access = 'sequential', form = 'formatted')
  ! Open Output file
  open (unit = output_unit, file = trim (fasta_dat), access = 'sequential', &
    form = 'formatted')
  ! Read in numbers_dat
  read (unit = numbers_unit, fmt = *) numstrain, length
  ! length is the number of nucleotides
  ! numstrain is the number of sequences
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
  ! Read in sequencesfasta_txt
  do jseq = 1, numstrain
    read (unit = sequences_unit, fmt = trim (strainname_format)) &
      (strainname(jseq, idigit), idigit = 1, max_strainname)
    read (unit = sequences_unit, fmt = trim (sequence_format)) &
      (sequence(jseq, jnuc), jnuc = 1, length)
    do mnuc = 1, length
      if (sequence(jseq, mnuc) .eq. 'A') sequence(jseq, mnuc) = 'a'
      if (sequence(jseq, mnuc) .eq. 'G') sequence(jseq, mnuc) = 'g'
      if (sequence(jseq, mnuc) .eq. 'T') sequence(jseq, mnuc) = 't'
      if (sequence(jseq, mnuc) .eq. 'C') sequence(jseq, mnuc) = 'c'
    end do
  end do
  ! ok, so we've read in all the sequences.  Now we get rid of
  ! all nucleotide sites in which even a single sequence has a '-'
  knuc = 0
  ! knuc is the number of nucleotides that have been found to be gap free
  do jnuc = 1, length
    badnuc = .false.
    do jseq = 1, numstrain
      if (sequence(jseq, jnuc) .ne. 'a' .and. &
          sequence(jseq, jnuc) .ne. 'c' .and. &
          sequence(jseq, jnuc) .ne. 'g' .and. &
          sequence(jseq, jnuc) .ne. 't') then
        badnuc = .true.
        exit
      end if
    end do
    if (badnuc) cycle
    knuc = knuc + 1
    do jseq = 1, numstrain
      sequence(jseq, knuc) = sequence(jseq, jnuc)
    end do
  end do
  ! Output to fasta_dat
  write (unit = output_unit, fmt = *) numstrain, knuc
  do jseq = 1, numstrain
    write (unit = output_unit, fmt = trim (strainname_format)) &
      (strainname(jseq, idigit), idigit = 1, max_strainname)
    write (unit = output_unit, fmt = trim (sequence_format)) &
      (sequence(jseq, jnuc), jnuc = 1, knuc)
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
  close (unit = sequences_unit)
  close (unit = numbers_unit)
  close (unit = output_unit)
  ! Successful termination of program.
  stop
end program removegaps
