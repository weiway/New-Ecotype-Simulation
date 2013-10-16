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
!> The correctpcr program corrects for PCR error by assuming that each error
!> would be a unique nucleotide at one site.  Of all the sites with a unique
!> nucleotide, we randomly choose some of them to be corrected.  Each strain
!> chosen to be corrected has the offending nucleotide modified to match that
!> of the strain most closely related.
!>
!> @pre  Requires that the population.dat and pcrerror.dat files exist and
!>         that they contain all of the variables necessary to perform the
!>         correction.
!>
!> @post Generates the correctpcr.out file that contains the corrected
!>         sequences.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program correctpcr
  use methods
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 1), allocatable :: sequence(:,:)
  character(len = 1), parameter   :: code(4) = (/'a', 't', 'g', 'c'/)
  character(len = 256)            :: population_dat
  character(len = 256)            :: pcrerror_dat
  character(len = 256)            :: correctpcr_out
  integer                         :: allocate_status
  integer                         :: numnuc(4)
  integer                         :: i
  integer                         :: iii
  integer                         :: num
  integer                         :: numtochange
  integer                         :: closest
  integer                         :: seqchange
  integer                         :: nucchange
  integer                         :: numsingletons
  integer                         :: iother
  integer                         :: jstrain
  integer                         :: knuc
  integer                         :: nucnuc
  integer                         :: jchange
  integer                         :: numstrain
  integer                         :: length
  integer, allocatable            :: singleton(:,:)
  integer, allocatable            :: tochange(:,:)
  integer, parameter              :: population_unit = 1
  integer, parameter              :: output_unit = 2
  integer, parameter              :: pcrerror_unit = 4
  logical                         :: file_exists
  real                            :: x
  real                            :: diff
  real                            :: xmindiff
  real                            :: pcrerror
  real, allocatable               :: divergematrix(:,:)
  ! Provide default file names to use.
  population_dat = 'population.dat'
  pcrerror_dat = 'pcrerror.dat'
  correctpcr_out = 'correctpcr.out'
  ! Read command line arguments.
  do i = 1, iargc()
    select case (i)
      case (1)
        call getArgument (1, population_dat)
      case (2)
        call getArgument (2, pcrerror_dat)
      case (3)
        call getArgument (3, correctpcr_out)
      case (4)
        call getArgument (4, debug)
      case default
        ! An unexpected number of arguments was supplied.
        write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
        write (unit = *, fmt = *) &
          "Expected: population.dat pcrerror.dat correctpcr.out"
        ! Error, exit the program.
        stop
    end select
  end do
  ! Verify that input files exist.
  inquire (file = trim (population_dat), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The population.dat file was not found at: ", &
      trim (population_dat)
    ! Error, exit the program.
    stop
  end if
  inquire (file = trim (pcrerror_dat), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The pcrerror.dat file was not found at: ", &
      trim (pcrerror_dat)
    ! Error, exit the program.
    stop
  end if
  ! Open input files
  open (unit = population_unit, file = trim (population_dat), &
    access = 'sequential', form = 'formatted')
  open (unit = pcrerror_unit, file = trim (pcrerror_dat), &
    access = 'sequential', form = 'formatted')
  ! Open output file
  open (unit = output_unit, file = trim (correctpcr_out), &
    access = 'sequential', form = 'formatted')
  ! Read in population_dat
  read (unit = population_unit, fmt = *) numstrain, length
  ! numstrain is the number of strains in the sample
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
  allocate (singleton(4 * numstrain, 2), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: singleton!"
    ! Error, exit the program.
    stop
  end if
  allocate (tochange(numstrain, 2), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory for: tochange!"
    ! Error, exit the program.
    stop
  end if
  ! Create the format for printing sequence data.
  write (unit = sequence_format, fmt = "(a1,i0,a3)") "(", length, "a1)"
  ! read in the sequences
  do jstrain = 1, numstrain
    read (unit = population_unit, fmt = trim(sequence_format)) &
      (sequence(jstrain, knuc), knuc = 1, length)
  end do
  ! Read in pcrerror_dat
  ! pcrerror is the per nucleotide site error in pcr; according to Akmaev
  ! and Wang, it is 1.37e-4, or 1/7300
  read (unit = pcrerror_unit, fmt = *) pcrerror
  ! iii is the odd random number seed
  read (unit = pcrerror_unit, fmt = *) iii
  ! Initialize the random number generator.
  call randomInitialize (iii)
  call diverge (numstrain, length, sequence, divergematrix)
  numsingletons = 0
  do knuc = 1, length
    ! numnuc(1) = number of a's at site knuc
    ! numnuc(2) = number of t's at site knuc
    ! numnuc(3) = number of g's at site knuc
    ! numnuc(4) = number of c's at site knuc
    numnuc(1) = 0
    numnuc(2) = 0
    numnuc(3) = 0
    numnuc(4) = 0
    do nucnuc = 1, 4
      do jstrain = 1, numstrain
        if (sequence(jstrain, knuc) .eq. code(nucnuc)) then
          numnuc(nucnuc) = numnuc(nucnuc) + 1
        end if
      end do
      if (numnuc(nucnuc) .eq. 1) then
        do jstrain = 1, numstrain
          if (sequence(jstrain, knuc) .eq. code(nucnuc)) then
            numsingletons = numsingletons + 1
            singleton(numsingletons, 1) = jstrain
            singleton(numsingletons, 2) = knuc
            exit
          endif
        end do
      endif
    end do
  end do
  ! now, how many singletons to correct for pcr error?
  numtochange = nint (pcrerror * length * numstrain)
  if (numtochange .gt. numsingletons) then
    numtochange = numsingletons
  end if
  ! change the singletons
  if (numtochange .gt. 0) then
    do jchange = 1, numtochange
      call randomNumber (x)
      num = int (x * numsingletons) + 1
      tochange(jchange, 1) = singleton(num, 1)
      tochange(jchange, 2) = singleton(num, 2)
      singleton(num, 1) = singleton(numsingletons, 1)
      singleton(num, 2) = singleton(numsingletons, 2)
      numsingletons = numsingletons - 1
    end do
    ! now, change each of the chosen singletons to the value of the most
    ! similar strain
    closest = 1
    do jchange = 1, numtochange
      xmindiff = 1.0
      seqchange = tochange(jchange, 1)
      nucchange = tochange(jchange, 2)
      do iother = 1, numstrain
        if (seqchange .eq. iother) cycle
        diff = divergematrix(seqchange, iother)
        if (diff .lt. xmindiff) then
          closest = iother
          xmindiff = diff
        endif
      end do
      sequence(seqchange, nucchange) = sequence(closest, nucchange)
    end do
  end if
  ! Output to correctpcr_out
  write (unit = output_unit, fmt = *) numstrain, length
  do jstrain = 1, numstrain
    write (unit = output_unit, fmt = trim(sequence_format)) &
      (sequence(jstrain, knuc), knuc = 1, length)
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
  deallocate (singleton, stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: singleton!"
    ! Error, exit the program.
    stop
  end if
  deallocate (tochange, stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to deallocate memory for: tochange!"
    ! Error, exit the program.
    stop
  end if
  ! Close the random number generator.
  call randomClose ()
  ! Close data files.
  close (unit = population_unit)
  close (unit = pcrerror_unit)
  close (unit = output_unit)
  ! Successful termination of program.
  stop
end program correctpcr
