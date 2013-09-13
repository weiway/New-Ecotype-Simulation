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
! We assume a PCR error rate (r) of *** (*ref).  
! The number of errors expected over the whole data set 
! is then rln, where l is the length of the segment in 
! nucleotides, and n is the number of sequences sampled.  
! We assumed that each PCR error would be unique nucleotide 
! at one site.  Of all the sites with a unique nucleotide, 
! we randomly chose rln of them to be corrected.  Each of 
! these chosen was corrected to the nucleotide state at that site 
! for the strain most closely related.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program correctpcr
  use Methods
  implicit none
  character(len = 20)             :: sequence_format
  character(len = 1), allocatable :: sequence(:,:)
  character(len = 1), parameter   :: code(4) = (/'a', 't', 'g', 'c'/)
  character(len = 256)            :: population_dat
  character(len = 256)            :: pcrerror_dat
  character(len = 256)            :: correctpcr_out
  integer                         :: allocate_status
  integer                         :: max_size
  integer                         :: numnuc(4)
  integer                         :: iii
  integer                         :: num
  integer                         :: numtochange
  integer                         :: closest
  integer                         :: seqchange
  integer                         :: nucchange
  integer                         :: numtouse
  integer                         :: numsingletons
  integer                         :: iother
  integer                         :: jstrain
  integer                         :: kstrain
  integer                         :: knuc
  integer                         :: nucnuc
  integer                         :: jchange
  integer                         :: numstrain
  integer                         :: length
  integer                         :: nuctype
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
  ! Read command line arguments.
  if (iargc() .ge. 3) then
    ! Expect filenames for population.dat, pcrerror.dat, and correctpcr.out.
    call getArgument (1, population_dat)
    call getArgument (2, pcrerror_dat)
    call getArgument (3, correctpcr_out)
  else if (iargc() .eq. 0) then
    ! No filenames supplied, so use defaults.
    write (unit = *, fmt = *) "Using default file names."
    population_dat = 'population.dat'
    pcrerror_dat = 'pcrerror.dat'
    correctpcr_out = 'correctpcr.out'
  else
    ! An unexpected number of arguments was supplied.
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: population.dat pcrerror.dat correctpcr.out"
    stop
  end if
  ! Check for optional debug command line argument.
  if (iargc() .eq. 4) then
    call getArgument (4, debug)
  else
    debug = .false.
  end if
  ! Verify that input files exist.
  inquire (file = trim (population_dat), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The population.dat file was not found at: ", trim (population_dat)
    stop
  end if
  inquire (file = trim (pcrerror_dat), exist = file_exists)
  if (file_exists .neqv. .true.) then
    write (unit = *, fmt = *) "The pcrerror.dat file was not found at: ", trim (pcrerror_dat)
    stop
  end if
  ! Open input files
  open (unit = population_unit, file = trim (population_dat), access = 'sequential', form = 'formatted')
  open (unit = pcrerror_unit, file = trim (pcrerror_dat), access = 'sequential', form = 'formatted')
  ! Open output file
  open (unit = output_unit, file = trim (correctpcr_out), access = 'sequential', form = 'formatted')
  ! Read in population_dat
  read (unit = population_unit, fmt = *) numstrain, length
  ! numstrain is the number of strains in the sample
  ! length is the number of nucleotide sites
  max_size = numstrain * length
  ! Allocate variables.
  allocate (sequence(numstrain, length), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (divergematrix(numstrain, numstrain), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (singleton(numstrain, max_size), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  allocate (tochange(numstrain, max_size), stat = allocate_status)
  if (allocate_status .gt. 0) then
    write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
  end if
  call initRandomSeed (iii)
  ! Create the format for printing sequence data.
  write (unit = sequence_format, fmt = "(a1,i0,a3)") "(", length, "a1)"
  ! read in the sequences
  do jstrain = 1, numstrain
    read (unit = population_unit, fmt = trim(sequence_format)) (sequence(jstrain, knuc), knuc = 1, length)
  end do
  ! Read in pcrerror_dat
  ! pcrerror is the per nucleotide site error in pcr; according to Akmaev and Wang, it is 1.37e-4, or 1/7300
  read (unit = pcrerror_unit, fmt = *) pcrerror
  ! iii is the odd random number seed
  read (unit = pcrerror_unit, fmt = *) iii
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
      call random_number (x)
      num = int (x * numsingletons) + 1
      tochange(jchange, 1) = singleton(num, 1)
      tochange(jchange, 2) = singleton(num, 2)
      singleton(num, 1) = singleton(numsingletons, 1)
      singleton(num, 2) = singleton(numsingletons, 2)
      numsingletons = numsingletons - 1
    end do
    ! now, change each of the chosen singletons to the value of the most similar strain
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
    write (unit = output_unit, fmt = trim(sequence_format)) (sequence(jstrain, knuc), knuc = 1, length)
  end do
  ! Close data files.
  close (unit = population_unit)
  close (unit = pcrerror_unit)
  close (unit = output_unit)
  stop
end program correctpcr
