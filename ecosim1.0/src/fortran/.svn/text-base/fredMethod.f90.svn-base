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
! On September 16, 2006, I now invent a new strain from the
! parent population when a new ecotype is created.  This new
! strain is the parent founder of the new ecotype.
! These are the programs that were changed between Sept. 21, 2005,
! and Sept. 16, 2006:
! donicheinv, choosepopnicheinv, coalescence, startpops
! Also, some of the CALL statements in the main programs were 
! changed.
! This is a non-ultrametric divergence design, with drift but no
! recombination. September 21, 2005
! note, this will soon be replaced by a quicker method that Danny
! is writing, which will not require the n-square binning algorithm,
! but rather an nlogn sorting routine
! note, the programs I've changed for nonultra are:
! main, whicheventandwhen, 
! divergenonultra (and deleted binning), binningdanny, coalescence,
! and startpops
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program fredMethod
  use Methods
  implicit none
  ! Local variables
  character(len = 256)    :: inputFile
  character(len = 256)    :: outputFile
  character(len = 256)    :: timeFile
  integer                 :: i
  integer                 :: indexxn
  integer                 :: indexsigma
  integer                 :: indexomega
  integer                 :: indexnpop
  integer, parameter      :: outputUnit = 2
  integer, parameter      :: timeUnit = 4
  real                    :: diffomega
  real                    :: diffsigma
  real                    :: diffnpop
  real                    :: diffxn
  real                    :: time1
  real                    :: time2
  real                    :: delta
  real                    :: avgsuccess(6)
  type(acinas_data)       :: acinas
  type(number_increments) :: numincs
  type(parameters_data)   :: bottom
  type(parameters_data)   :: top
  type(parameters_data)   :: parameters
  ! Read in command line arguments
  ! Expect path and filename for: acinas.dat acinas.out time.out
  if (iargc() .ge. 3) then
    call getArgument (1, inputFile)
    call getArgument (2, outputFile)
    call getArgument (3, timeFile)
    debug = .false.
  else
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: acinas.dat acinas.out time.out"
    stop
  end if
  ! Check for optional debug command line argument.
  if (iargc() .eq. 4) then
    call getArgument (4, debug)
  else
    debug = .false.
  end if
  ! Open the output files.
  open (unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  open (unit = timeUnit, file = trim (timeFile), access = 'sequential', form = 'formatted')
  call cpu_time (time1)
  ! Read in the input file.
  call readAcinas (trim (inputFile), acinas, bottom, top, numincs)
  diffomega = log10 (top%omega) - log10 (bottom%omega)
  diffsigma = log10 (top%sigma) - log10 (bottom%sigma)
  diffnpop = log10 (real (top%npop)) - log10 (real (bottom%npop))
  diffxn = log10 (top%xn) - log10 (bottom%xn)
  if (numincs%omega .eq. 0) then
    top%omega = bottom%omega + 1.0
    diffomega = 2.0 * (log10 (top%omega) - log10 (bottom%omega))
    numincs%omega = 1
  endif
  if (numincs%sigma .eq. 0) then
    top%sigma = bottom%sigma + 1.0
    diffsigma = 2.0 * (log10 (top%sigma) - log10 (bottom%sigma))
    numincs%sigma = 1
  endif
  if (numincs%npop .eq. 0) then
    top%npop = bottom%npop + 1
    diffnpop = 2.0 * (log10 (real (top%npop)) - log10 (real (bottom%npop)))
    numincs%npop = 1
  endif
  if (numincs%xn .eq. 0) then
    top%xn = bottom%xn + 1.0
    diffxn = 2.0 * (log10 (top%xn) - log10 (bottom%xn))
    numincs%xn = 1
  endif
  indexomega = 0
  do while (indexomega .lt. numincs%omega)
    parameters%omega = 10.0 ** (log10 (bottom%omega) + indexomega * (diffomega / numincs%omega) * 0.999)
    indexsigma = 0
    do while (indexsigma .lt. numincs%sigma)
      parameters%sigma = 10.0 ** (log10 (bottom%sigma) + indexsigma * (diffsigma / numincs%sigma) * 0.999)
      indexnpop = 0
      do while (indexnpop .lt. numincs%npop)
        parameters%npop = round (10.0 ** (log10 (real (bottom%npop)) + indexnpop * (diffnpop / numincs%npop) * 0.999))
        if (parameters%npop .gt. acinas%nu) parameters%npop = acinas%nu
        indexxn = 0
        do while (indexxn .lt. numincs%xn)
          parameters%xn = 10.0 ** (log10 (bottom%xn) + indexxn * (diffxn / numincs%xn) * 0.999)
          ! Run the simulation.
          call simulation (acinas, parameters, avgsuccess)
          ! If the simulation is successful, write to the output file.
          if (avgsuccess(1) .gt. 0.0) then
            write (unit = outputUnit, fmt = "(f13.5,',',f13.5,',',i5,',',1e9.2,6(',',1e10.3))") &
              parameters%omega, parameters%sigma, parameters%npop, parameters%xn, (avgsuccess(i), i = 1, 6)
              ! avgsuccess is the number of results that are within a tolerance for a particular set of parameter values.
          end if
          ! Display debug data if requested.
          if (debug) then
            call displayDebug (parameters)
            call displayDebug (avgsuccess)
            write (unit = *, fmt = *)
          end if
          indexxn = indexxn + 1
        end do
        indexnpop = indexnpop + 1
      end do
      indexsigma = indexsigma + 1
    end do
    indexomega = indexomega + 1
  end do
  ! delta is the time elapsed
  call cpu_time (time2)
  delta = time2 - time1
  write (unit = timeUnit, fmt = *) delta, '    time'
  ! Close the output files.
  close (unit = outputUnit);
  close (unit = timeUnit);
  stop
end program fredMethod
