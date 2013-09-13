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
! Hill-climbing version of non-ultrametric design
! This is a non-ultrametric divergence design, with drift but no
! recombination. September 21, 2003
! note, this will soon be replaced by a quicker method that Danny
! is writing, which will not require the n-square binning algorithm,
! but rather an nlogn sorting routine
! note, the programs I've changed for nonultra are:
! main, whicheventandwhen, 
! divergenonultra (and deleted binning), binningdanny, coalescence,
! and startpops
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program hillclimb
  use Methods
  use NelderMead
  implicit none
  character(len = 256)                     :: inputFile
  character(len = 256)                     :: outputFile
  integer, parameter                       :: nparams = 4
  double precision                         :: simp
  double precision                         :: step(nparams)
  double precision                         :: stopcr
  double precision                         :: var(nparams)
  double precision                         :: params(nparams)
  double precision                         :: yvalue
  integer                                  :: ier
  integer                                  :: iprint
  integer                                  :: iquad
  integer                                  :: maxf
  integer                                  :: nloop
  integer, parameter                       :: outputUnit = 3
  type(acinas_data)                        :: acinas
  type(number_increments)                  :: numincs
  type(parameters_data)                    :: bottom
  type(parameters_data)                    :: top
  procedure(MinimizationFunction), pointer :: functn
  functn => fredprogram
  ! Read in command line arguments
  ! Expect path and filename for: hClimbIn.dat hClimbOut.dat
  if (iargc() .ge. 2) then
    call getArgument (1, inputFile)
    call getArgument (2, outputFile)
  else
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: hClimbIn.dat hClimbOut.dat"
    stop
  end if
  ! Check for optional debug command line argument.
  if (iargc() .eq. 3) then
    call getArgument (3, debug)
  else
    debug = .false.
  end if
  open (unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  call readAcinas (trim (inputFile), acinas, bottom, top, numincs)
  ! Set max. no. of function evaluations = MAXF, print every IPRINT.
  maxf = 100
  iprint = 1
  ! Set value for stopping criterion.   Stopping occurs when the
  ! standard deviation of the values of the objective function at
  ! the points of the current simplex < stopcr.
  stopcr = 1.0d-1
  nloop = 8
  ! Fit a quadratic surface to be sure a minimum has been found.
  iquad = 0
  ! As function value is being evaluated in DOUBLE PRECISION, it
  ! should be accurate to about 15 decimals.   If we set simp = 1.d-6,
  ! we should get about 9 dec. digits accuracy in fitting the surface.
  simp = 1.0d-6
  params(1) = log (bottom%omega)
  params(2) = log (bottom%sigma)
  params(3) = bottom%npop
  params(4) = log (bottom%xn)
  step(1) = log (bottom%omega) / 2.0
  step(2) = log (bottom%sigma) / 2.0
  step(3) = bottom%npop / 2.0
  step(4) = log (bottom%xn) / 2.0
  yvalue = 0.0
  ! Danny, yvalue is just a dummy variable that gets passed to my
  ! program in the CALL statement, but the meaningful passing of 
  ! variables occurs through the COMMON block
  ! yvalue is the value returned for success level determined by jwhichxavg
  call Minimization (params, step, nparams, yvalue, maxf, iprint, stopcr, nloop, iquad, simp, var, functn, ier, &
    outputUnit, acinas%probthreshold)
  close (unit = outputUnit)
  stop

  contains

  subroutine fredprogram (params, yvalue)
    ! Subroutine parameters
    double precision, intent(inout) :: params(:)
    double precision, intent(out)   :: yvalue
    ! Local variables
    real :: avgsuccess(6)
    type(parameters_data) :: parameters
    ! Make sure that npop is in the right range.
    if (params(3) .lt. 2) params(3) = 2
    if (params(3) .gt. acinas%nu) params(3) = acinas%nu
    parameters%omega = exp (params(1))
    parameters%sigma = exp (params(2))
    parameters%npop = nint (params(3))
    parameters%xn = exp (params(4))
    ! avgsuccess(1) is the average number of results that are within 500% tolerance
    ! for a particular set of parameter values, etc. 
    call simulation (acinas, parameters, avgsuccess)
    yvalue = -1.0 * avgsuccess(acinas%whichavg)
    if (debug) then
      call displayDebug (parameters)
      call displayDebug (yvalue)
      write (unit = *, fmt = *)
    end if
    return
  end subroutine fredprogram

end program hillclimb
