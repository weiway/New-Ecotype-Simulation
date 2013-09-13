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
! This is a hillclimbing version of finding the npop confidence interval.  For each value of npop tested,
! xn stays fixed at infinity, and the optimal values of omega and sigma are found for that npop.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program npopCI
  use Methods
  use NelderMead
  implicit none
  ! Local variables
  character(len = 256)                     :: inputFile
  character(len = 256)                     :: outputFile
  integer                                  :: ier
  integer                                  :: iprint
  integer                                  :: iquad
  integer                                  :: maxf
  integer                                  :: nloop
  integer                                  :: indexnpop
  integer, parameter                       :: outputUnit = 2
  integer, parameter                       :: nparams = 2
  double precision                         :: simp
  double precision                         :: stopcr
  double precision                         :: yvalue
  double precision                         :: params(nparams)
  double precision                         :: step(nparams)
  double precision                         :: var(nparams)
  real                                     :: probthreshold
  type(acinas_data)                        :: acinas
  type(number_increments)                  :: numincs
  type(parameters_data)                    :: bottom
  type(parameters_data)                    :: top
  type(parameters_data)                    :: parameters
  procedure(MinimizationFunction), pointer :: functn
  functn => fredprogram
  ! Read in command line arguments
  ! Expect path and filename for: npopIn.dat npopOut.dat
  if (iargc() .ge. 2) then
    call getArgument (1, inputFile)
    call getArgument (2, outputFile)
  else
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: npopIn.dat npopOut.dat"
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
  iprint = -1
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
  parameters%xn = bottom%xn
  parameters%omega = bottom%omega
  parameters%sigma = bottom%sigma
  do indexnpop = bottom%npop, top%npop, numincs%npop
    parameters%npop = indexnpop
    ! Note, for npop values besides the original one, 
    ! we will start with the omega and sigma values
    ! calculated for the previous npop value. 
    params(1) = log (parameters%omega)
    params(2) = log (parameters%sigma)
    step(1) = log (parameters%omega) / 2.0
    step(2) = log (parameters%sigma) / 2.0
    if (log (parameters%omega) .lt. 0.3 .and. log (parameters%omega) .gt. -0.3) step(1) = 0.15
    if (log (parameters%sigma) .lt. 0.3 .and. log (parameters%sigma) .gt. -0.3) step(2) = 0.15
    yvalue = 0.0
    call Minimization (params, step, nparams, yvalue, maxf, iprint, stopcr, nloop, iquad, simp, var, functn, ier, &
     outputUnit, probthreshold)
    write (unit = outputUnit, fmt = *) parameters%omega, parameters%sigma, parameters%npop, parameters%xn, yvalue
    if (-yvalue .lt. probthreshold) exit
  end do
  ! Close the output file.
  close (unit = outputUnit)
  stop

  contains

  subroutine fredprogram (params, yvalue)
    ! Subroutine parameters
    double precision, intent(inout) :: params(:)
    double precision, intent(out)   :: yvalue
    ! Local variables
    real :: avgsuccess(6)
    parameters%omega = exp (params(1))
    parameters%sigma = exp (params(2))
    if (parameters%omega .lt. 1.0e-7) parameters%omega = 1.0e-7
    if (parameters%sigma .lt. 1.0e-7) parameters%sigma = 1.0e-7
    ! avgsuccess%x500 is the average number of results that are within 500% tolerance
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

end program npopCI
