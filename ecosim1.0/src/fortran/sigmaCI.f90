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
! As of Feb. 3, 2007, this is a hillclimbing version of
! finding the npop confidence interval.  For each value 
! of npop tested, xn stays fixed at infinity, and
! the optimal values of omega and sigma are found
! for that npop.
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
! params(1)=omega, params(2)=sigma
! npop and xn are fixed, and are passed through the
! common block called 'parameters'
! nparams is the number of parameters, in this case, 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program sigmaCI
  use Methods
  use NelderMead
  implicit none
  ! Local variables
  character(len = 256)                     :: inputFile
  character(len = 256)                     :: outputFile
  integer                                  :: xincs
  integer                                  :: ier
  integer                                  :: iprint
  integer                                  :: iquad
  integer                                  :: maxf
  integer                                  :: nloop
  integer                                  :: indexsigma
  integer, parameter                       :: nparams = 2
  integer, parameter                       :: outputUnit = 2
  double precision                         :: simp
  double precision                         :: stopcr
  double precision                         :: yvalue
  double precision                         :: params(nparams)
  double precision                         :: step(nparams)
  double precision                         :: var(nparams)
  real                                     :: probthreshold
  real                                     :: diffsigma
  type(acinas_data)                        :: acinas
  type(number_increments)                  :: numincs
  type(parameters_data)                    :: bottom
  type(parameters_data)                    :: top
  type(parameters_data)                    :: parameters
  procedure(MinimizationFunction), pointer :: functn
  functn => fredprogram
  ! Read in command line arguments
  ! Expect path and filename for: sigmaIn.dat sigmaOut.dat
  if (iargc() .ge. 2) then
    call getArgument (1, inputFile)
    call getArgument (2, outputFile)
  else
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: sigmaIn.dat sigmaOut.dat"
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
  ! Fixed factor for drift:
  parameters%xn = bottom%xn
  ! Starting values for omega and npop:
  parameters%omega = bottom%omega
  parameters%npop = bottom%npop
  ! if xincs=0, it means that only one set of omega values will be used
  xincs = 1
  diffsigma = log (top%sigma) - log (bottom%sigma)
  if (numincs%sigma .eq. 0) then
    top%sigma = bottom%sigma * 10.0
    diffsigma = log (top%sigma) - log (bottom%sigma)
    numincs%sigma = 1
    xincs = 0
  endif
  do indexsigma = 0, numincs%sigma
    parameters%sigma = exp (log (bottom%sigma) + indexsigma * (diffsigma / numincs%sigma) * 0.999)
    params(1) = log (parameters%omega)
    params(2) = parameters%npop
    step(1) = log (parameters%omega) / 2.0
    step(2) = parameters%npop / 2.0
    if (log (parameters%omega) .lt. 0.30 .and. log (parameters%omega) .gt. -0.30) step(1) = 0.15
    yvalue = 0.0
    ! yvalue is the negative of the likelihood
    call Minimization (params, step, nparams, yvalue, maxf, iprint, stopcr, nloop, iquad, simp, var, functn, ier, &
      outputUnit, probthreshold)
    write (unit = outputUnit, fmt = *) parameters%omega, parameters%sigma, parameters%npop, parameters%xn, yvalue
    if (-yvalue .lt. probthreshold .or. xincs .eq. 0) exit
  end do
  ! Close output file.
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
    parameters%npop = nint (params(2))
    if (parameters%omega .lt. 1.0e-7) parameters%omega = 1.0e-7
    if (parameters%npop .lt. 1) parameters%npop = 1
    if (parameters%npop .gt. acinas%nu) parameters%npop = acinas%nu
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

end program sigmaCI
