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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program demarcationsCI
  use Methods
  implicit none
  ! Local variables
  character(len = 256)                     :: inputFile
  character(len = 256)                     :: outputFile
  integer                                  :: istep
  integer                                  :: indexnpop
  integer, parameter                       :: nparams = 4
  integer, parameter                       :: outputUnit = 2
  double precision                         :: params(nparams)
  double precision                         :: yvalue
  type(acinas_data)                        :: acinas
  type(number_increments)                  :: numincs
  type(parameters_data)                    :: bottom
  type(parameters_data)                    :: top
  type(parameters_data)                    :: parameters
  ! Read in command line arguments
  ! Expect path and filename for: demarcationsIn.dat demarcationsOut.dat
  if (iargc() .ge. 2) then
    call getArgument (1, inputFile)
    call getArgument (2, outputFile)
  else
    write (unit = *, fmt = *) "Invalid number of paramenters supplied!"
    write (unit = *, fmt = *) "Expected: demarcationsIn.dat demarcationsOut.dat"
    stop
  end if
  ! Check for optional debug command line argument.
  if (iargc() .eq. 3) then
    call getArgument (3, debug)
  else
    debug = .false.
  end if
  open (unit = outputUnit, file = trim (outputFile), access = 'sequential', form = 'formatted')
  !Print *, "CALL READACINAS"
  call readAcinas (trim (inputFile), acinas, bottom, top, numincs)
  !  As function value is being evaluated in DOUBLE PRECISION, it
  !  should be accurate to about 15 decimals.   If we set simp = 1.d-6,
  !  we should get about 9 dec. digits accuracy in fitting the surface.
  !Print *, "FINISHED READACINAS"
  parameters%omega = bottom%omega
  parameters%sigma = bottom%sigma
  parameters%xn = bottom%xn
  istep = numincs%npop
  if (istep .le. 0) istep = 1
  do indexnpop = bottom%npop, top%npop, istep
    parameters%npop = indexnpop
    params(1) = log (parameters%omega)
    params(2) = log (parameters%sigma)
    params(3) = parameters%npop
    params(4) = log (parameters%xn)
    yvalue = 0.0
    ! yvalue is the negative of the likelihood
	!Print *, "FREDPROGRAM"
    call fredprogram (params, yvalue)
	!Print *, "FINISHED FREDPROGRAM"
    write (unit = outputUnit, fmt = *) params(3), ',', yvalue
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
    ! Make sure that npop is in the right range.
    if (params(3) .lt. 1.0) params(3) = 1.0
    if (params(3) .gt. acinas%nu) params(3) = acinas%nu
    parameters%omega = exp (params(1))
    parameters%sigma = exp (params(2))
    parameters%npop = nint (params(3))
    parameters%xn = exp (params(4))
    ! avgsuccess%x500 is the average number of results that are within 500% tolerance
    ! for a particular set of parameter values, etc. 
	!Print *, "SIMULATION"
    call simulation (acinas, parameters, avgsuccess)
	!Print *, "FINISHED SIMULATION"
    yvalue = -1.0 * avgsuccess(acinas%whichavg)
    if (debug) then
      call displayDebug (parameters)
      call displayDebug (yvalue)
      write (unit = *, fmt = *)
    end if
    return
  end subroutine fredprogram

end program demarcationsCI
