!     A program for function minimization using the simplex method.
!
!     For details, see Nelder and Mead, The Computer Journal, January 1965
!
!     Programmed by D.E. Shaw,
!     CSIRO, Division of Mathmatics and Statistics
!     P.O. Box 218, Lindfield, N.S.W. 2070
!
!     With amendments by R.W.M. Wedderburn
!     Rothamsted Experimental Station
!     Harpenden, Hertfordshire, England
!
!     Further amended by Alan Miller
!     CSIRO Division of Mathematics and Statistics
!     Private Bag 10, Clayton, Victoria 3168
!
!     Conversion to Fortran 2003 by Jason Wood
!     Montana State University
!     Bozeman, MT 59715

module NelderMead
  implicit none
  private

  public :: Minimization
  public :: MinimizationFunction

  interface
     subroutine MinimizationFunction (params, yvalue)
      double precision, intent(inout) :: params(:)
      double precision, intent(out)   :: yvalue
    end subroutine MinimizationFunction
  end interface

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Params:
  !    p()        = Input. Starting values of parameters.
  !                 Output. Final values of parameters.
  !    step()     = Input. Initial step sizes.
  !    nop        = Input. Number of parameters, including any to be held fixed.
  !    func       = Output. The function value corresponding to the final parameter values.
  !    max        = Input. The maximum number of function evaluations allowed. Say, 20 times the number of parameters.
  !    iprint     = Input. Print control parameter.
  !                 < 0 no printing.
  !                 = 0 printing of parameter values and the function value after initial evidence of convergence.
  !                 > 0 as for iprint = 0 plus progress reports after every iprint evaluations, plus printing
  !                 for the initial simplex.
  !    stopcr     = Input. Stopping criterion. The criterion is applied to the standard deviation of
  !                 the values of func at the points of the simplex.
  !    nloop      = Input. The stopping rule is applied after every nloop function evaluations.
  !                 Normally nloop should be slightly greater than nop, say nloop = 2 * nop.
  !    iquad      = Input.
  !                 = 1 if fitting of a quadratic surface is required.
  !                 = 0 if not.
  !                 N.B. The fitting of a quadratic surface is strongly recommended, provided that the fitted
  !                 function is continuous in the vicinity of the minimum. It is often a good indicator of
  !                 whether a premature termination of the search has occurred.
  !    simp       = Input. Criterion for expanding the simplex to overcome rounding errors before fitting the
  !                 quadratic surface. The simplex is expanded so that the function values at the points of 
  !                 the simplex exceed those at the supposed minimum by at least an amount simp.
  !    var()      = Output. Contains the diagonal elements of the inverse of the information matrix.
  !    functn     = Input. Name of the user's subroutine, which returns the function value for a given set of
  !                 parameter values in array p.
  !    ifault     = Output. 
  !                 = 0 for successful termination.
  !                 = 1 if maximum number of function evaluations exceeded.
  !                 = 2 if information matrix is not +VE semi-definite.
  !                 = 3 if nop < 1.
  !                 = 4 if nloop < .1
  !    lout       = Input. The file handle to output to.
  !    probthresh = Input. Stop when the solution reaches this level.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Minimization (p, step, nop, func, max, iprint, stopcr, nloop, iquad, simp, var, functn, ifault, lout, probthresh)
    integer, intent(in)                      :: nop
    integer, intent(in)                      :: max
    integer, intent(in)                      :: iprint
    integer, intent(in)                      :: nloop
    integer, intent(in)                      :: iquad
    integer, intent(out)                     :: ifault
    integer, intent(in)                      :: lout
    double precision, intent(out)            :: func
    double precision, intent(in)             :: stopcr
    double precision, intent(in)             :: simp
    double precision, intent(out)            :: var(:)
    double precision, intent(inout)          :: p(:)
    double precision, intent(inout)          :: step(:)
    real, intent(in)                         :: probthresh
    procedure(MinimizationFunction), pointer :: functn
    ! Local Variables
    double precision            :: bmat(nop * (nop + 1) / 2)
    double precision            :: vc(nop * (nop + 1) / 2)
    double precision            :: g(nop + 1, nop)
    double precision            :: h(nop + 1)
    double precision            :: pbar(nop)
    double precision            :: pstar(nop)
    double precision            :: pstst(nop)
    double precision            :: aval(nop)
    double precision            :: pmin(nop)
    double precision            :: a0
    double precision            :: calcprob
    double precision            :: hmin
    double precision            :: hmax
    double precision            :: hmean
    double precision            :: hstar
    double precision            :: hstd
    double precision            :: hstst
    double precision            :: rmax
    double precision            :: savemn
    double precision            :: test
    double precision            :: ymin
    double precision, parameter :: a = 1.0d0 ! Reflection Coefficient.
    double precision, parameter :: b = 0.5d0 ! Contraction Coefficient.
    double precision, parameter :: c = 2.0d0 ! Expansion Coefficient.
    integer                     :: i
    integer                     :: i1
    integer                     :: i2
    integer                     :: iflag
    integer                     :: ii
    integer                     :: ij
    integer                     :: ijk
    integer                     :: imin
    integer                     :: imax
    integer                     :: irank
    integer                     :: irow
    integer                     :: j
    integer                     :: j1
    integer                     :: jj
    integer                     :: k
    integer                     :: l
    integer                     :: loop
    integer                     :: nap
    integer                     :: neval
    integer                     :: nmore
    integer                     :: np1
    integer                     :: nullty
    ! Prepare formats
1000 format (' Progress Report every',I4,' function evaluations'/,' EVAL.   FUNC.VALUE.',10X,'PARAMETER VALUES')
1010 format (/1X, I4, 2X, G12.5, 2X, 5G11.4, 3(/21X, 5G11.4))
1020 format (' No. of function evaluations > ',I5)
1030 format (' RMS of function values of last simplex =',G14.6)
1040 format (' Centroid of last simplex =',4(/1X,6G13.5))
1050 format (' Function value at centroid =',G14.6)
1060 format (/' EVIDENCE OF CONVERGENCE')
1070 format (//' Minimum found after',I5,' function evaluations')
1080 format (' Minimum at',4(/1X,6G13.6))
1090 format (' Function value at minimum =',G14.6)
1110 format (/' Fitting quadratic surface about supposed minimum'/)
1120 format (/' MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.'/' MINIMUM PROBABLY NOT FOUND'/)
1130 format (/10X, 'Search restarting'/)
1140 format (' Minimum of quadratic surface =',G14.6,' at',4(/1X,6G13.5))
1150 format (' IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED', &
        & 1X,'FROM THE MINIMIZATION,'/ &
        & ' THE MINIMUM MAY BE FALSE &/OR THE INFORMATION MATRIX MAY BE', &
        & 1X,'INACCURATE'/)
1160 format (' Rank of information matrix =',I3/' Inverse of information matrix:-')
1170 format (/' If the function minimized was -LOG(LIKELIHOOD),'/ &
        & ' this is the covariance matrix of the parameters.'/ &
        & ' If the function was a sum of squares of residuals,'/ &
        & ' this matrix must be multiplied by twice the estimated', &
        & 1X, 'residual variance'/' to obtain the covariance matrix.'/)
1190 format (' INFORMATION MATRIX:-'/)
1200 format (//' CORRELATION MATRIX:-')
1210 format (/' A further',I4,' function evaluations have been used'/)
1230 format (1X,6G13.5)
1240 format (/)

    ! If progress reports have been requested, print the heading.
    if (iprint .gt. 0) write (unit = lout, fmt = 1000) iprint
    ! Check input arguments.
    ifault = 0
    if (nop .le. 0) ifault = 3
    if (nloop .le. 0) ifault = 4
    if (ifault .ne. 0) return
    ! Set nap = number of parameters to be varied, i.e. with step not equal to 0.
    nap = 0
    neval = 0
    loop = 0
    iflag = 0
    do i = 1, nop
      if (step(i) .ne. 0.d0) nap = nap + 1
    end do
    ! If nap = 0 evaluate function at the starting point and return.
    if (nap .eq. 0) then
      call functn (p, func)
      return
    end if
    ! Set up the initial simplex.
30  continue
    do i = 1, nop
      g(1, i) = p(i)
    end do
    irow = 2
    do i = 1, nop
      if (step(i) .eq. 0.d0) cycle
      do j = 1, nop
        g(irow, j) = p(j)
      end do
      g(irow, i) = p(i) + step(i)
      irow = irow + 1
    end do
    np1 = nap + 1
    do i = 1, np1
      do j = 1, nop
        p(j) = g(i, j)
      end do
        call functn (p, h(i))
        ! added July 2007.
	calcprob = -h(i)
	if (calcprob .ge. probthresh) then
	  func = -calcprob
	  return
	endif
        neval = neval + 1
        if (iprint .gt. 0) then
          write (unit = lout, fmt = 1010) neval, h(i), (p(j), j = 1, nop)
        end if
    end do
    ! Start of main cycle.
    ! Find max and min values for current simplex (hmax and hmin).
100 continue
    loop = loop + 1
    imax = 1
    imin = 1
    hmax = h(1)
    hmin = h(1)
    do i = 2, np1
      if (h(i) .gt. hmax) then
        imax = i
        hmax = h(i)
      else if (h(i) .lt. hmin) then
        imin = i
        hmin = h(i)
      end if
    end do
    ! Find the centroid of the vertices other than p(imax).
    do i = 1, nop
      pbar(i) = 0.d0
    end do
    do  i = 1, np1
      if (i .eq. imax) cycle
      do j = 1, nop
        pbar(j) = pbar(j) + g(i, j)
      end do
    end do
    do j = 1, nop
      pbar(j) = pbar(j) / float (nap)
    end do
    ! Reflect maximum through pbar to pstar.
    ! hstar = function value at pstar.
    do i = 1, nop
      pstar(i) = a * (pbar(i) - g(imax, i)) + pbar(i)
    end do
    call functn (pstar, hstar)
    ! added July 2007.
    calcprob = -hstar
    if (calcprob .ge. probthresh) then
      func = -calcprob
      return
    endif
    neval = neval + 1
    if (iprint .gt. 0) then
      if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) neval, hstar, (pstar(j), j = 1, nop)
    end if
    ! If hstar < hmin, reflect pbar through pstar.
    ! hstst = function value at pstst.
    if (hstar .lt. hmin) then
      do i = 1, nop
        pstst(i) = c * (pstar(i) - pbar(i)) + pbar(i)
      end do
      call functn (pstst, hstst)
      ! added July 2007.
      calcprob = -hstst
      if (calcprob .ge. probthresh) then
        func = -calcprob
        return
      endif
      neval = neval + 1
      if (iprint .gt. 0) then
        if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) neval, hstst, (pstst(j), j = 1, nop)
      end if
      ! If hstst < hmin replace current maximum point by pstst and hmax by hstst, then test for convergence.
      if (hstst .ge. hmin) goto 320
      do i = 1, nop
        if (step(i) .ne. 0.d0) g(imax, i) = pstst(i)
      end do
      h(imax) = hstst
      goto 340
    end if
    ! hstar is not < hmin.
    ! Test whether it is < function value at some point other than p(imax).
    ! If it is, replace p(imax) by pstar and hmax by hstar.
    do i = 1, np1
        if (i .eq. imax) cycle
        if (hstar .lt. h(i)) goto 320
    end do
    ! hstar > all function values except possibly hmax.
    ! If hstar <= hmax, replace p(imax) by pstar and hmax by hstar.
    if (hstar .le. hmax) then
      do i = 1, nop
        if (step(i) .ne. 0.d0) g(imax, i) = pstar(i)
      end do
      hmax = hstar
      h(imax) = hstar
    end if
    ! Contracted step to the point pstst.
    ! hstst = function value at pstst.
    do i = 1, nop
      pstst(i) = b * g(imax, i) + (1.0 - b) * pbar(i)
    end do
    call functn (pstst, hstst)
   ! added July 2007.
    calcprob = -hstst
    if (calcprob .ge. probthresh) then
      func = -calcprob
      return
    endif
    neval = neval + 1
    if (iprint .gt. 0) then
      if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) neval, hstst, (pstst(j), j = 1, nop)
    end if
    ! If hstst < hmax replace p(imax) by pstst and hmax by hstst.
    if (hstst .le. hmax) then
      do i = 1, nop
        if (step(i) .ne. 0.d0) g(imax, i) = pstst(i)
      end do
      h(imax) = hstst
      goto 340
    end if
    ! hstst > hmax.
    ! Shrink the simplex by replacing each point, other than the current minimum, by a point midway between
    ! its current position and the minimum.
    do i = 1, np1
      if (i .eq. imin) cycle
      do j = 1, nop
        if (step(j) .ne. 0.d0) g(i, j) = (g(i, j) + g(imin, j)) * 0.5
        p(j) = g(i, j)
      end do
      call functn (p, h(i))
      ! added July 2007.
      calcprob = -h(i)
      if (calcprob .ge. probthresh) then
        func = -calcprob
        return
      endif
      neval = neval + 1
      if (iprint .gt. 0) then
        if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) neval, h(i), (p(j), j = 1, nop)
      end if
    end do
    goto 340
320 continue
    ! Replace maximum point by pstar and h(imax) by hstar.
    do i = 1, nop
      if (step(i) .ne. 0.d0) g(imax, i) = pstar(i)
    end do
    h(imax) = hstar
    ! If loop = nloop test for convergence, otherwise repean main cycle.
340 continue
    if (loop .lt. nloop) goto 100
    ! Calculate mean and standard deviation of function values for the current simplex.
    hstd = 0.d0
    hmean = 0.d0
    do i = 1, np1
      hmean = hmean + h(i)
    end do
    hmean = hmean / float (np1)
    do i = 1, np1
      hstd = hstd + (h(i) - hmean) ** 2
    end do
    hstd = sqrt (hstd / float (np1))
    ! If the rms > stopcr, set iflag and loop to zero and goto the start of the main cycle again.
    if (hstd .gt. stopcr .or. neval .eq. max) then
      iflag = 0
      loop = 0
      goto 100
    end if
    ! Find the centroid of the current simplex and the function value there.
    do i = 1, nop
      if (step(i) .eq. 0.d0) cycle
      p(i) = 0.d0
      do j = 1, np1
        p(i) = p(i) + g(j, i)
      end do
        p(i) = p(i) / float (np1)
    end do
    call functn (p, func)
    ! added July 2007.
    calcprob = -func
    if (calcprob .ge. probthresh) then
      func = -calcprob
      return
    endif
    neval = neval + 1
    if (iprint .gt. 0) then
      if (mod (neval, iprint) .eq. 0) write (unit = lout, fmt = 1010) neval, func, (p(j), j = 1, nop)
    end if
    ! Test whether the number of function values allowed, max, has been overrun; if so, exit with ifault = 1.
    if (neval .gt. max) then
      ifault = 1
      if (iprint .gt. 0) then
        write (unit = lout, fmt = 1020) max
        write (unit = lout, fmt = 1030) hstd
        write (unit = lout, fmt = 1050) func
        write (unit = lout, fmt = 1040) (p(i), i = 1, nop)
      end if
      return
    end if
    ! Convergence criterion satisfied.
    ! If iflag = 0, set iflag and save hmean.
    ! If iflag = 1 and change in hmean <= stopcr then search is complete.
    if (iprint .gt. 0) then
      write (unit = lout, fmt = 1060)
      write (unit = lout, fmt = 1040) (p(i), i = 1, nop)
      write (unit = lout, fmt = 1050) func
    end if
    if (iflag .gt. 0) goto 450
    iflag = 1
440 continue
    savemn = hmean
    loop = 0
    goto 100
450 continue
    if (abs (savemn - hmean) .ge. stopcr) goto 440
    if (iprint .gt. 0) then
      write (unit = lout, fmt = 1070) neval
      write (unit = lout, fmt = 1080) (p(i), i = 1, nop)
      write (unit = lout, fmt = 1090) func
    end if
    if (iquad .le. 0) return
    !------------------------------------------------------------------
    ! Quadratic surface fitting.
    if (iprint .ge. 0) write (unit = lout, fmt = 1110)
    ! Expand the final simplex, if necessary, to overcome rounding errors.
    hmin = func
    nmore = 0
    do i = 1, np1
      test = abs (h(i) - func)
      if (test .ge. simp) cycle
      do j = 1, nop
        if (step(j) .ne. 0.d0) g(i, j) = (g(i, j) - p(j)) + g(i, j)
        pstst(j) = g(i, j)
      end do
      call functn (pstst, h(i))
      ! added July 2007.
      calcprob = -h(i)
      if (calcprob .ge. probthresh) then
        func = -calcprob
        return
      endif
      nmore = nmore + 1
      neval = neval + 1
      if (h(i) .ge. hmin) cycle
      hmin = h(i)
      if (iprint .ge. 0) write (unit = lout, fmt = 1010) neval, hmin, (pstst(j), j = 1, nop)
    end do
    ! Function values are calculated at an additional nap points.
    do i = 1, nap
      i1 = i + 1
      do j = 1, nop
        pstar(j) = (g(1, j) + g(i1, j)) * 0.5
      end do
      call functn (pstar, aval(i))
      ! added July 2007.
      calcprob = -aval(i)
      if (calcprob .ge. probthresh) then
        func = -calcprob
        return
      endif
      nmore = nmore + 1
      neval = neval + 1
    end do
    ! The matrix of estimated second derivatives is calculated and its lower triangle stored in bmat.
    a0 = h(1)
    do i = 1, nap
      i1 = i - 1
      i2 = i + 1
      if (i1 .lt. 1) cycle
      do j = 1, i1
        j1 = j + 1
        do k = 1, nop
          pstst(k) = (g(i2, k) + g(j1, k)) * 0.5
        end do
        call functn (pstst, hstst)
        ! added July 2007.
        calcprob = -hstst
       if (calcprob .ge. probthresh) then
          func = -calcprob
          return
        endif
        nmore = nmore + 1
        neval = neval + 1
        l = i * (i - 1) / 2 + j
        bmat(l) = 2.0 * (hstst + a0 - aval(i) - aval(j))
      end do
    end do
    l = 0
    do i = 1, nap
      i1 = i + 1
      l = l + i
      bmat(l) = 2.0 * (h(i1) + a0 - 2.0 * aval(i))
    end do
    ! The vector of estimated first derivatives is calculated and stored in aval.
    do i = 1, nap
      i1 = i + 1
      aval(i) = 2.0 * aval(i) - (h(i1) + 3.0 * a0) * 0.5
    end do
    ! The matrix q of Nelder and Mead is calculated and stored in g.
    do i = 1, nop
      pmin(i) = g(1, i)
    end do
    do i = 1, nap
      i1 = i + 1
      do j = 1, nop
        g(i1, j) = g(i1, j) - g(1, j)
      end do
    end do
    do i = 1, nap
      i1 = i + 1
      do j = 1, nop
        g(i, j) = g(i1, j)
      end do
    end do
    ! Invert bmat.
    call SymInv (bmat, nap, bmat, nullty, ifault, rmax)
    if (ifault .eq. 0) then
      irank = nap - nullty
    else
      if (iprint .ge. 0) write (unit = lout, fmt = 1120)
      ifault = 2
      if (neval .gt. max) return
      write (unit = lout, fmt = 1130)
      do i = 1, nop
        step(i) = 0.5 * step(i)
      end do
      goto 30
    end if
    ! bmat * a / 2 is calculated and stored in h.
    do i = 1, nap
      h(i) = 0.d0
      do j = 1, nap
        if (j .le. i) then
          l = i * (i - 1) / 2 + j
        else
          l = j * (j - 1) / 2 + i
        end if
          h(i) = h(i) + bmat(l) * aval(j)
      end do
    end do
    ! Find the position, pmin, and value, ymin, of the minimum of the quadratic.
    ymin = 0.d0
    do i = 1, nap
      ymin = ymin + h(i) * aval(i)
    end do
    ymin = a0 - ymin
    do i = 1, nop
      pstst(i) = 0.d0
      do j = 1, nap
        pstst(i) = pstst(i) + h(j) * g(j, i)
      end do
    end do
    do i = 1, nop
      pmin(i) = pmin(i) - pstst(i)
    end do
    if (iprint .gt. 0) then
      write (unit = lout, fmt = 1140) ymin, (pmin(i), i = 1, nop)
      write (unit = lout, fmt = 1150)
    end if
    ! q * bmat * q' / 2 is calculated and its lower triangle stored in vc.
    do i = 1, nop
      do j = 1, nap
        h(j) = 0.d0
        do k = 1, nap
          if (k .le. j) then
            l = j * (j - 1) / 2 + k
          else
            l = k * (k - 1) / 2 + j
          end if
          h(j) = h(j) + bmat(l) * g(k, i) * 0.5
        end do
      end do
      do j = i, nop
        l = j * (j - 1) / 2 + i
        vc(l) = 0.d0
        do k = 1, nap
          vc(l) = vc(l) + h(k) * g(k, j)
        end do
      end do
    end do
    ! The diagonal elements of vc are copied into var.
    j = 0
    do i = 1, nop
      j = j + i
      var(i) = vc(j)
    end do
    if (iprint .lt. 0) return
    write (unit = lout, fmt = 1160) irank
    ijk = 1
    goto 880
790 continue
    write (unit = lout, fmt = 1170)
    call SymInv (vc, nap, bmat, nullty, ifault, rmax)
    ! bmat now contains the information matrix.
    write (unit = lout, fmt = 1190)
    ijk = 3
    goto 880
800 continue
    ijk = 2
    ii = 0
    ij = 0
    do i = 1, nop
      ii = ii + i
      if (vc(ii) .gt. 0.d0) then
        vc(ii) = 1.0 / sqrt (vc(ii))
      else
        vc(ii) = 0.d0
      end if
      jj = 0
      do j = 1, i - 1
        jj = jj + j
        ij = ij + 1
        vc(ij) = vc(ij) * vc(ii) * vc(jj)
      end do
      ij = ij + 1
    end do
    write (unit = lout, fmt = 1200)
    ii = 0
    do i = 1, nop
      ii = ii + i
      if (vc(ii) .ne. 0.d0) vc(ii) = 1.d0
    end do
    goto 880
    ! Exit, on successful termination.
860 continue
    write (unit = lout, fmt = 1210) nmore
    return
880 continue
    l = 1
890 continue
    if (l .gt. nop) goto (790, 860, 800), ijk
    ii = l * (l - 1) / 2
    do i = l, nop
      i1 = ii + l
      ii = ii + i
      i2 = min (ii, i1 + 5)
      if (ijk .ne. 3) then
       write (unit = lout, fmt = 1230) (vc(j), j = i1, i2)
       cycle
      end if
      write (unit = lout, fmt = 1230) (bmat(j), j = i1, i2)
    end do
    write (unit = lout, fmt = 1240)
    l = l + 6
    goto 890
    return
  end subroutine Minimization

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Algorithm AS7, Applied Statistics, volume 17, 1968.
  !
  !  Params:
  !    a()    = Input. The symmetric matrix to be inverted, stored in lower triangular form.
  !    n      = Input. Order of the matrix.
  !    c()    = Output. The inverse of a (a generalized inverse if c is singular), also stored
  !             in lower triangular form. c and a may occupy the same locations.
  !    nullty = Output. The rank deficiency of a.
  !    ifault = Output. Error indicator.
  !             = 1 if n < 1
  !             = 2 if a is not +VE semi-definite
  !             = 0 otherwise
  !    rmax   = Output. Approximate bound on the accuracy of the diagonal elements of c.
  !             e.g. if rmax = 1.e-04 then the diagonal elements of c will be accurate to about 4 decimal digits.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine SymInv (a, n, c, nullty, ifault, rmax)
    double precision, intent(in)  :: a(:)
    double precision, intent(out) :: c(:)
    double precision, intent(out) :: rmax
    integer, intent(in)           :: n
    integer, intent(out)          :: nullty
    integer, intent(out)          :: ifault
    ! Local Variables
    double precision :: w(n)
    double precision :: x
    integer          :: i
    integer          :: icol
    integer          :: irow
    integer          :: j
    integer          :: jcol
    integer          :: k
    integer          :: l
    integer          :: ndiag
    integer          :: nn
    integer          :: nrow
    integer          :: mdiag
    logical          :: equal
    nrow = n
    ifault = 1
    if (nrow .le. 0) return
    ifault = 0
    ! Cholesky factorization of a, result in c.
    call Chola (a, nrow, c, nullty, ifault, rmax, w)
    if (ifault .ne. 0) return
    ! Invert c and form the product (cinv)' * cinv, where cinv is the inverse of c, row by row starting with the last row.
    ! irow = The row number.
    ! ndiag = Location of last element in the row.
    nn = nrow * (nrow + 1) / 2
    irow = nrow
    ndiag = nn
    do
      if (c(ndiag) .ne. 0.d0) then
        l = ndiag
        do i = irow, nrow
          w(i) = c(l)
          l = l + i
        end do
        icol = nrow
        jcol = nn
        mdiag = nn
        do
          l = jcol
          x = 0.d0
          if (icol .eq. irow) x = 1.d0 / w(irow)
          k = nrow
          do
            if (k .eq. irow) exit
            x = x - w(k) * c(l)
            k = k - 1
            l = l - 1
            if (l .gt. mdiag) l = l - k + 1
          end do
          c(l) = x / w(irow)
          if (icol .eq. irow) goto 80
          mdiag = mdiag - icol
          icol = icol - 1
          jcol = jcol - 1
        end do
      end if
      l = ndiag
      do j = irow, nrow
        c(l) = 0.d0
        l = l + j
      end do
80    continue
      ndiag = ndiag - irow
      irow = irow - 1
      if (irow .eq. 0) exit
    end do
    return
  end subroutine SymInv

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Algorithm AS6, Applied Statistics, Volume 17, 1968, with modifications by A.J. Miller
  !
  !  Params:
  !    a()    = Input. A +VE definite matrix stored in lower-triangular form.
  !    n      = Input. The order of a
  !    u()    = Output. A lower triangular matrix such that u * u' = a. a and u may occupy the same locations.
  !    nullty = Output. The rank deficiency of a.
  !    ifault = Output. Error indicator
  !             = 1 if n < 1
  !             = 2 if a is not +VE semi-definite
  !             = 0 otherwise
  !    rmax   = Output. An estimate of the relative accuracy of the diagonal elements of u.
  !    r()    = Output. Array containing bounds on the relative accuracy of each diagonal element of u.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Chola (a, n, u, nullty, ifault, rmax, r)
    double precision, intent(in)  :: a(:)
    double precision, intent(out) :: u(:)
    double precision, intent(out) :: rmax
    double precision, intent(out) :: r(:)
    integer, intent(in)           :: n
    integer, intent(out)          :: nullty
    integer, intent(out)          :: ifault
    ! Local Variables
    double precision, parameter :: eta = epsilon(1.d0)
    double precision            :: rsq
    double precision            :: w
    integer                     :: i
    integer                     :: icol
    integer                     :: irow
    integer                     :: j
    integer                     :: k
    integer                     :: l
    integer                     :: m
    ifault = 1
    if (n .le. 0) return
    ifault = 2
    nullty = 0
    rmax = eta
    r(1) = eta
    j = 1
    k = 0
    ! Factorize column by column, icol = column number.
    do icol = 1, n
      l = 0
      ! irow = row number within column icol.
      do irow = 1, icol
        k = k + 1
        w = a(k)
        if (irow .eq. icol) rsq = (w * eta) ** 2
        m = j
        do i = 1, irow
          l = l + 1
          if (i .eq. irow) exit
          w = w - u(l) * u(m)
          if (irow .eq. icol) rsq = rsq + (u(l) ** 2 * r(i)) ** 2
          m = m + 1
        end do
        if (irow .eq. icol) goto 50
        if (u(l) .eq. 0.d0) goto 30
        u(k) = w / u(l)
        cycle
30      continue
        u(k) = 0.d0
        if (abs (w) .gt. abs (rmax * a(k))) return
      end do
      ! End of row, estimate relative accuracy of diagonal element.
50    continue
      rsq = sqrt (rsq)
      if (abs (w) .le. 5.0 * rsq) goto 60
      if (w .lt. 0.d0) return
      u(k) = sqrt (w)
      r(i) = rsq / w
      if (r(i) .gt. rmax) rmax = r(i)
      goto 70
60    continue
      u(k) = 0.d0
      nullty = nullty + 1
70    continue
      j = j + icol
    end do
    ifault = 0
    return
  end subroutine Chola

end module NelderMead
