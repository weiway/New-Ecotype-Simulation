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
!  Methods module contains all of the common methods used in Ecotype Simulation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Methods
  implicit none
  private

  ! Acinas data type.
  type acinas_data
    integer              :: numcrit
    integer              :: nu
    integer              :: nrep
    integer              :: lengthseq
    integer              :: whichavg
    integer, allocatable :: realdata(:)
    real                 :: probthreshold
    real, allocatable    :: crit(:)
  end type acinas_data

  ! Number of increments data type.
  type number_increments
    integer :: npop
    integer :: omega
    integer :: sigma
    integer :: xn
  end type number_increments

  ! Parameters data type.
  type parameters_data
    integer :: npop
    real    :: omega
    real    :: sigma
    real    :: xn
  end type parameters_data

  ! Declare public methods.
  public :: displayDebug
  public :: diverge
  public :: getArgument
  public :: initRandomSeed
  public :: readAcinas
  public :: round
  public :: simulation

  ! Declare public types.
  public :: acinas_data
  public :: number_increments
  public :: parameters_data

  ! Declare public variables.
  logical, public :: debug

  ! Declare global parameters.
  integer, parameter :: EVENT_NICHE_INVASION     = 100
  integer, parameter :: EVENT_PERIODIC_SELECTION = 101
  integer, parameter :: EVENT_POPULATION_DRIFT   = 102

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Displays debug data for the AverageSuccess, Parameters, and Yvalue variable types.
  !
  !  Params:
  !    debugData             : Data to display.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface displayDebug
    module procedure displayAverageSuccess
    module procedure displayParameters
    module procedure displayYValue
  end interface displayDebug

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets an argument from the command line.
  !
  !  Params:
  !    arg                   : Argument number to get.
  !    out                   : Value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface getArgument
    module procedure getLogicalArgument
    module procedure getIntegerArgument
    module procedure getRealArgument
    module procedure getStringArgument
  end interface getArgument

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Complete linkage hierarchical clustering based sequence identity.
  !
  !  Params:
  !    nclustersarray(i)     : The number of clusters (bins) at identity level sought.
  !    nu                    : The number of sequences.
  !    numcrit               : The number of criteria for making cluster bins.
  !    crit(i,j)             : List of bin identity levels (between 0.0 and 1.0).
  !    identitymatrix(i,j)   : Matrix of clusters.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine binningdanny (nclustersarray, nu, numcrit, crit, identitymatrix)
    integer, intent(inout) :: nclustersarray(:)
    integer, intent(in)    :: nu
    integer, intent(in)    :: numcrit
    real, intent(in)       :: crit(:)
    real, intent(in)       :: identitymatrix(:,:)
    ! Local variables.
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: xi
    integer :: xj
    integer :: n_clusters
    integer :: cluster_size(nu)
    integer :: cluster(nu, nu)
    real    :: x
    real    :: identity_level
    real    :: cluster_dist(nu, nu)
    ! Loop: for each identity level find number of bins
    do l = 1, numcrit
      ! Initialize variables
      ! identity_level = max identity difference in cluster
      ! cluster(i,j) = jth element of ith cluster
      ! cluster_size(i) = size of ith cluster
      ! cluster_dist(i,j) = min percent identity seq in cluster i vs j
      ! n_clusters = number of clusters
      identity_level = crit(l)
      do i = 1, nu
        do j = 1, nu
          cluster(i, j) = 0
          cluster_dist(i, j) = identitymatrix(i, j)
        end do
        cluster(i, 1) = i
        cluster_size(i) = 1
      end do
      n_clusters = nu
      ! Loop: find closest pair of clusters, combine, recompute distances
      do
        ! Find closest pair, x = distance between pair xi,xj; xi < xj
        x = -1.0
        xi = 1
        xj = 2
        do i = 1, nu
          do j = i + 1, nu
            if (cluster_dist(i, j) .lt. x) cycle
              x = cluster_dist(i, j)
              xi = i
              xj = j
          end do
        end do
        ! If closest pair further apart than identity_level then done
        if (x .lt. identity_level) exit
        ! Otherwise, combine the closest pair into xi and eliminate xj
        do i = cluster_size(xi) + 1, cluster_size(xi) + cluster_size(xj)
          cluster(xi, i) = cluster(xj, i - cluster_size(xi))
        end do
        cluster_size(xi) = cluster_size(xi) + cluster_size(xj)
        n_clusters = n_clusters - 1
        do i = 1, nu
          cluster_dist(xj, i) = -1.0
          cluster_dist(i, xj) = -1.0
          cluster(xj, i) = 0
        end do
        cluster_size(xj) = 0
        ! Recalculate distance from xi to all other active clusters and repeat
        do i = 1, nu
          if (cluster_size(i) .eq. 0) cycle
          x = 1.0
          do j = 1, cluster_size(xi)
            do k = 1, cluster_size(i)
              if (identitymatrix(cluster(xi, j), cluster(i, k)) .lt. x) then
                x = identitymatrix(cluster(xi, j), cluster(i, k))
              end if
            end do
          end do
          cluster_dist(xi, i) = x
          cluster_dist(i, xi) = x
        end do
      end do
      nclustersarray(l) = n_clusters
    end do
    return
  end subroutine binningdanny

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Canonical
  !
  !  Params:
  !    npop                  : The number of ecotypes assumed to be in the environmental DNA sample.
  !    nu                    : The number of sequences.
  !    numstrain(i)          :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine canonical (npop, nu, numstrain)
    integer, intent(in)    :: npop
    integer, intent(in)    :: nu
    integer, intent(inout) :: numstrain(:)
    ! Local variables.
    integer :: upto
    integer :: numberofstrains
    integer :: istrain
    integer :: ipop
    real    :: x
    numstrain(1) = nu
    if (npop .gt. 1) then
      ! for the upto-1 population, choose a random population (ipop)
      do upto = 2, npop
        do
          call random_number (x)
          ipop = int ((upto - 1) * x) + 1
          numberofstrains = numstrain(ipop)
          if (numberofstrains .ne. 1) exit
        end do
        ! next, pick a random breakpoint
        call random_number (x)
        istrain = int ((numberofstrains - 1) * x) + 1
        numstrain(ipop) = istrain
        numstrain(upto) = numberofstrains - istrain
      end do
    end if
    return
  end subroutine canonical

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Coalescence
  !
  !  Params:
  !    lengthseq             : The length in nucleotides of the sequence analyzed.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namecoalesce(i)       :
  !    namedesc(i,j)         :
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numcoalesce           :
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    div(i)                :
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine coalescence (lengthseq, nameanc, namecoalesce, namedesc, numanc, numanctot, numcoalesce, numdesc, ancdist, div, time)
    integer, intent(in)    :: lengthseq
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(in)    :: namecoalesce(:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(in)    :: numcoalesce
    integer, intent(inout) :: numdesc(:)
    real, intent(inout)    :: ancdist(:,:)
    real, intent(inout)    :: div(:)
    real, intent(in)       :: time
    ! Local variables.
    integer :: strainname
    integer :: idivaddactual
    integer :: jcoal
    integer :: iancnow
    integer :: nameofdesc
    integer :: namelastancest
    integer :: js
    integer :: numofanc
    integer :: numancest
    integer :: jdesc
    real    :: divaddexpect
    real    :: divaddactual
    numanctot = numanctot + 1
    ! we add one to the total number of ancestors NUMANCTOT
    ! next, for the new ancestor being created (number numanctot),
    ! we must list all its descendants
    numdesc(numanctot) = 0
    ! numdesc(numanctot) is the number of descendants of the ancestor
    ! currently being described
    do jcoal = 1, numcoalesce
      strainname = namecoalesce(jcoal)
      numancest = numanc(strainname)
      namelastancest = nameanc(strainname, numancest)
      ! namelastancest is the name of strainname's most distant ancestor,
      ! until now
      do js = 1, numdesc(namelastancest)
        numdesc(numanctot) = numdesc(numanctot) + 1
        nameofdesc = namedesc(namelastancest, js)
        namedesc(numanctot, numdesc(numanctot)) = nameofdesc
        ! next, give nameofdesc its oldest ancestor
        numanc(nameofdesc) = numanc(nameofdesc) + 1
        nameanc(nameofdesc, numanc(nameofdesc)) = numanctot
      end do
    end do
    div(numanctot) = time / lengthseq
    ! note do-loop below is an addition for nonultra
    do jcoal = 1, numcoalesce
      strainname = namecoalesce(jcoal)
       ! note, this penultimate ancestor is:
      iancnow = nameanc(strainname, numanc(strainname) - 1)
      divaddexpect = (div(numanctot) - div(iancnow))
      ! now get poisson--note we have to deal with integer number
      ! of mutations, so I multiply the per nucleotide divergence
      ! by lengthseq, and then divide later
      call poisson(divaddexpect * lengthseq, idivaddactual)
      divaddactual = real (idivaddactual) / lengthseq
      ! divaddactual is the actual number of substitutions
      ! in the node going to the last ancestor from a lineage,
      ! compared to the previous ancestor of strain strainname
      ! now we deal with all the descendants of ancestor iancnow
      do jdesc = 1, numdesc(iancnow)
        nameofdesc = namedesc(iancnow, jdesc)
        numofanc = numanc(nameofdesc)
        ancdist(nameofdesc, numofanc) = divaddactual + ancdist(nameofdesc, numofanc - 1)
      end do
    end do
    return
  end subroutine coalescence

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  displayAverageSuccess
  !
  !  Part of the displayData interface.
  !
  !  Params:
  !    debugData             : Data to display.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine displayAverageSuccess (debugData)
    real, intent(in) :: debugData(6)
    ! Local variables
    integer :: i
    write (unit = *, fmt = *) 'avgSuccess = ', (debugData(i), i = 1, 6)
  end subroutine displayAverageSuccess

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  displayParameters
  !
  !  Part of the displayData interface.
  !
  !  Params:
  !    debugData             : Data to display.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine displayParameters (debugData)
    type(parameters_data), intent(in) :: debugData
    write (unit = *, fmt = *) '     omega = ', debugData%omega
    write (unit = *, fmt = *) '     sigma = ', debugData%sigma
    write (unit = *, fmt = *) '      npop = ', debugData%npop
    write (unit = *, fmt = *) '        xn = ', debugData%xn     
  end subroutine displayParameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  displayYValue
  !
  !  Part of the displayData interface.
  !
  !  Params:
  !    debugData             : Data to display.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine displayYValue (debugData)
    double precision, intent(in) :: debugData
    write (unit = *, fmt = *) '   Y Value = ', debugData
  end subroutine displayYValue

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Diverge
  !
  !  Params:
  !    numstrains    :
  !    length        :
  !    seq(i,j)           :
  !    divergematrix(i,j) :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diverge (numstrains, length, seq, divergematrix)
    character(len = 1), intent(in) :: seq(:,:)
    integer, intent(in)            :: numstrains
    integer, intent(in)            :: length
    real, intent(inout)            :: divergematrix(:,:)
    ! Local variables.
    integer :: jstrain, kstrain, knuc
    real    :: diff
    do jstrain = 1, numstrains
      do kstrain = 1, numstrains
        diff = 0.0
        do knuc = 1, length
          if (seq(jstrain, knuc) .ne. seq(kstrain, knuc)) then
            diff = diff + 1.0
          end if
        end do
        divergematrix(jstrain, kstrain) = diff / length
      end do
    end do
    return
  end subroutine diverge

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Divergenonultra
  !
  !  Params:
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    nu                    : The number of sequences.
  !    numanc                : The number of ancestors for the ith strain.
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    identitymatrix(i,j)   : Matrix of clusters.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine divergenonultra (nameanc, nu, numanc, ancdist, identitymatrix)
    integer, intent(in) :: nameanc(:,:)
    integer, intent(in) :: nu
    integer, intent(in) :: numanc(:)
    real, intent(in)    :: ancdist(:,:)
    real, intent(inout) :: identitymatrix(:,:)
    ! Local variables.
    integer :: jstrain
    integer :: kstrain
    integer :: janc
    integer :: kanc
    real    :: xjukesdivergence
    real    :: divergence
    logical :: equal
    do jstrain = 1, nu
      do kstrain = 1, nu
        do janc = 1, numanc(jstrain)
          equal = .false.
          do kanc = 1, numanc(kstrain)
            if (nameanc(jstrain, janc) .eq. nameanc(kstrain, kanc)) then
              equal = .true.
              exit
            end if 
          end do
          if (equal) exit
        end do
        divergence = ancdist(jstrain, janc) + ancdist(kstrain, kanc)
        xjukesdivergence = -0.75 * (exp ((-4.0 / 3.0) * divergence) - 1.0)
        identitymatrix(jstrain, kstrain) = 1.0 - xjukesdivergence
      end do
    end do
    return
  end subroutine divergenonultra

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  doNicheInvasion
  !
  !  Params:
  !    activename(i)         : The activename is the name from the list of still-active populations.
  !    activepop             : The number of active populations.
  !    lengthseq             : The length in nucleotides of the sequence analyzed.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namedesc(i,j)         :
  !    namestrain(i,j)       :
  !    ntotalpop             :
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numcoalesce           :
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    numstrain(i)          :
  !    popfornascent         :
  !    realname(i)           : The realname of a population is its original number.
  !    strainparent          :
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    div(i)                :
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doNicheInvasion (activename, activepop, lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, &
    numanctot, numcoalesce, numdesc, numstrain, popfornascent, realname, strainparent, ancdist, div, time)
    integer, intent(inout) :: activename(:)
    integer, intent(inout) :: activepop
    integer, intent(in)    :: lengthseq
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(in)    :: namestrain(:,:)
    integer, intent(out)   :: ntotalpop
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(in)    :: numcoalesce
    integer, intent(inout) :: numdesc(:)
    integer, intent(in)    :: numstrain(:)
    integer, intent(in)    :: popfornascent
    integer, intent(inout) :: realname(:)
    integer, intent(in)    :: strainparent
    real, intent(inout)    :: ancdist(:,:)
    real, intent(inout)    :: div(:)
    real, intent(in)       :: time
    ! Local variables.
    integer :: namecoalesce(numcoalesce)
    integer :: jpop
    ! Now do the coalescence of the whole nascent population (popfornascent) and strainparent
    namecoalesce = namestrain(popfornascent, 1:numstrain(popfornascent))
    namecoalesce(numcoalesce) = strainparent
    call coalescence (lengthseq, nameanc, namecoalesce, namedesc, numanc, numanctot, numcoalesce, numdesc, ancdist, div, time)
    ! Next, delete the nascent population
    realname(activename(popfornascent)) = realname(activepop)
    activename(realname(activepop)) = activename(popfornascent)
    activepop = activepop - 1
    ntotalpop = 0
    do jpop = 1, activepop
      ntotalpop = ntotalpop + numstrain(realname(jpop))
    end do
    return
  end subroutine doNicheInvasion

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  doPeriodicSelection
  !
  !  Params:
  !    activepop             : The number of active populations.
  !    lengthseq             : The length in nucleotides of the sequence analyzed.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namedesc(i,j)         :
  !    namestrain(i,j)       :
  !    ntotalpop             :
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numcoalesce           : 
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    numstrain(i)          :
  !    popforps              :
  !    realname(i)           : The realname of a population is its original number.
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    div(i)                :
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doPeriodicSelection (activepop, lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, numanctot, &
    numcoalesce, numdesc, numstrain, popforps, realname, ancdist, div, time)
    integer, intent(inout) :: activepop
    integer, intent(in)    :: lengthseq
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(in)    :: namestrain(:,:)
    integer, intent(inout) :: ntotalpop
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(in)    :: numcoalesce
    integer, intent(inout) :: numdesc(:)
    integer, intent(inout) :: numstrain(:)
    integer, intent(in)    :: popforps
    integer, intent(in)    :: realname(:)
    real, intent(inout)    :: ancdist(:,:)
    real, intent(inout)    :: div(:)
    real, intent(in)       :: time
    ! Local variables.
    integer :: namecoalesce(numcoalesce)
    integer :: jpop
    ! First, make an array of the strains that will be coalesced
    namecoalesce = namestrain(popforps, 1:numstrain(popforps))
    call coalescence (lengthseq, nameanc, namecoalesce, namedesc, numanc, numanctot, numcoalesce, numdesc, ancdist, div, time)
    numstrain(popforps) = 1
    ntotalpop = 0
    do jpop = 1, activepop
      ntotalpop = ntotalpop + numstrain(realname(jpop))
    end do
    return
  end subroutine doPeriodicSelection

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  doPopulationDrift
  !
  !  Params:
  !    lengthseq             : The length in nucleotides of the sequence analyzed.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namedesc(i,j)         :
  !    namestrain(i,j)       :
  !    ntotalpop             :
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numcoalesce           :
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    numstrain(i)          :
  !    popfordrift           :
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    div(i)                :
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine doPopulationDrift (lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, numanctot, numcoalesce, &
    numdesc, numstrain, popfordrift, ancdist, div, time)
    integer, intent(in)    :: lengthseq
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(inout) :: namestrain(:,:)
    integer, intent(out)   :: ntotalpop
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(in)    :: numcoalesce
    integer, intent(inout) :: numdesc(:)
    integer, intent(inout) :: numstrain(:)
    integer, intent(in)    :: popfordrift
    real, intent(inout)    :: ancdist(:,:)
    real, intent(inout)    :: div(:)
    real, intent(in)       :: time
    ! Local variables.
    integer :: driftee
    integer :: drifter
    integer :: namecoalesce(numcoalesce)
    real    :: x
   ! Pick two strains from population popfordrift
    call random_number (x)
    driftee = int (x * numstrain(popfordrift)) + 1
    do
      call random_number (x)
      drifter = int (x * numstrain(popfordrift)) + 1
      if (drifter .ne. driftee) exit
    end do
    namecoalesce(1) = namestrain(popfordrift, drifter)
    namecoalesce(2) = namestrain(popfordrift, driftee)
    call coalescence (lengthseq, nameanc, namecoalesce, namedesc, numanc, numanctot, numcoalesce, numdesc, ancdist, div, time)
    namestrain(popfordrift, driftee) = namestrain(popfordrift, numstrain(popfordrift))
    numstrain(popfordrift) = numstrain(popfordrift) - 1
    ntotalpop = ntotalpop - 1
    return
  end subroutine doPopulationDrift

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Returns the number of populations eligible to be the parent ecotypes, as
  !  well as eligible for periodic selection. for a population to be eligible to be a nascent
  !  ecotype, it can have any number of strain members. For a population to be eligible for
  !  periodic selection, it must have more than one strain member.
  !
  !  Params:
  !    activepop             : The number of active populations.
  !    eligibleparent        :
  !    eligibleps            :
  !    numstrain(i)          :
  !    realname(i)           : The realname of a population is its original number.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eligible (activepop, eligibleparent, eligibleps, numstrain, realname)
    integer, intent(in)  :: activepop
    integer, intent(out) :: eligibleparent
    integer, intent(out) :: eligibleps
    integer, intent(in)  :: numstrain(:)
    integer, intent(in)  :: realname(:)
    ! Local variables.
    integer :: jpop
    eligibleparent = activepop
    if (activepop .eq. 1) then
      eligibleparent = 0
    endif
    eligibleps = 0
    do jpop = 1, activepop
      if (numstrain(realname(jpop)) .ne. 1) then
        eligibleps = eligibleps + 1
      end if
    end do
    return
  end subroutine eligible

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  eventNicheInvasion
  !
  !  Params:
  !    activename(i)         : The activename is the name from the list of still-active populations.
  !    activepop             : The number of active populations.
  !    lengthseq             : The length in nucleotides of the sequence analyzed.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namedesc(i,j)         :
  !    namestrain(i,j)       :
  !    ntotalpopwithinvented :
  !    nu                    : The number of sequences.
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    numstrain(i)          :
  !    realname(i)           : The realname of a population is its original number.
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    div(i)                :
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eventNicheInvasion (activename, activepop, lengthseq, nameanc, namedesc, namestrain, ntotalpop, &
    ntotalpopwithinvented, nu, numanc, numanctot, numdesc, numstrain, realname, ancdist, div, time)
    integer, intent(inout) :: activename(:)
    integer, intent(inout) :: activepop
    integer, intent(in)    :: lengthseq
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(inout) :: namestrain(:,:)
    integer, intent(inout) :: ntotalpop
    integer, intent(out)   :: ntotalpopwithinvented
    integer, intent(in)    :: nu
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(inout) :: numdesc(:)
    integer, intent(inout) :: numstrain(:)
    integer, intent(inout) :: realname(:)
    real, intent(inout)    :: ancdist(:,:)
    real, intent(inout)    :: div(:)
    real, intent(in)       :: time
    ! Local variables.
    integer :: chosen
    integer :: numeligible
    integer :: numcoalesce
    integer :: popfornascent
    integer :: popforparent
    integer :: jpop
    integer :: eligibleparentpop(activepop - 1)
    real    :: x
    ! first, choose the population for the nascent population
    call random_number (x)
    chosen = int (x * activepop) + 1
    popfornascent = realname(chosen)
    ! above, we use the first strain of the chosen population as the strain
    ! next, choose the population for the parent population
    numeligible = 0
    do jpop = 1, activepop
      if (realname(jpop) .eq. popfornascent) cycle
      numeligible = numeligible + 1
      eligibleparentpop(numeligible) = realname(jpop)
    end do
    call random_number (x)
    chosen = int (x * numeligible) + 1
    popforparent = eligibleparentpop(chosen)
    ! now we create a parent strain from population popforparent
    numstrain(popforparent) = numstrain(popforparent) + 1
    chosen = numstrain(popforparent)
    ntotalpopwithinvented = ntotalpopwithinvented + 1
    namestrain(popforparent, chosen) = ntotalpopwithinvented
    numanctot = numanctot + 1
    numanc(ntotalpopwithinvented) = 1
    nameanc(ntotalpopwithinvented, 1) = numanctot
    ancdist(ntotalpopwithinvented, 1) = 0
    numdesc(numanctot) = 1
    namedesc(numanctot, 1) = ntotalpopwithinvented
    div(numanctot) = 0
    numcoalesce = numstrain(popfornascent) + 1
    call doNicheInvasion (activename, activepop, lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, &
      numanctot, numcoalesce, numdesc, numstrain, popfornascent, realname, ntotalpopwithinvented, ancdist, div, time)
    return
  end subroutine eventNicheInvasion

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  eventPeriodicSelection
  !
  !  Params:
  !    activepop             : The number of active populations.
  !    lengthseq             : The length in nucleotides of the sequence analyzed.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namedesc(i,j)         :
  !    namestrain(i,j)       :
  !    ntotalpop             :
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    numstrain(i)          :
  !    realname(i)           : The realname of a population is its original number.
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    div(i)                :
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eventPeriodicSelection (activepop, lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, numanctot, &
    numdesc, numstrain, realname, ancdist, div, time)
    integer, intent(inout) :: activepop
    integer, intent(in)    :: lengthseq
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(in)    :: namestrain(:,:)
    integer, intent(inout) :: ntotalpop
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(inout) :: numdesc(:)
    integer, intent(inout) :: numstrain(:)
    integer, intent(in)    :: realname(:)
    real, intent(inout)    :: ancdist(:,:)
    real, intent(inout)    :: div(:)
    real, intent(in)       :: time
    ! Local variables.
    integer :: chosen
    integer :: eligiblepspop(activepop)
    integer :: jpop
    integer :: numeligible
    integer :: numcoalesce
    integer :: popforps
    real    :: x
    numeligible = 0
    do jpop = 1, activepop
      if (numstrain(realname(jpop)) .le. 1) cycle
      numeligible = numeligible + 1
      eligiblepspop(numeligible) = realname(jpop)
    end do
    call random_number (x)
    chosen = int (x * numeligible) + 1
    popforps = eligiblepspop(chosen)
    numcoalesce = numstrain(popforps)
    call doPeriodicSelection (activepop, lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, numanctot, &
      numcoalesce, numdesc, numstrain, popforps, realname, ancdist, div, time)
    return
  end subroutine eventPeriodicSelection

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  eventPopulationDrift
  !
  !  Params:
  !    activepop             : The number of active populations.
  !    lengthseq             : The length in nucleotides of the sequence analyzed.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namedesc(i,j)         :
  !    namestrain(i,j)       :
  !    ntotalpop             :
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    numstrain(i)          :
  !    realname(i)           : The realname of a population is its original number.
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    div(i)                :
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !    xn                    : The actual population size in nature.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine eventPopulationDrift (activepop, lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, numanctot, &
    numdesc, numstrain, realname, ancdist, div, time, xn, popfordrift)
    integer, intent(in)    :: activepop
    integer, intent(in)    :: lengthseq
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(inout) :: namestrain(:,:)
    integer, intent(out)   :: ntotalpop
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(inout) :: numdesc(:)
    integer, intent(inout) :: numstrain(:)
    integer, intent(in)    :: realname(:)
    integer, intent(out)   :: popfordrift
    real, intent(inout)    :: ancdist(:,:)
    real, intent(inout)    :: div(:)
    real, intent(in)       :: time
    real, intent(in)       :: xn
    ! Local variables.
    integer :: jpop
    real    :: x
    real    :: totaldrift
    real    :: probdrift(4000)
    totaldrift = 0.0
    do jpop = 1, activepop
      probdrift(realname(jpop)) = totaldrift + (1.0 / xn) * numstrain(realname(jpop)) * (numstrain(realname(jpop)) - 1.0) * 0.5
      totaldrift = probdrift(realname(jpop))
    end do
    if (totaldrift .eq. 0) return
    do jpop = 1, activepop
      probdrift(realname(jpop)) = probdrift(realname(jpop)) / totaldrift
    end do
    call random_number (x)
    do jpop = 1, activepop
      if (x .lt. probdrift(realname(jpop))) then
        popfordrift = realname(jpop)
        return
      endif
    end do
    return
  end subroutine eventPopulationDrift

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Perform the fred method for the given values.
  !
  !  Params:
  !    acinas                : The acinas data for this simulation.
  !    parameters            : The parameters for this simulation (omega, sigma, npop, xn).
  !    success               : The success status of this simulation.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fredMethod (acinas, parameters, success)
    type(acinas_data), intent(in)     :: acinas
    type(parameters_data), intent(in) :: parameters
    logical, intent(out)              :: success(6)
    ! Local variables.
    integer :: event
    integer :: activepop
    integer :: ntotalpop
    integer :: ntotalpopwithinvented
    integer :: numanctot
    integer :: nclustersarray(acinas%numcrit)
    integer :: numstrain(parameters%npop)
    integer :: realname(parameters%npop)
    integer :: nameanc(4000, 4000)                      ! bl1 common block
    integer :: namedesc(4000, 4000)                     ! bl1 common block
    integer :: namestrain(4000, 4000)                   ! bl1 common block
    integer :: numanc(4000)
    integer :: numdesc(4000)
    integer :: activename(4000)                         ! candidate for removal.  initialized in startpops as index = value (1..npop), only used in doNicheInvasion
    integer :: numcoalesce
    integer :: popfordrift
    real    :: time
    real    :: ancdist(2000, 2000)                      ! bl10 common block
    real    :: div(4000)
    real    :: identitymatrix(acinas%nu, acinas%nu)     ! bl10 common block
    ! Start the simulation.
    call startpops (activename, activepop, nameanc, namedesc, namestrain, parameters%npop, ntotalpop, ntotalpopwithinvented, &
      acinas%nu, numanc, numanctot, numdesc, numstrain, realname, ancdist, time)
    do
      ! subroutine WhichEventAndWhen returns the name of the key event, and gives the amount of time since the previous key event
      call whichEventAndWhen (parameters, activepop, event, numstrain, realname, time)
      select case (event)
        case (EVENT_NICHE_INVASION)
          call eventNicheInvasion (activename, activepop, acinas%lengthseq, nameanc, namedesc, namestrain, ntotalpop, &
            ntotalpopwithinvented, acinas%nu, numanc, numanctot, numdesc, numstrain, realname, ancdist, div, time)
        case (EVENT_PERIODIC_SELECTION)
          call eventPeriodicSelection (activepop, acinas%lengthseq, nameanc, namedesc, namestrain, ntotalpop, &
            numanc, numanctot, numdesc, numstrain, realname, ancdist, div, time)
        case (EVENT_POPULATION_DRIFT)
          call eventPopulationDrift (activepop, acinas%lengthseq, nameanc, namedesc, namestrain, ntotalpop, &
            numanc, numanctot, numdesc, numstrain, realname, ancdist, div, time, parameters%xn, popfordrift)
          numcoalesce = 2
          call doPopulationDrift (acinas%lengthseq, nameanc, namedesc, namestrain, ntotalpop, numanc, numanctot, numcoalesce, &
            numdesc, numstrain, popfordrift, ancdist, div, time)
      end select
      if (ntotalpop .le. 1) exit
    end do
    call divergenonultra (nameanc, acinas%nu, numanc, ancdist, identitymatrix)
    ! next, do binning
    call binningdanny (nclustersarray, acinas%nu, acinas%numcrit, acinas%crit, identitymatrix)
    ! Initialize success values.
    success = .false.
    ! Check whether a replicate is within a factor of (5.00, 2.00, 1.50, 1.25, 1.10, 1.05) for all bins
    success(1) = testForSuccessFit (nclustersarray, acinas%numcrit, acinas%realdata, 5.00)
    if (success(1)) success(2) = testforsuccessfit (nclustersarray, acinas%numcrit, acinas%realdata, 2.00)
    if (success(2)) success(3) = testforsuccessfit (nclustersarray, acinas%numcrit, acinas%realdata, 1.50)
    if (success(3)) success(4) = testforsuccessfit (nclustersarray, acinas%numcrit, acinas%realdata, 1.25)
    if (success(4)) success(5) = testforsuccessfit (nclustersarray, acinas%numcrit, acinas%realdata, 1.10)
    if (success(5)) success(6) = testforsuccessfit (nclustersarray, acinas%numcrit, acinas%realdata, 1.05)
    return
  end subroutine fredMethod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets an integer argument from the command line.
  !
  !  Params:
  !    arg                   : Argument number to get.
  !    out                   : Integer value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getIntegerArgument (arg, out)
    integer, intent(in)  :: arg
    integer, intent(out) :: out
    ! Local variables.
    character(len = 100) :: buffer
    integer              :: error
    call getarg (arg, buffer)
    read (unit = buffer, fmt = *, iostat = error) out
    if (error .ne. 0) then
      write (unit = *, fmt = *) "ERROR: Integer number not supplied with argument ", arg
      stop
    end if
    return
  end subroutine getIntegerArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets a logical argument from the command line.
  !
  !  Params:
  !    arg                   : Argument number to get.
  !    out                   : Logical value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getLogicalArgument (arg, out)
    integer, intent(in)  :: arg
    logical, intent(out) :: out
    character(len = 100) :: buffer
    call getarg (arg, buffer)
    ! Unless we receive an appropriate true string, return false.
    out = .false.
    if (trim(buffer) .eq. 'true' .or. trim(buffer) .eq. 't') then
      out = .true.
    else if (trim(buffer) .eq. 'True' .or. trim(buffer) .eq. 'T') then
      out = .true.
    end if
    return
  end subroutine getLogicalArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets a real argument from the command line.
  !
  !  Params:
  !    arg                   : Argument number to get.
  !    out                   : Real value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getRealArgument (arg, out)
    integer, intent(in) :: arg
    real, intent(out)   :: out
    ! Local variables.
    character(len = 100) :: buffer
    integer              :: error
    call getarg (arg, buffer)
    read (unit = buffer, fmt = *, iostat = error) out
    if (error .ne. 0) then
      write (unit = *, fmt = *) "ERROR: Real number not supplied with argument ", arg
      stop
    end if
    return
  end subroutine getRealArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Gets a string argument from the command line.
  !
  !  Params:
  !    arg                   : Argument number to get.
  !    out                   : String value to return.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getStringArgument (arg, out)
    integer, intent(in)             :: arg
    character(len = *), intent(out) :: out
    call getarg (arg, out)
    return
  end subroutine getStringArgument

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Initialize the random number generator with the given value of iii.
  !
  !  Params:
  !    iii                   : The seed.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initRandomSeed (iii)
    integer, intent(in) :: iii
    ! Local variables.
    integer              :: i
    integer              :: n
    integer              :: allocate_status
    integer, allocatable :: seed(:)
    call random_seed (size = n)
    ! Allocate local variables.
    allocate (seed(n), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
    end if
    seed = (/ (iii - i, i = 1, n) /)
    call random_seed (put = seed)
  end subroutine initRandomSeed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Takes the expected number of substitutions, and gives back the actual number of substitutions using the Poisson distribution.
  !
  !  Params:
  !    xmean                 :
  !    nsubs                 :
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson (xmean, nsubs)
    integer, intent(out) :: nsubs
    real, intent(in)     :: xmean
    ! Local variables.
    integer :: jmut
    real    :: x
    real    :: accumprob
    real    :: prob
    real    :: expect
    expect = xmean
    call random_number (x)
    accumprob = 0.0
    prob = exp (-1.0 * expect)
    do jmut = 0, 100
      if (jmut .ne. 0) then
        prob = (prob * expect) / jmut
      end if
      accumprob = accumprob + prob
      if (x .lt. accumprob) then
        nsubs = jmut
        return
      end if
    end do
    nsubs = int (expect)
    return
  end subroutine poisson

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Read in the data for this simulation.
  !
  !  Params:
  !    fname                 : File name to load.
  !    acinas                : The acinas data for this simulation.
  !    bottom                : Bottom parameter values.
  !    top                   : Top parameter values.
  !    numincs               : The number of increments for the parameters.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readAcinas (fname, acinas, bottom, top, numincs)
    character(len = *), intent(in)       :: fname
    type(acinas_data), intent(out)       :: acinas
    type(parameters_data), intent(out)   :: bottom
    type(parameters_data), intent(out)   :: top
    type(number_increments), intent(out) :: numincs
    ! Local variables.
    integer            :: allocate_status
    integer            :: iii
    integer            :: jcrit
    integer            :: io_status
    integer, parameter :: io_unit = 101
    ! Open the file
    open (unit = io_unit, file = fname, action = 'read', access = 'sequential', form = 'formatted')
    ! numcrit is the number of criteria for making cluster bins
    read (unit = io_unit, fmt = *) acinas%numcrit
    ! Allocate variables.
    allocate (acinas%crit(acinas%numcrit), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
    end if
    allocate (acinas%realdata(acinas%numcrit), stat = allocate_status)
    if (allocate_status .gt. 0) then
      write (unit = *, fmt = *) "Failed to allocate memory!", allocate_status
    end if
    do jcrit = 1, acinas%numcrit
      ! realdata(jcrit) is the number of bins in the real data at criterion jcrit / 1000.0
      read (unit = io_unit, fmt = *) acinas%realdata(jcrit)
    end do
    ! realdata(jcrit) is the actual number of bins for the jcrit th criterion
    ! crit(jcrit) is the jcrit th criterion value
    do jcrit = 1, acinas%numcrit
      read (unit = io_unit, fmt = *) acinas%crit(jcrit)
    end do
    ! omega is the rate of niche invasion per eligible parental population
    ! measured as niche invasions per nucleotide substitution in a given gene
    ! omegabot and omegatop are the range
    read (unit = io_unit, fmt = *) bottom%omega, top%omega
    ! sigma is the rate of periodic selection per eligible population,
    ! measured as periodic selection events per population per nucleotide substitution in
    ! a given gene
    read (unit = io_unit, fmt = *) bottom%sigma, top%sigma
    ! npop is the number of ecotypes assumed to be in the environmental DNA sample
    read (unit = io_unit, fmt = *) bottom%npop, top%npop
    ! note that omega, sigma, and npop are the three parameters we will estimate by
    ! maximum likelihood
    ! xn is the actual population size in nature
    read (unit = io_unit, fmt = *) bottom%xn, top%xn
    ! numincs values are the numbers of increments to be investigated
    read (unit = io_unit, fmt = *) numincs%omega, numincs%sigma, numincs%npop, numincs%xn
    ! nu is the number of homologous gene sequences in the environmental sample
    ! following Acinas et al., this should be in the thousands.  The program is
    ! currently written to allow 10000 samples
    read (unit = io_unit, fmt = *) acinas%nu
    ! nrep is the number of replicate simulations for a given set of sigma, omega, and npop
    read (unit = io_unit, fmt = *) acinas%nrep
    ! iii is the odd random number seed (up to nine digits)
    read (unit = io_unit, fmt = *) iii
    call initRandomSeed (iii)
    ! lengthseq is the length in nucleotides of the sequence analyzed
    read (unit = io_unit, fmt = *) acinas%lengthseq
    ! whichavg indicates success level considered:
    !   1 represents 5.00 times
    !   2 represents 2.00 times
    !   3 represents 1.50 times
    !   4 represents 1.25 times
    !   5 represents 1.10 times
    !   6 represents 1.05 times
    read (unit = io_unit, fmt = *, iostat = io_status) acinas%whichavg
    if (io_status .ne. 0) then
      acinas%whichavg = 1
    end if
    ! get the number of individuals in each ecotype
    read (unit = io_unit, fmt = *, iostat = io_status) acinas%probthreshold
    if (io_status .ne. 0) then
      acinas%probthreshold = 1.0
    end if
    ! Close the file
    close (unit = io_unit)
    return
  end subroutine readAcinas

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Round the given real number to the closest integer.
  !
  !  Params:
  !    x                     : The real number to round.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function round (x)
    real, intent(in) :: x
    ! Local variables.
    integer :: i
    i = int (x)
    if (x - i .ge. 0.5) then
      i = i + 1
    end if
    round = i
  end function round

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Perform nrep repetitions of the fred method, and calculate the average number of results
  !  within (500%, 200%, 150%, 125%, 110%, and 105%) tolerance for a particular set of parameter values.
  !
  !  Params:
  !    acinas                : The acinas data for this simulation.
  !    parameters            : The parameters for this simulation (omega, sigma, npop, xn).
  !    avgsuccess            : The average success for this run.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine simulation (acinas, parameters, avgsuccess)
    type(acinas_data), intent(in)      :: acinas
    type(parameters_data), intent(in)  :: parameters
    real, intent(out)                  :: avgsuccess(6)
    ! Local variables.
    integer :: i
    integer :: irep
    logical :: success(6)
    ! Initialize avgsuccess array.
    avgsuccess = 0.0
    do irep = 1, acinas%nrep
      call fredMethod (acinas, parameters, success)
      ! Count success at each level.
      do i = 1, 6
        if (success(i)) avgsuccess(i) = avgsuccess(i) + 1.0
      end do
    end do
    ! Convert avgsuccess counts to averages.
    avgsuccess = avgsuccess / acinas%nrep
    return
  end subroutine simulation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Startpops
  !
  !  Params:
  !    activename(i)         : The activename is the name from the list of still-active populations.
  !    activepop             : The number of active populations.
  !    nameanc(i,j)          : The name of the jth ancestor of the ith strain.
  !    namedesc(i,j)         :
  !    namestrain(i,j)       :
  !    npop                  : The number of ecotypes assumed to be in the environmental DNA sample.
  !    ntotalpop             :
  !    ntotalpopwithinvented :
  !    nu                    : The number of sequences.
  !    numanc(i)             : The number of ancestors for the ith strain.
  !    numanctot             : The total number of ancestors, over all strains.
  !    numdesc(i)            : The number of descendants of the ith ancestor.
  !    numstrain(i)          :
  !    realname(i)           : The realname of a population is its original number.
  !    ancdist(i,j)          : The actual divergence between contemporary strain i and its jth ancestor.
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine startpops (activename, activepop, nameanc, namedesc, namestrain, npop, ntotalpop, ntotalpopwithinvented, nu, &
    numanc, numanctot, numdesc, numstrain, realname, ancdist, time)
    integer, intent(inout) :: activename(:)
    integer, intent(out)   :: activepop
    integer, intent(inout) :: nameanc(:,:)
    integer, intent(inout) :: namedesc(:,:)
    integer, intent(inout) :: namestrain(:,:)
    integer, intent(in)    :: npop
    integer, intent(out)   :: ntotalpop
    integer, intent(out)   :: ntotalpopwithinvented
    integer, intent(in)    :: nu
    integer, intent(inout) :: numanc(:)
    integer, intent(out)   :: numanctot
    integer, intent(inout) :: numdesc(:)
    integer, intent(inout) :: numstrain(:)
    integer, intent(inout) :: realname(:)
    real, intent(inout)    :: ancdist(:,:)
    real, intent(out)      :: time
    ! Local variables.
    integer :: strainname
    integer :: jpop
    integer :: jstrain
    ! Time is the time of the last ancestor described; it starts at 0
    activepop = npop
    call canonical (npop, nu, numstrain)
    ! activepop is the number of active populations (note some will disappear
    ! in the backwards simulation as they become invented in reverse)
    time = 0.0
    ntotalpop = nu
    ntotalpopwithinvented = nu
    ! here we assign strains to the various ecotypes
    strainname = 0
    do jpop = 1, npop
      realname(jpop) = jpop
      activename(jpop) = jpop
      ! the realname of a population is its original number;
      ! the activename is the name from the list of still-active populations
      do jstrain = 1, numstrain(jpop)
        strainname = strainname + 1
        namestrain(jpop, jstrain) = strainname
        ! here we set up each strain as its own ancestor;
        ! that is, the first ancestor of each strain is itself
        numanc(strainname) = 1
        nameanc(strainname, 1) = strainname
        ! numanc(i) is the number of ancestors for the ith strain
        ! nameanc(i,j) is the name of the jth ancestor of the ith strain
        ! each ancestor has all its descendants listed, as below
        ancdist(strainname, 1) = 0.0
        ! ancdist(i,k) is the actual divergence between contemporary strain i
        ! and its kth ancestor.  Note, the first ancestor of each strain is
        ! itself.
        numdesc(strainname) = 1
        ! numdesc (i) is the number of descendants of the ith ancestor
        ! note that the 1 through nu terminal nodes have their first
        ! ancestor as 1 through nu, respectively
        namedesc(strainname, 1) = strainname
      end do
    end do
    numanctot = nu
    ! numanctot is the total number of ancestors, over all strains
    return
  end subroutine startpops

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Test the data for a successful fit at the given success level.
  !
  !  Params:
  !    nclustersarray        : The number of clusters (bins) at identity level sought.
  !    numcrit               : The number of criteria for making cluster bins.
  !    realdata              : The number of bins in the real data at criterion (1..numcrit).
  !    successlevel          : Success level to use.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function testForSuccessFit (nclustersarray, numcrit, realdata, successlevel)
    integer, intent(in) :: nclustersarray(:)
    integer, intent(in) :: numcrit
    integer, intent(in) :: realdata(:)
    real, intent(in)    :: successlevel
    ! Local variables.
    integer :: jcrit
    logical :: success
    real    :: xbin
    real    :: xreal
    success = .true.
    do jcrit = 1, numcrit
      xreal = realdata(jcrit)
      xbin = nclustersarray(jcrit)
      if (xreal / xbin .gt. successlevel .or. xbin / xreal .gt. successlevel) then
        success = .false.
        exit
      end if
    end do
    testforsuccessfit = success
  end function testForSuccessFit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Returns the name of the key event and gives the time until the event.
  !
  !  Params:
  !    parameters            : The parameters for this simulation (omega, sigma, npop, xn).
  !    activepop             : The number of active populations.
  !    eligibleparent        :
  !    eligibleps            :
  !    event                 : The event.
  !    numstrain(i)          :
  !    realname(i)           : The realname of a population is its original number.
  !    time                  : The total time (in substitutions) before the present for a coalescence event.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine whichEventAndWhen (parameters, activepop, event, numstrain, realname, time)
    type(parameters_data), intent(in) :: parameters
    integer, intent(in)               :: activepop
    integer, intent(out)              :: event
    integer, intent(in)               :: numstrain(:)
    integer, intent(in)               :: realname(:)
    real, intent(out)                 :: time
    ! Local variables.
    integer :: jpop
    integer :: eligibleparent
    integer :: eligibleps
    real    :: effectiveOmega
    real    :: effectiveSigma
    real    :: effectiveDrift
    real    :: x
    real    :: timeWait
    real    :: rateKey
    ! subroutine eligible returns the number of populations eligible to be the parent ecotypes, as well as eligible for periodic selection
    call eligible (activepop, eligibleparent, eligibleps, numstrain, realname)
    ! note that omega is a per population rate of niche invasion,
    ! so we multiply by the number of eligible parent populations
    effectiveOmega = eligibleparent * parameters%omega
    effectiveSigma = eligibleps * parameters%sigma
    effectiveDrift = 0.0
    do jpop = 1, activepop
      effectiveDrift = effectiveDrift + (1.0 / parameters%xn) * numstrain(realname(jpop)) * (numstrain(realname(jpop)) - 1.0) * 0.5
    end do
    ! rateKey is the total rate of all key events.
    rateKey = effectiveOmega + effectiveSigma + effectiveDrift
    ! The following uses 0.01 / rateKey as the length of time over which the probability of getting a key event is 0.01.
    ! The full formula for timeWait is:  timeWait = (-1 * log (x) / 0.01) * (0.01 / rateKey)
    call random_number (x)
    timeWait = -1 * log (x) / rateKey
    ! time is incremented by the timeWait
    ! time is the total time (in substitutions) before the present for a coalescence event
    time = time + timeWait
    ! now, what kind of key event
    call random_number (x)
    if (x .lt. effectiveOmega / rateKey) then
      event = EVENT_NICHE_INVASION
    else if (x .lt. (effectiveOmega + effectiveSigma) / rateKey) then
      event = EVENT_PERIODIC_SELECTION
    else if (x .lt. (effectiveOmega + effectiveSigma + effectiveDrift) / rateKey) then
      event = EVENT_POPULATION_DRIFT
    endif
    return
  end subroutine whichEventAndWhen

end module Methods
