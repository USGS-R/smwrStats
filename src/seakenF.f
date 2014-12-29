      subroutine seakenf(x,n,nseas,results)
c     Seasonal Kendall test for trend.
c
c     Arguments:
c        x is the data vector in chronological order with one observation
c           (or missing value = -99999.0) per time period (season or month).
c        n is the length (total number of time periods) of x, it should be
c           be a multiple of nseas.
c        nseas is the number of equally sized time periods per year.
c        results is the output vector containing
c           results(1) = Kendall's tau,
c           results(2) = p-level without serial correlation correction,
c           results(3) = p-level with serial correlation correction,
c           results(4) = slope estimate,
c           results(5) = median value of time, and
c           results(6) = median value of the data.
c           results(7) = S
c           results(8) = VarS
c           results(9) = CovarS
c        Each p-level is the attained (two-sided) significance level of the
c           test [e.g., if (p-level.le.0.05) reject the null hypothesis
c           of no trend at the 5% level].
c
c     Limitations:
c        Tere must be at least two years of data.
c
c     References:
c     1) A Nonparametric Trend Test for Seasonal Data with Serial
c          Dependence, R.M.Hirsch & J.R.Slack, Water Resources Research,
c          Vol.20, No.6, Pages 727-732, June 1984.
c     2) Techniques of Trend Analysis for Monthly Water Quality Data,
c          R.M.Hirsch, J.R.Slack, & R.A.Smith, Water Resources Research,
c          Vol.18, No.1, Pages 107-121, February 1982.
c
c     Variables names correspond approximately to the variables used in
c        the first reference. Equation numbers cited in the code comments
c        refer to the first reference.
c
c     Subprograms called:
c        kens to compute the basic statistics for a season.
c        rank to compute the ranks of the data.
c        sgn to compute the sign of a value.
c        cdfn to compute the cumulative Normal distribution function.
c        vssort to sort a vector
c
c     Code history:
c        28 May 85  JRSlack  Initial coding.
c        18 Dec 85  JRSlack  Trap to handle each season being all ties.
c        18 Oct 00  JRSlack  Converted to DLL for S-Plus by moving error
c                               checking to the S function.
c        20 Dec 02  JRSlack  Converted arrays to automatic to eliminate 
c                               upper limit on seasons and years.
c        27 Aug 07 DLLorenz  extended results to 9 to return S, varS, covS
c
      real x(*), results(*)
      real, dimension(n/nseas, nseas) :: xjg, rjg
      real, dimension(max(n,n*(n/nseas-1)/2)) :: xdiff
      real xcheck/-99998.0/   ! Check for missing values against this.
c
c       Compute the s value for each of the nseas seasons,
c          compute the diagonal elements of the covariance matrix,
c          and add them up.
c
      ny = n / nseas  ! Integral number of years.
      index = 1       ! Next available location in array xdiff.
      sumcomp = 0.0   ! Number of comparisons.
      sums = 0.0      ! Sum of S's.
      sumsiggh = 0.0  ! Sum of sigma sub gh.
      sumngh = 0.0    ! Sum for 3rd part of equation 14.
      sumngh2 = 0.0   ! Sum of squares for 3rd part of equation 14.
c
c        Move the data into a new array to make subroutine calls easier.
c
      do 20 iseas = 1, nseas
      ng = 0
      do 10 iy = 1 , ny   ! Each time period corresponding to this season.
      xjg(iy,iseas) = x((iy-1) * nseas + iseas)
      if (xjg(iy,iseas) .gt. xcheck) ng = ng + 1
   10 continue
      sumngh = sumngh + ng+1                ! Note addition of 1 for eq. 14.
      sumngh2 = sumngh2 + (ng+1)**2         !   "
      call kens(xjg(:,iseas),ny,s,var,ncomp,xdiff(index))
c
      index = index + ncomp
      sums = sums + s
      sumsiggh = sumsiggh + var             ! Variance part of eq. 5.
      sumcomp = sumcomp + ncomp
      call rank(xjg(:,iseas),rjg(:,iseas),ny)  ! Ranks for RigRih.
   20 continue
c
c     Check for the possibility that each season is all ties;i.e., was a
c        constant giving zero variance in each season.
c
      if (sumsiggh.le.0.0) then
         results(1)=0.0
         results(2)=1.0
         results(3)=1.0
         results(4)=0.0
         go to 70
      endif
c
c     First results.
c
      tau = sums / sumcomp                  ! Kendall's tau.
c
      results(1)=tau
c
      results(7)=sums
      if (sums.gt.0.0) sums = sums - 1.0    ! Continuity correction.
      if (sums.lt.0.0) sums = sums + 1.0    !   "
      results(8)=sumsiggh
      z = sums / sqrt(sumsiggh)             ! Approximately Normal deviate.
      if (z.le.0.0) plevel = 2.0 * cdfn(z)
      if (z.gt.0.0) plevel = 2.0 * (1.0 - cdfn(z))
c
c                                           ! p-level without correction.
c
                             results(2)=plevel
c
      sumngh=sumngh**2-sumngh2              ! Sum over g&h of (ng+1)(nh+1).
c
c     Now calculate the serial correlation correction.
c       Compute covariances.
c
      sumkgh=0.0
      sumrgh=0.0
      do 60 j=1,ny
c
c     Sum RigRih part.
c
      sumr=0.0
      sumr2=0.0
      do 30 iseas=1,nseas
      sumr=sumr+rjg(j,iseas)
      sumr2=sumr2+rjg(j,iseas)**2
   30 continue
      sumrgh=sumrgh+sumr**2-sumr2
      do 50 i=1,j-1   ! This requires zero-trip DO loop.
c
c     Sum Kgh part.
c
      sumx=0.0
      sumx2=0.0
      do 40 iseas=1,nseas
      xj=xjg(j,iseas)
      xi=xjg(i,iseas)
      if (xi.le.xcheck.or.xj.le.xcheck) go to 40  ! Skip missing values.
      s=sgn(xj-xi)
      sumx=sumx+s
      sumx2=sumx2+s**2
   40 continue
      sumkgh=sumkgh+sumx**2-sumx2
   50 continue
   60 continue
c
c     Add covariance to eq. 5.
c
      sumsiggh=sumsiggh+(sumkgh+4.0*sumrgh-ny*sumngh)/3.0
      z = sums / sqrt(sumsiggh)             ! Approximately normal deviate.
      results(9)=sumsiggh - results(8)
      if (z.le.0.0) plevel = 2.0 * cdfn(z)
      if (z.gt.0.0) plevel = 2.0 * (1.0 - cdfn(z))
c                                           ! p-level with correction.
c
                             results(3)=plevel
c
c     Calculte the slope estimate.
c
      index = index - 1                     ! Adjust index to actual used.
c
      call vssort(xdiff,index)               ! Sort the difference vector.
      slope = (xdiff((index+1)/2) + xdiff((index+2)/2)) / 2.0
c                                           ! The median is the slope.
c
                             results(4)=slope
c
c     Pick up here on zero variance seasons since the seasons may
c        each be constant but not the same constant.
c
   70 continue
c
c     The trend line is the line with given slope going thru the
c        point (tmed,dmed) where
c           tmed is the median in time (x-axis) and
c           dmed is the median of the data (y-axis).
c
      tmed = ny / 2.0                       ! Median vaule of time.
c
                             results(5)=tmed
c
      index = 0                             ! Last used location in xdiff.
      do 80 i = 1 , ny*nseas
      if (x(i).le.xcheck) go to 80          ! Skip if missing.
      index = index + 1                     ! Next location.
      xdiff(index) = x(i)                   ! Save good value.
   80 continue
c
      call vssort(xdiff,index)               ! Sort the data.
c                                           ! Median of the data.
      dmed = (xdiff((index+1)/2) + xdiff((index+2)/2)) / 2.0
c
                             results(6)=dmed
c
      return
c
c     debugging section.
C  400 write(6,9000) x
C 9000 format(1x,12f10.1)
C      write(6,9001) sumsiggh
C 9001 format(//1x,f15.2,'  sumsiggh')
C      plevel = 1.00
C      return
      end

      subroutine kens(x,n,s,var,ncomp,xdiff)
c
c     Computes Kendall's S statistic and other information for
c        Kendall's tau test for trend. Sample usage of this
c        routine may be seen in routine kentau.
c
c     Arguments:
c        x is the time series vector of equally spaced observation (set
c           missing values to -99999.0).
c        n is the number of observations in x.
c        s is returned as Kendall's S statistic.
c        var is returned as the variance of S adjusted for ties.
c        ncomp is returned as the number of comparisons made in
c           computing S.
c        xdiff is returned as an array containing the values of each
c           comparison adjusted for the distance apart.
c
c     Reference:
c     Rank Correlation Methods, M.G.Kendall, Charles Griffin,
c        London, 1975
c
c     Subprograms called:
c        None.
c
c     Code history:
c        28 May 85  JRSlack  Initial coding.
c        20 Sep 88  JRSlack  Clarify the calculation of a tie.
c        19 Dec 02  JRSlack  Removed the limit on n.
c
      real x(*),xdiff(*)
      logical, dimension(n) :: wastie ! Was a value previously in a tie?
      real xcheck/-99998.0/           ! Check for missimg values against this.
      ncomp = 0                       ! Set number of comparisons to zero.
      if(n.gt.1) go to 10
c
c         Case of n <= 1.
c
      s = 0.0           ! Nil S.
      var = 0.0         ! Nil variance.
      return
c
c        Case of n > 1.
c
   10 nplus = 0         ! Number of up ticks.
      nminus = 0        ! Number of down ticks.
      fixvar = 0.0      ! Variance correction for ties.
c
      do 20 i = 1, n
      wastie(i) = .false.   ! Clear wastie array.
   20 continue
c
      do 40 istart = 1, n-1              ! Pick an observation.
      if (x(istart).le.xcheck) go to 40   ! Skip missing values.
      ntie=1                             ! A value is always tied with itself.
c
      do 30 iend = istart+1, n           ! Test each later observation.
      if (x(iend).le.xcheck) go to 30     ! Skip missing values.
c
c         Valid pair. Compare them.
c
      ncomp = ncomp + 1                              ! Note comparison.
c     Explanation of the 20 Sep 88 change:
c      Any of the x values may be a calculated (rather than read in) value;
c      e.g., (0.1+0.7)/2.0.  Such a calculated value may not be stored
c      internally the same as its indicated value (0.4) read directly in.
c      So the arithmetic may not determine a tie when there is one.  The IF
c      statement below takes care of this.  Specifically, a tie is declared
c      if two numbers differ by less than 1.0E-5 (10 to the minus 5th) times
c      their average.
c      The following line is commented out
c     yy = (x(iend) - x(istart)) / (iend - istart)   ! Adjust for separation.
c      and replaced by the following six lines of code.
      yy = x(iend) - x(istart)                       ! Calculte the difference.
      ave = (x(iend) + x(istart)) / 2.0              ! Calculate the average.
      if (yy.ne.0.0 .and. ave.ne.0.0) then           ! Want non-zero values.
         if (abs(yy/ave) .lt. 1.0e-5) yy = 0.0       ! Zero spurious
      endif                                          ! differences.
      yy = yy / (iend - istart)                      ! Adjust for separation.
c      End of 20 Sep 88 change.
      if (yy.gt.0.0) nplus = nplus + 1               ! Up tick.
      if (yy.lt.0.0) nminus = nminus + 1             ! Down tick.
      if (yy.eq.0.0) ntie = ntie + 1                 ! Tie.
      if (yy.eq.0.0) wastie(iend) = .true.           ! Mark ties.
      xdiff(ncomp) = yy                              ! Save differences.
   30 continue
c
c        Update variance correction if tie occured and tie was not counted
c           before.
c
      if (ntie.ne.1.and..not.wastie(istart)) fixvar = fixvar +
     1   ntie * (ntie-1.0) * (2.0*ntie+5.0) / 18.0
   40 continue
c
c         Compute the variance of the test statistic S.
c
      nok = (1.0 + sqrt(1.0+8.0*ncomp)) / 2.0       ! Number of actual values.
      var = nok * (nok-1.0) * (2.0*nok+5.0) / 18.0  ! Simple variance.
      var = var - fixvar                            ! Adjust for ties.
c
c       Compute the test statistic S.
c
      s = nplus - nminus
      return
      end

      subroutine rank(x,r,n)
c
c       Computes the ranks of the values in vector x.
c       x is not re-arranged on return.
c       In cases of ties, mid-ranks are used.
c       Missing values are assigned overall mid-rank.
c
c     Code history:
c        23 Sep 88  JRSlack  Skip sort & adjustment if all missing values.
c        19 Dec 02  JRSlack  Removed the limit on n.
c
      real x(*), r(*)
      real, dimension(n) :: y
      integer, dimension(n) :: ord
      real xcheck/-99998.0/
c
c       Initialize ord and y.
c
      m = 0
      do 10 i = 1 , n
      if (x(i).le.xcheck) go to 10
      m = m + 1
      ord(m) = i
      y(m) = x(i)
   10 continue
c
c       Rearrange y in ascending order.
c
      if (m.ne.0) call vsordr(y,m,ord)
c
c       Assign overall mid-rank to all
c       to take care of missing values.
c       Rank is 1/2 if all missing values - 23 Sep 88.
c
      if (m.eq.n) go to 30
      rm = (m + 1) / 2.0
      do 20 i = 1 , n
      r(i) = rm
   20 continue
      if (m.eq.0) return       ! Done if all missing values - 23 Sep 88.
   30 continue
c
c       Initial ranking.
c
      do 40 i = 1 , m
      r(ord(i)) = i
   40 continue
c
c       Adjust ranks for ties.
c
      i = 1
c
c       Value always tied with itself.
c
      iend = i
   50 iend = iend + 1
      if (iend.le.m.and.y(iend).eq.y(i)) go to 50
c
c       Compute average rank.
c
      ave = (iend*(iend-1) - i*(i-1)) / (2.0*(iend-i))
      do 60 j = i , iend-1
      r(ord(j)) = ave
   60 continue
      i = iend
      if (i.lt.m) go to 50
      return
      end

      function cdfn(x)
c     cumulative distribution function for the
c        normal zero-one distribution.
c
c     Primary reference is: Abramowitz & Stegun,
c        NBS Handbook of Mathematical Functions, equation 26.2.19
c
      if (x) 10,20,30
c
c     Negative argument.
c
   10 continue
      if (x.lt.-6.0) go to 40
      t=-x
      cdfn=      0.5/(1.0+0.0498673470*t+0.0211410061*t**2
     1   +0.0032776263*t**3+0.380036e-4*t**4+0.488906e-4*t**5
     1   +0.53830e-5*t**6)**16
      return
c
c     Zero argument.
c
   20 continue
      cdfn=0.5
      return
c
c     Positive argument.
c
   30 continue
      if (x.gt.6.0) go to 50
      cdfn=1.0-0.5/(1.0+0.0498673470*x+0.0211410061*x**2
     1   +0.0032776263*x**3+0.380036e-4*x**4+0.488906e-4*x**5
     1   +0.53830e-5*x**6)**16
      return
c
c     Outside the range +-6 the approximation is useless.
c
   40 continue
      cdfn=0.0
      return
   50 continue
      cdfn=1.0
      return
      end

      function sgn(x)
c
c     Computes the sign of the argument.
c
      sgn=1.0
      if (x.eq.0.0) sgn=0.0
      if (x.lt.0.0) sgn=-1.0
      return
      end
