*	These functions are renamed to avoid any conflict with the same 
*	functions in the USGSqw package
*
      SUBROUTINE PHIMVUE(AW4, NU2, BCF, N)
      double precision aw4(N), nu2, bcf(N)
      double precision mvue_phi
      integer*4 N, i
      do 10 i=1,N
         bcf(i) = mvue_phi(aw4(i), nu2)
 10   continue
      return
      end
**********************************************************************
*
*     Function MVUE_PHI                            Called by: PHIMVUE
*
*     calculate the function PHI given in Likes (1980)
*
*     Arguments
*     ---------
*     AW4     aW/4 - a is from Likes and W is sum of squares
*     NU2     degrees of freedom/2
*
*     local vars
*     ----------
*     ARG1     DLOGAM(K+1)
*     ARG2     DLOGAM(NU2)
*     ARG3     DLOGAM(NU2+K)
*     IK       I-1
*     RK       DBLE(I-1)
*     TOL      convergence criteria (convergence when term < tol)
*     WORK     amount of MVUEPHI contributed by a single term
*
**********************************************************************
      DOUBLE PRECISION FUNCTION MVUE_PHI(AW4,NU2)
*
*     subroutine arguments
*
      DOUBLE PRECISION AW4,NU2
*
*     local variables
*
      INTEGER*4 I,IK
      DOUBLE PRECISION ARG1,ARG2,ARG3,RK,TOL,WORK
*
*     function declaration
*
      DOUBLE PRECISION D_LOGAM
*
*     define criteria for convergence
*
      PARAMETER (TOL=1.D-9)
*
*     for AW4 equal to zero, PHI equals 1
*
      IF (AW4 .EQ. 0.D0) THEN
         MVUE_PHI = 1.D0
         RETURN
      ENDIF
*
*     For AW4 not equal to zero, calculate PHI by iteration.  PHI may
*     be computed using factorial and gamma functions, e.g.:
*
*     WORK = (1.D0/FAC(RK))*(GAMMA(NU2)/GAMMA(NU2+RK))*GAMMA(AW4)**IK
*
*     instead of GAMMA(X) use EXP(LN_GAMMA(X)) to avoid machine
*     overflow.  This allows log GAMMA terms to be added/subtracted
*     rather than multiplying/dividing.  Use EXP(LOGAM(K+1)) for
*     FACTORIAL(K).
*
      MVUE_PHI = 0.D0
      DO 10 I=1,200
         IK = I-1
         RK = DBLE(I-1)
         ARG1 = D_LOGAM((RK+1.D0))
         ARG2 = D_LOGAM(NU2)
         ARG3 = D_LOGAM(NU2+RK)
         IF (AW4**IK.GT.2.D0**1000) THEN
            MVUE_PHI = -1.d+0
            RETURN
         ENDIF
         WORK = DEXP(-ARG1+ARG2-ARG3)*(AW4**IK)
         MVUE_PHI = MVUE_PHI+WORK
         IF (I.GT.2.AND.DABS(WORK).LT.TOL) RETURN
 10   CONTINUE
*
*     the function has not converged in 200 iterations; terminate
*     execution
*
      MVUE_PHI = -1.d+0
      RETURN
      END
