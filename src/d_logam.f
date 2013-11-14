*     Function D_LOGAM                       Called by: PHI_MVUE
*
*
*
C     ALGORITHM ACM 291, COMM. ACM. (1966) VOL. 9, P. 686
C
C     EVALUATES NATURAL LOGARITHM OF GAMMA(X) FOR X GREATER THAN ZERO
*     This code does not test for (X .LE. 0.D0) as it shouldn't ever
*     happen
C
C     FORTRAN TRANSLATION PRESENTED IN PIKE, C.M., AND HILL, I.D., 1985,
C     LOGARITHM OF THE GAMMA FUNCTION, IN GRIFFITHS, P., AND HILL, I.D.,
C     1985, APPLIED STATISTICS ALGORITHMS: LONDON, ROYAL STATISTICAL
C     SOCIETY
*
*
**********************************************************************
      DOUBLE PRECISION FUNCTION D_LOGAM(X)
*
*     subroutine arg
*
      DOUBLE PRECISION X
*
*     local vars
*
      DOUBLE PRECISION A1,A2,A3,A4,A5,F,Y,Z,ZLOG
*
*     constants for DLOG(2 PI)/2, 1/1680, 1/1260, 1/360, AND 1/12
*
      DATA A1,A2,A3,A4,A5 / 0.918938533204673D0, 0.000595238095238D0,
     &                      0.000793650793651D0, 0.002777777777778D0,
     &                      0.083333333333333D0 /

      ZLOG(F) = DLOG(F)

      Y = X
      F = 0.D0
      IF (Y.LT.7.D0) THEN
         F = Y
 10      Y = Y+1.D0
         IF (Y.LT.7.D0) THEN
            F = F*Y
            GOTO 10
         ENDIF
         F = -ZLOG(F)
      ENDIF
      Z = 1.D0/(Y*Y)

      D_LOGAM = F+(Y-0.5D0)*ZLOG(Y)-Y+A1 +(((-A2*Z+A3)*Z-A4)*Z+A5)/Y

      RETURN
      END
