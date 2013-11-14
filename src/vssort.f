************************************************************************
*
*     Subroutine VSSORT                           Called by: seakencode 
*                VSORDR
*
*     Fast routine for sorting a double precision vector in ascending
*     order. Based on Algorithm 271, Quickersort, (R.S. Scowen,
*     collected algorithms of the ACM).
*
*     subroutine arguments:
*
*     A   - vector to be sorted (on input); sorted vector (on output)
*     N   - number of values in vector
*         vsordr only:
*     IA  - integer vector recording the order
*
*     The method continually splits the vector into parts such that all
*     elements of one part are less than all elements of the other, with
*     a third part in the middle consisting of one element.  An element
*     with value T is chosen arbitrarily (here we choose the middle
*     element). I and J give the lower and upper limits of the segment
*     being split. After the split a value Q will have been found such
*     that A(Q)=T and A(L)<=T<=A(M) for all I<=L<Q<M<=J. The routine
*     then performs operations on the two segments (I,Q-1) and (Q+1,J)
*     as follows. The smaller segment is split and the position of the
*     larger segment is stored in the LT and UT vectors. If the segment
*     to be split contains two or fewer elements, it is sorted and
*     another segment is obtained from the LT and UT vectors. When no
*     more segments remain, the vector is completely sorted.
*
*     Vectors can be sorted in descending order by changing some of 
*     the LE and LT tests - see original code.
*
*     local vars
*     ----------
*
************************************************************************
      SUBROUTINE VSSORT(A,N)
*     
*     subroutine arguments
*
      real A(*)
      INTEGER*4 N
*
*     local vars (the dimensions for LT and UT have to be at least
*     log-base2 N. 17 was chosen to handle N<131,073.)
*
      real T,X
      INTEGER LT(17),UT(17),I,J,K,M,P,Q
*
*     initialize
*
      M = 1
      I = 1
      J = N
*
*     if this segment has more than two elements, split it
*
 10   IF ((J-I-1).LT.0) THEN
         GOTO 100
      ELSEIF ((J-I-1).EQ.0) THEN
         GOTO 90
      ELSEIF ((J-I-1).GT.0) THEN
         GOTO 15
      ENDIF
*
*     P is the position of an arbitrary element in the segment we choose
*     the middle element. Under certain circumstances it may be
*     advantageous to choose P at random.
*
 15   P = (J+I)/2
      T = A(P)
      A(P) = A(I)
*
*     starting at the beginning of the segment, search for K such that
*     A(K)>T
*
      Q = J
      K = I
 20   K = K+1
      IF (K.GT.Q) GOTO 60
      IF (A(K).LE.T) GOTO 20
*
*     such an element has now been found.  Now search for a Q such that
*     A(Q)<T starting at the end of the segment.
*
 30   IF (A(Q).LT.T) GOTO 40
      Q = Q-1
      IF (Q.GT.K) GOTO 30
      GOTO 50
*
*     A(Q) has now been found -- interchange A(Q) and A(K)
*
 40   X = A(K)
      A(K) = A(Q)
      A(Q) = X
*
*     update Q and search for another pair to interchange
*
      Q = Q-1
      GOTO 20
 50   Q = K-1
 60   CONTINUE
*
*     the upwards search has now met the downwards search
*
      A(I) = A(Q)
      A(Q) = T
*
*     the segment is now divided in three parts: (I,Q-1),(Q),(Q+1,J).
*     Store the position of the largest segment in LT and UT
*
      IF (2*Q.LE.I+J) GOTO 70
      LT(M) = I
      UT(M) = Q-1
      I = Q+1
      GOTO 80
 70   LT(M) = Q+1
      UT(M) = J
      J = Q-1
*
*     update M and split the new smaller segment
*
 80   M = M+1
      GOTO 10
*
*     arrive here if the segment has 2 elements.  Test to see if the
*     segment is properly ordered.  If not, perform an interchange
*
 90   IF (A(I).LE.A(J)) GOTO 100
      X = A(I)
      A(I) = A(J)
      A(J) = X
*
*     if LT and UT contain more segments to be sorted repeat process
*
 100  M = M-1
      IF (M.LE.0) RETURN
      I = LT(M)
      J = UT(M)
      GOTO 10
      END

************************************************************************
      SUBROUTINE VSORDR(A,N,IA)
*     
*     subroutine arguments
*
      real A(*)
      INTEGER*4 N, IA(*)
*
*     local vars (the dimensions for LT and UT have to be at least
*     log-base2 N. 17 was chosen to handle N<131,073.)
*
      real T,X
      INTEGER LT(17),UT(17),I,J,K,M,P,Q, IT,IX
*
*     initialize
*
      M = 1
      I = 1
      J = N
*
*     if this segment has more than two elements, split it
*
 10   IF ((J-I-1).LT.0) THEN
         GOTO 100
      ELSEIF ((J-I-1).EQ.0) THEN
         GOTO 90
      ELSEIF ((J-I-1).GT.0) THEN
         GOTO 15
      ENDIF
*
*     P is the position of an arbitrary element in the segment we choose
*     the middle element. Under certain circumstances it may be
*     advantageous to choose P at random.
*
 15   P = (J+I)/2
      T = A(P)
      A(P) = A(I)
      IT = IA(P)
      IA(P) = IA(I)
*
*     starting at the beginning of the segment, search for K such that
*     A(K)>T
*
      Q = J
      K = I
 20   K = K+1
      IF (K.GT.Q) GOTO 60
      IF (A(K).LE.T) GOTO 20
*
*     such an element has now been found.  Now search for a Q such that
*     A(Q)<T starting at the end of the segment.
*
 30   IF (A(Q).LT.T) GOTO 40
      Q = Q-1
      IF (Q.GT.K) GOTO 30
      GOTO 50
*
*     A(Q) has now been found -- interchange A(Q) and A(K), and IA
*
 40   X = A(K)
      A(K) = A(Q)
      A(Q) = X
      IX = IA(K)
      IA(K) = IA(Q)
      IA(Q) = IX
*
*     update Q and search for another pair to interchange
*
      Q = Q-1
      GOTO 20
 50   Q = K-1
 60   CONTINUE
*
*     the upwards search has now met the downwards search
*
      A(I) = A(Q)
      A(Q) = T
      IA(I) = IA(Q)
      IA(Q) = IT
*
*     the segment is now divided in three parts: (I,Q-1),(Q),(Q+1,J).
*     Store the position of the largest segment in LT and UT
*
      IF (2*Q.LE.I+J) GOTO 70
      LT(M) = I
      UT(M) = Q-1
      I = Q+1
      GOTO 80
 70   LT(M) = Q+1
      UT(M) = J
      J = Q-1
*
*     update M and split the new smaller segment
*
 80   M = M+1
      GOTO 10
*
*     arrive here if the segment has 2 elements.  Test to see if the
*     segment is properly ordered.  If not, perform an interchange
*
 90   IF (A(I).LE.A(J)) GOTO 100
      X = A(I)
      A(I) = A(J)
      A(J) = X
      IX = IA(I)
      IA(I) = IA(J)
      IA(J) = IX
*
*     if LT and UT contain more segments to be sorted repeat process
*
 100  M = M-1
      IF (M.LE.0) RETURN
      I = LT(M)
      J = UT(M)
      GOTO 10
      END
