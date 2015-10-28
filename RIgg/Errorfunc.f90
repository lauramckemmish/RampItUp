
        SUBROUTINE calcerrorfun(X,ERR,EXPX2,pts)

!       =========================================
!       Purpose: Compute error function erf(x)
!       Input:   x   --- Argument of erf(x)
!       Output:  ERR --- erf(x)
!       =========================================

 implicit none
 real*8 :: eps, pi, er, c0, x2, r
 integer :: i,pts, k
 real*8 :: x(1:2*pts+1),err(1:2*pts+1), expX2(1:2*pts+1)

        EPS=1.0D-15


        PI=3.141592653589793D0
 do i=1,2*pts+1   
     X2=X(i)*X(i)
        IF (DABS(X(i)).LT.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(R).LE.DABS(ER)*EPS) GO TO 15
10         CONTINUE
15         C0=2.0D0/DSQRT(PI)*X(i)*EXPX2(i);
           ERR(i)=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO 20 K=1,12
              R=-R*(K-0.5D0)/X2
20            ER=ER+R
           C0=EXPX2(i)/(DABS(X(i))*DSQRT(PI))
           ERR(i)=1.0D0-C0*ER
           IF (X(i).LT.0.0) ERR=-ERR
        ENDIF
 end do
        RETURN
        END

