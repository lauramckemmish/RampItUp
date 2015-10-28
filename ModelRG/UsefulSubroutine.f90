real*8 function IRR(m,n)
 implicit none
 integer, intent(in) :: m,n
 real*8  :: pi
 real*8  :: Nhat

 pi = dacos(-1d0);

 IRR = Nhat(m)*Nhat(n)*(8d0*pi)/((1+m+n)*(2+m+n)*(3+m+n));

end function IRR


!*********

real*8 function NHat(m)
 implicit none
 integer, intent(in) :: m
 real*8 :: pi

 pi = dacos(-1d0)

 Nhat = sqrt((2*m+1d0)*(2*m+2d0)*(2*m+3d0)/(8d0*pi))


end function Nhat


!***************

real*8 function IR1(n)
 implicit none
 integer, intent(in) :: n
 real*8 :: den, pi, num 

 pi = dacos(-1d0);
 num = sqrt((2*n+1)*(2*n+2)*(2*n+3)*2d0*pi)
 den = (n+1d0)*(n+2d0)
 IR1 = num/den;

end function IR1

!***************

real*8 function IR2(n)
 implicit none
 integer, intent(in) :: n
 real*8 :: den, pi, num 

 pi = dacos(-1d0);
 num = sqrt((2*n+1)*(2*n+2)*(2*n+3)*8d0*pi)
 den = (n+1d0)*(n+2d0)*(n+3d0)
 IR2 = num/den;

end function IR2


!***************

real*8 function IR3(n)
 implicit none
 integer, intent(in) :: n
 real*8 :: den, pi, num 

 pi = dacos(-1d0);
 num = sqrt((2*n+1)*(2*n+2)*(2*n+3)*72d0*pi)
 den = (n+1d0)*(n+2d0)*(n+3d0)*(n+4d0)
 IR3 = num/den;

end function IR3

!************************



real*8 function gpf(alpha)
 implicit none
 real*8, intent(in) :: alpha
 real*8 :: pi

pi = dacos(-1d0);
gpf = (2*alpha/pi)**(0.75d0)

end function gpf

!********************


subroutine calcintRampGaus4pirK(intRampGaus4pirK ,n,alpha,K, pts,type)
  implicit none
  integer, intent(in)      :: pts, type
  integer                  :: i, n, k, v
  real*8 :: intRampGaus4pirK
  real*8                   :: absc(1:2*pts+1), weights(1:2*pts+1) 
  real*8                   :: alpha, integrand(1:2*pts+1), otherpart(1:2*pts+1)
  real*8  :: Nhat, gpf, pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates (S|s) for concentric Ss using an asymptotic expression
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getquad(weights,absc,pts,type);
pi = dacos(-1d0);
  do v=1,2*pts+1; 
       integrand(v) = 4d0*pi*absc(v)**K*weights(v);  
        otherpart(v) = absc(2*pts+1-v+1)**n*exp(-alpha*absc(v)**2);
     end do
     intRampGaus4pirK = sum(integrand*otherpart)*Nhat(n)*gpf(alpha); 

end subroutine calcintRampGaus4pirK





