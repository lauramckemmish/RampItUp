subroutine calcggssol(OLggss,noggss,ggss)
 implicit none
  integer, intent(in) :: noggss
  real*8, intent(in)  :: ggss(1:30,1:noggss)
  real*8              :: olggss(1:noggss,1:4)
  real*8              :: pi3r2, zeta, alpha, beta, RAB2
  integer             :: i

 pi3r2 = dacos(-1d0)**(3d0/2d0);

  do i=1,noggss
     zeta = ggss(8,i)
     OLggss(i,1) = pi3r2/(zeta**(3d0/2d0)); 
     OLggss(i,2) = ggss(1,i)
     OLggss(i,3) = ggss(2,i)
     OLggss(i,4) = ggss(9,i)
  end do

end subroutine calcggssol

!***************************

subroutine calcggspol(OLggsp,noggsp,ggsp)
 implicit none
  integer, intent(in) :: noggsp
  real*8, intent(in)  :: ggsp(1:30,1:noggsp)
  real*8              :: olggsp(1:3*noggsp,1:4)
  integer             :: i
  real*8              :: ABx, ABy, ABz, pi3r2, alpha, zeta

 pi3r2 = dacos(-1d0)**(3d0/2d0);

  do i=1,noggsp
     ABx = ggsp(11,i); ABy = ggsp(12,i); ABz = ggsp(13,i);
     alpha = ggsp(6,i); zeta = ggsp(8,i);
     OLggsp(3*i-2    ,1) = ABx*pi3r2*alpha/(zeta**(5d0/2d0)); 
     OLggsp(3*i-1    ,1) = ABy*pi3r2*alpha/(zeta**(5d0/2d0)); 
     OLggsp(3*i      ,1) = ABz*pi3r2*alpha/(zeta**(5d0/2d0)); 
     OLggsp(3*i-2:3*i,2) = ggsp(1,i)
     OLggsp(3*i-2    ,3) = ggsp(2,i)
     OLggsp(3*i-1    ,3) = ggsp(2,i)+1
     OLggsp(3*i      ,3) = ggsp(2,i)+2
     OLggsp(3*i-2:3*i,4) = ggsp(9,i)
  end do

end subroutine calcggspol


!***************************

subroutine calcggppol(OLggpp,noggpp,ggpp)
 implicit none
  integer, intent(in) :: noggpp
  real*8, intent(in)  :: ggpp(1:30,1:noggpp)
  real*8              :: olggpp(1:9*noggpp,1:4)
  integer             :: i
  real*8              :: ABx, ABy, ABz, pi3r2, alpha, beta, zeta

 pi3r2 = dacos(-1d0)**(3d0/2d0);

  do i=1,noggpp
     ABx = ggpp(11,i); ABy = ggpp(12,i); ABz = ggpp(13,i);
     alpha = ggpp(6,i); beta = ggpp(7,i); zeta = ggpp(8,i);
     OLggpp(9*i-8      ,1) = pi3r2*0.5d0/(zeta**(7d0/2d0))*(alpha+beta-2*ABx**2*alpha*beta);
     OLggpp(9*i-7      ,1) = -pi3r2*ABx*ABy/(zeta**(7d0/2d0))*alpha*beta
     OLggpp(9*i-6      ,1) = -pi3r2*ABx*ABz/(zeta**(7d0/2d0))*alpha*beta
     OLggpp(9*i-5      ,1) = -pi3r2*ABx*ABy/(zeta**(7d0/2d0))*alpha*beta
     OLggpp(9*i-4      ,1) = pi3r2*0.5d0/(zeta**(7d0/2d0))*(alpha+beta-2*ABy**2*alpha*beta);
     OLggpp(9*i-3      ,1) = -pi3r2*ABy*ABz/(zeta**(7d0/2d0))*alpha*beta
     OLggpp(9*i-2      ,1) = -pi3r2*ABx*ABz/(zeta**(7d0/2d0))*alpha*beta
     OLggpp(9*i-1      ,1) = -pi3r2*ABz*ABy/(zeta**(7d0/2d0))*alpha*beta
     OLggpp(9*i-0      ,1) = pi3r2*0.5d0/(zeta**(7d0/2d0))*(alpha+beta-2*ABz**2*alpha*beta); 

     OLggpp(9*i-8:9*i-6,2) = ggpp(1,i)
     OLggpp(9*i-5:9*i-3,2) = ggpp(1,i)+1
     OLggpp(9*i-2:9*i-0,2) = ggpp(1,i)+2
     OLggpp(9*i-8      ,3) = ggpp(2,i)
     OLggpp(9*i-7      ,3) = ggpp(2,i)+1
     OLggpp(9*i-6      ,3) = ggpp(2,i)+2
     OLggpp(9*i-5      ,3) = ggpp(2,i)
     OLggpp(9*i-4      ,3) = ggpp(2,i)+1
     OLggpp(9*i-3      ,3) = ggpp(2,i)+2
     OLggpp(9*i-2      ,3) = ggpp(2,i)
     OLggpp(9*i-1      ,3) = ggpp(2,i)+1
     OLggpp(9*i        ,3) = ggpp(2,i)+2
     OLggpp(9*i-8:9*i  ,4) = ggpp(9,i)
  end do

end subroutine calcggppol
