subroutine calcggsske(keggss,noggss,ggss)
 implicit none
  integer, intent(in) :: noggss
  real*8, intent(in)  :: ggss(1:30,1:noggss)
  real*8              :: keggss(1:noggss,1:4)
  real*8              :: pi3r2, zeta, alpha, beta, RAB2
  integer             :: i

 pi3r2 = dacos(-1d0)**(3d0/2d0);

  do i=1,noggss
     RAB2 = ggss(16,i);
     beta = ggss(15,i);
     zeta = ggss(8,i);
     alpha = zeta - beta;
     keggss(i,1) = alpha*beta/(zeta**(5d0/2d0))*(3-2*alpha*beta/zeta*RAB2)*pi3r2;
     keggss(i,2) = ggss(1,i)
     keggss(i,3) = ggss(2,i)
     keggss(i,4) = ggss(9,i)
  end do

end subroutine calcggsske

!***************************

subroutine calcggspke(keggsp,noggsp,ggsp)
 implicit none
  integer, intent(in) :: noggsp
  real*8, intent(in)  :: ggsp(1:30,1:noggsp)
  real*8              :: keggsp(1:3*noggsp,1:4)
  integer             :: i
  real*8              :: ABx, ABy, ABz, pi3r2, alpha, beta, zeta, RAB2

 pi3r2 = dacos(-1d0)**(3d0/2d0);

  do i=1,noggsp
     ABx = ggsp(11,i); ABy = ggsp(12,i); ABz = ggsp(13,i);
     RAB2 = ABx**2+ABy**2+ABz**2;
     alpha = ggsp(6,i); beta = ggsp(7,i); zeta = ggsp(8,i);
     keggsp(3*i-2    ,1) = -alpha**2*beta/(zeta**(9d0/2d0))*pi3r2*ABx*&
        (2*RAB2*alpha*beta-5*zeta); 
     keggsp(3*i-1    ,1) = -alpha**2*beta/(zeta**(9d0/2d0))*pi3r2*ABy*&
        (2*RAB2*alpha*beta-5*zeta); 
     keggsp(3*i      ,1) = -alpha**2*beta/(zeta**(9d0/2d0))*pi3r2*ABz*&
        (2*RAB2*alpha*beta-5*zeta); 
     keggsp(3*i-2:3*i,2) = ggsp(1,i)
     keggsp(3*i-2    ,3) = ggsp(2,i)
     keggsp(3*i-1    ,3) = ggsp(2,i)+1
     keggsp(3*i      ,3) = ggsp(2,i)+2
     keggsp(3*i-2:3*i,4) = ggsp(9,i)
  end do

end subroutine calcggspke


!***************************

subroutine calcggppke(keggpp,noggpp,ggpp)
 implicit none
  integer, intent(in) :: noggpp
  real*8, intent(in)  :: ggpp(1:30,1:noggpp)
  real*8              :: keggpp(1:9*noggpp,1:4)
  integer             :: i
  real*8              :: ABx, ABy, ABz, pi3r2, alpha, beta, zeta, RAB2

 pi3r2 = dacos(-1d0)**(3d0/2d0);

  do i=1,noggpp
     ABx = ggpp(11,i); ABy = ggpp(12,i); ABz = ggpp(13,i);
     RAB2 = ABx**2+ABy**2+ABz**2;
     alpha = ggpp(6,i); beta = ggpp(7,i); zeta = ggpp(8,i);

     keggpp(9*i-8      ,1) = 0.5d0*pi3r2*alpha*beta/(zeta**(11d0/2d0))*&
           (4*ABx**2*RAB2*alpha**2*beta**2-2*alpha*beta*zeta*(7*ABx**2+RAB2)+5*zeta**2); 
     keggpp(9*i-7      ,1) = ABx*ABy*pi3r2*alpha**2*beta**2*(2*RAB2*alpha*beta-7*zeta)/(zeta**(11d0/2d0));
     keggpp(9*i-6      ,1) = ABx*ABz*pi3r2*alpha**2*beta**2*(2*RAB2*alpha*beta-7*zeta)/(zeta**(11d0/2d0));
     keggpp(9*i-5      ,1) = ABx*ABy*pi3r2*alpha**2*beta**2*(2*RAB2*alpha*beta-7*zeta)/(zeta**(11d0/2d0));
     keggpp(9*i-4      ,1) = 0.5d0*pi3r2*alpha*beta/(zeta**(11d0/2d0))*&
           (4*ABy**2*RAB2*alpha**2*beta**2-2*alpha*beta*zeta*(7*ABy**2+RAB2)+5*zeta**2); 
     keggpp(9*i-3      ,1) = ABz*ABy*pi3r2*alpha**2*beta**2*(2*RAB2*alpha*beta-7*zeta)/(zeta**(11d0/2d0));
     keggpp(9*i-2      ,1) = ABx*ABz*pi3r2*alpha**2*beta**2*(2*RAB2*alpha*beta-7*zeta)/(zeta**(11d0/2d0));
     keggpp(9*i-1      ,1) = ABz*ABy*pi3r2*alpha**2*beta**2*(2*RAB2*alpha*beta-7*zeta)/(zeta**(11d0/2d0));
     keggpp(9*i-0      ,1) = 0.5d0*pi3r2*alpha*beta/(zeta**(11d0/2d0))*&
           (4*ABz**2*RAB2*alpha**2*beta**2-2*alpha*beta*zeta*(7*ABz**2+RAB2)+5*zeta**2); 
     keggpp(9*i-8:9*i-6,2) = ggpp(1,i)
     keggpp(9*i-5:9*i-3,2) = ggpp(1,i)+1
     keggpp(9*i-2:9*i-0,2) = ggpp(1,i)+2
     keggpp(9*i-8      ,3) = ggpp(2,i)
     keggpp(9*i-7      ,3) = ggpp(2,i)+1
     keggpp(9*i-6      ,3) = ggpp(2,i)+2
     keggpp(9*i-5      ,3) = ggpp(2,i)
     keggpp(9*i-4      ,3) = ggpp(2,i)+1
     keggpp(9*i-3      ,3) = ggpp(2,i)+2
     keggpp(9*i-2      ,3) = ggpp(2,i)
     keggpp(9*i-1      ,3) = ggpp(2,i)+1
     keggpp(9*i        ,3) = ggpp(2,i)+2
     keggpp(9*i-8:9*i  ,4) = ggpp(9,i)
  end do


end subroutine calcggppke
