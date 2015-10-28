!************************************************************************************

subroutine genrgSske(kergSs,no,RGSs,pts)
  implicit none
  integer, intent(in)      :: no, pts
  real*8, intent(in)       :: RGSs(1:30,1:no)
  real*8, intent(out)      :: kergSs(1:no,1:4)
  real*8                   :: beta, BAx, BAy, BAz, BA2, BA
  integer                  :: n, i, v
  real*8                   :: absc(1:2*pts+1), weights(1:2*pts+1) 
  real*8                   :: integrand(1:2*pts+1), otherpart(1:2*pts+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine calculates the (S|s) kinetic energy using quadrature
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  kergSs = 0d0;

  call getquad(weights,absc,pts); 

  do v=1,2*pts+1;   integrand(v) = -dacos(-1d0)*weights(v);   end do

  do i=1,no

      beta = RgSs(8,i)  
      n    = RgSs(7,i) 
      BAx = RGSs(4,i);   BAy = RGSs(5,i) ;   BAz = RGSs(6,i) 
      BA2 =(BAx**2+BAy**2+BAz**2); BA = sqrt(BA2);

     kergSs(i,2) = RGSs(1,i)
     kergSs(i,3) = RGSs(2,i)
     kergSs(i,4) = RGSs(9,i)

if (BA2 .lt. 1d-10) then
    do v=1,2*pts+1;
       otherpart(v) = 2*exp(-beta*absc(v)**2)*absc(v)*n*absc(2*pts+1-v+1)**(n-2)*(-2+absc(v)*(n+1));
    end do
else
     do v=1,2*pts+1;
        otherpart(v) = (1d0/(BA*beta))*absc(2*pts+1-v+1)**(n-2)*(-2+absc(v)*(n+1))*n   &
                 *exp(-beta*(absc(v)**2+BA2))*sinh(2*beta*BA*absc(v));
     end do
end if

     kergSs(i,1) = sum(integrand*otherpart); 
  end do            
  return
end subroutine genrgSske

