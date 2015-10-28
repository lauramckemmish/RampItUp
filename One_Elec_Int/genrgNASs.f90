!************************************************************************************

subroutine genRgSsna(naRgSs,no,RgSs,noatoms,pts)
  implicit none
  integer, intent(in)      :: no, noatoms, pts
  real*8, intent(in)       :: RgSs(1:30,1:no)
  real*8, intent(out)      :: naRgSs(1:no,1:4,1:noatoms)
  real*8                   :: beta,  BAx, BAy, BAz, BA, BA2
  integer                  :: rampcenter, n, i, v
  real*8                   :: absc(1:2*pts+1), weights(1:2*pts+1)
  real*8                   :: integrand(1:2*pts+1), otherpart(1:2*pts+1)
  
  naRgSs = 0d0;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call getquad(weights,absc,pts); 

  do v=1,2*pts+1;   integrand(v) = 2d0*dacos(-1d0)*weights(v);   end do

    do i=1,no
      beta = RgSs(8,i)  
      n    = RgSs(7,i) 
      BAx = RGSs(4,i);   BAy = RGSs(5,i) ;   BAz = RGSs(6,i) 
      BA2 =(BAx**2+BAy**2+BAz**2); BA = sqrt(BA2);
      rampcenter = int(RgSs(3,i))
    if (BA2 .lt. 1d-10) then
     do v=1,2*pts+1;
       otherpart(v) = 2*exp(-beta*absc(v)**2)*absc(v)*absc(2*pts+1-v+1)**n
     end do
    else
     do v=1,2*pts+1;
        otherpart(v) = (1d0/(BA*beta))*absc(2*pts+1-v+1)**n*exp(-beta*(absc(v)**2+BA2))*sinh(2*beta*BA*absc(v));
     end do
    end if

      naRgSs(i,1:4,rampcenter) = (/sum(integrand*otherpart), RgSs(1,i), RgSs(2,i), RgSs(9,i)/);
    end do

  
  return
end subroutine genRgSsna
