  !************************************************************************************
  subroutine newgenRgSsol(OLrgSs,no,RGSs,pts)
    implicit none
    integer, intent(in)          :: no, pts
    real*8, intent(in)           :: RGSs(1:30,1:no)
    real*8, intent(out)          :: OLrgSs(1:no,1:4)
    integer                      :: n, i, v
    real*8                       :: beta,  BAx, BAy, BAz, BA2, BA
    real*8                       :: absc(1:2*pts+1), weights(1:2*pts+1), integrand(1:2*pts+1), otherpart(1:2*pts+1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program calculates (S|s) using COLD-like algorithm (truncation of infinite expansion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call getquad(weights,absc,pts); 

  do v=1,2*pts+1;   integrand(v) = 2d0*dacos(-1d0)*absc(v)*weights(v);   end do

    do i=1,no
      beta = RgSs(8,i)  
      n    = RgSs(7,i) 
      BAx = RGSs(4,i);   BAy = RGSs(5,i) ;   BAz = RGSs(6,i) 
      BA2 =(BAx**2+BAy**2+BAz**2); BA = sqrt(BA2);

if (BA2 .lt. 1d-10) then
    do v=1,2*pts+1;
       otherpart(v) = 2*exp(-beta*absc(v)**2)*absc(v)*absc(2*pts+1-v+1)**n
    end do
else
     do v=1,2*pts+1;
        otherpart(v) = (1d0/(BA*beta))*absc(2*pts+1-v+1)**n&
                 *exp(-beta*(absc(v)**2+BA2))*sinh(2*beta*BA*absc(v));
     end do
end if

      OLrgSs(i,1) = sum(integrand*otherpart);
      OLrgSs(i,2) = RGSs(1,i);
      OLrgSs(i,3) = RGSs(2,i);
      OLrgSs(i,4) = RGSs(9,i);
    end do

    return
  end subroutine newgenRgSsol
