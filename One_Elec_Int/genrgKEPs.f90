  !************************************************************************************
  
  subroutine calcgenRgPsKE(KErgPs,no,RGPs,pts)
    implicit none
    integer, intent(in)          :: no, pts
    real*8, intent(in)           :: RGPs(1:30,1:no)
    real*8, intent(out)          :: KErgPs(1:3*no,1:4)
    real*8                       :: beta,  BAx, BAy, BAz, BA2, BA
    integer                      :: n, i, v
    real*8                       :: absc(1:2*pts+1), weights(1:2*pts+1), pi, r, base
    real*8                       :: integrand(1:2*pts+1), otherpartx(1:2*pts+1)
    real*8                       :: otherparty(1:2*pts+1), otherpartz(1:2*pts+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
pi = dacos(-1d0)
  call getquad(weights,absc,pts); 


    do i=1,no
      beta = RgPs(8,i)  
      n    = RgPs(7,i) 
      BAx = RGPs(4,i);   BAy = RGPs(5,i) ;   BAz = RGPs(6,i) 
      BA2 =(BAx**2+BAy**2+BAz**2); BA = sqrt(BA2);

if (BA2 .lt. 1d-10) then
 otherpartx = 0d0; 
 otherparty = 0d0; 
 otherpartz = 0d0; 

else
     do v=1,2*pts+1;
        r = absc(v)

  base = exp(-beta*(r+BA)**2)*(1+2*r*BA*beta) + exp(-beta*(r-BA)**2)*(-1+2*r*BA*beta)
  base = base*(-n*sqrt(3d0*pi))/(8d0*BA2*beta**2)
  base = base*(1-r)**(n-2)*(-4+(3+n)*r)

  otherpartx(v) = base*BAx/BA;
  otherparty(v) = base*BAy/BA;
  otherpartz(v) = base*BAz/BA;
     end do
end if

      KErgPs(3*i-2,1) = sum(weights*otherpartx);
      KErgPs(3*i-1,1) = sum(weights*otherparty);
      KErgPs(3*i  ,1) = sum(weights*otherpartz);

      KErgPs(3*i-2      ,2) = RGPs(1,i);
      KErgPs(3*i-1      ,2) = RGPs(1,i)+1;
      KErgPs(3*i        ,2) = RGPs(1,i)+2;
      KErgPs(3*i-2:3*i  ,3) = RGPs(2,i);
      KErgPs(3*i-2:3*i  ,4) = RGPs(9,i);

    end do

    return
  end subroutine calcgenRgPsKE


