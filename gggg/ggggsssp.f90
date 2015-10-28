subroutine calcggggsssp(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggss,ggsp,noggss,noggsp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
 implicit none
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun),ggss(1:20,1:noggss),ggsp(1:20,1:noggsp)
  real*8, intent(in)      :: basis(1:nobasisfun,1:30)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer  :: i, j, bf1, bf2, bf3, bf4, printopt, i1, out,compare, lengthtwoelec, noggss, noggsp, nobasisfun
  real*8   :: alpha, beta, gamma, delta, zeta, eta, pf, f0eval,RPQ2
  real*8   :: tempint(1:3), u, uket2(1:3), expeval, uketint, nc1, nc2, CD(1:3), PQ(1:3)
  real*8   :: normtwoelecmat2(0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  real*8   :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  character(len=30), intent(in)           :: moleculename

#include "DEBUG.h"

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
call readQchemtwoeecompare(normtwoelecmat2,nobasisfun,basis,moleculename);
print *, "in ggggsssp";   compare = 0;
else;  compare = 0;  end if
#endif

  pf = (dacos(-1d0))**(5d0/2d0);
  do i=1,noggss
    alpha = ggss(14,i); beta = ggss(15,i); zeta = ggss(8,i); 
    bf1 = int(ggss(1,i)); bf2 = int(ggss(2,i));    nc1 = ggss(9,i);
    do j=1,noggsp
      gamma = ggsp(6,j); delta = ggsp(7,j); eta = ggsp(8,j);
      CD = ggsp(11:13,j);
      bf3 = int(ggsp(1,j)); bf4 = int(ggsp(2,j));      nc2 = ggsp(9,j);
      PQ = ggss(4:6,i)-ggsp(3:5,j);      RPQ2 = (PQ(1))**2+(PQ(2))**2+(PQ(3))**2;

      call calcf0(f0eval,zeta*eta/(zeta+eta)*RPQ2);      expeval = exp(-RPQ2*zeta*eta/(zeta+eta));

 if (RPQ2 .lt. 1d-10) then; 
   tempint = 2*CD*pf*gamma/(zeta*eta**2*sqrt(zeta+eta));
 else
   u = 2*pf/(zeta*eta*sqrt(zeta+eta))*f0eval;
   uketint = 2*pf/(RPQ2*zeta*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);  uket2 = uketint*PQ*delta; 
   tempint = CD*gamma*u/(eta)+uket2/(2d0*delta)
 end if
printopt = 0;  compare=0;
tempint = tempint*nc1*nc2

do i1=1,3;
 call digestintsnewcompare(tempint(i1), bf1  , bf2  , bf3  , bf4+i1-1  , Galpha,Gbeta,Ptot,Palpha,Pbeta, &
  nobasisfun, printopt,2, normtwoelecmat2,out,compare)
#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+i1-1, tempint(i1),2)
end if
#endif
end do
   end do
 end do



#ifdef INTPRINT
if (abs(Palpha(1,1)) .lt. 1d-10) then
!print *, "TESTING DIFFERENCES"
do bf1=1,nobasisfun; do bf2=bf1,nobasisfun; do bf3=bf1,nobasisfun; do bf4=bf2,nobasisfun;
   if (abs(newtwo(bf1,bf2,bf3,bf4,1)) .gt. 1d-10) then;
 if ( abs(newtwo(bf1,bf2,bf3,bf4,1)- normtwoelecmat2(bf1,bf2,bf3,bf4,1)) .gt. 1d-9) then
!     print *, bf1,bf2,bf3,bf4, newtwo(bf1,bf2,bf3,bf4,1), normtwoelecmat2(bf1,bf2,bf3,bf4,1)
 end if
   end if 
end do; end do; end do; end do;
!print *, "END OF TESTING DIFFERENCES"
end if
#endif

end subroutine calcggggsssp
