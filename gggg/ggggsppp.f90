subroutine calcggggsppp(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggsp,ggpp,noggsp,noggpp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
 implicit none
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: ggsp(1:20,1:noggsp),ggpp(1:20,1:noggpp), basis(1:nobasisfun,1:30)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: noggsp, noggpp, nobasisfun
  integer :: out,compare, i, j, bf1, bf2, bf3, bf4, printopt, i1, i2, i3, lengthtwoelec
  real*8 :: alpha, beta, gamma, delta, zeta, eta, pf, f0eval, RPQ2, nc1, nc2
  real*8 :: tintegral(1:3,1:3,1:3), expeval, AB(1:3), CD(1:3), PQ(1:3), ubraket1(1:3,1:3),ubraket2(1:3,1:3)
  real*8 :: ubraketket(1:3,1:3,1:3), ubraket(1:3,1:3), uketket(1:3,1:3), ubraint, uketint
  real*8 :: u, uij, uii, uijk, uiij, uiii, uijkl, uiiij, uiipart, ubra(1:3),uket1(1:3),uket2(1:3)
  real*8 :: normtwoelecmat2(0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  character(len=30), intent(in)           :: moleculename

#include "DEBUG.h"

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
call readQchemtwoeecompare(normtwoelecmat2,nobasisfun,basis,moleculename)
print *, "in ggggsppp"
compare = 0;
else
compare = 0;
end if
#endif

  pf = (dacos(-1d0))**(5d0/2d0);

  do i=1,noggsp
    alpha = ggsp(6,i); beta = ggsp(7,i); zeta = ggsp(8,i);
    AB(1:3) = ggsp(11:13,i)
    bf1 = int(ggsp(1,i)); bf2 = int(ggsp(2,i));    nc1 = ggsp(9,i);
    do j=1,noggpp
      gamma = ggpp(6,j); delta = ggpp(7,j); eta = ggpp(8,j);
      bf3 = int(ggpp(1,j)); bf4 = int(ggpp(2,j));      nc2 = ggpp(9,j);
      CD(1:3) = ggpp(11:13,j);      PQ = ggsp(3:5,i)-ggpp(3:5,j);
      RPQ2 = PQ(1)**2+PQ(2)**2+PQ(3)**2;    
tintegral = 0d0;

 if (RPQ2 .lt. 1d-10 ) then
 
do i1=1,3; do i2=1,3; do i3=1,3;
 tintegral(i1,i2,i3) = tintegral(i1,i2,i3) - 2*AB(i1)*CD(i2)*CD(i3)*pf*alpha*gamma*delta/(zeta**2*eta**3*sqrt(zeta+eta));
 if (i2 == i3) then;
   tintegral(i1,i2,i3) = tintegral(i1,i2,i3) + AB(i1)*pf*alpha*(2*zeta+3*eta)/(3*zeta**2*eta**2*(zeta+eta)**1.5d0); end if
 if (i1 == i3) then;   tintegral(i1,i2,i3) = tintegral(i1,i2,i3) - CD(i2)*pf*delta/(3*zeta*eta**2*(zeta+eta)**1.5d0); end if
 if (i1 == i2) then;   tintegral(i1,i2,i3) = tintegral(i1,i2,i3) + CD(i3)*pf*gamma/(3*zeta*eta**2*(zeta+eta)**1.5d0); end if
end do; end do; end do;
!*********************************
! RPQ2 .ne. 0 SECTION
!*********************************
   else
      call calcf0(f0eval,zeta*eta/(zeta+eta)*RPQ2);      expeval = exp(-RPQ2*zeta*eta/(zeta+eta));

  u = 2*pf/(zeta*eta*sqrt(zeta+eta))*f0eval;

  ubraint = 2*pf/(RPQ2*eta*zeta**2*sqrt(zeta+eta))*(expeval-f0eval);
  uketint = 2*pf/(RPQ2*zeta*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);
  ubra = ubraint*PQ*beta;   uket1 = uketint*PQ*gamma;   uket2 = uketint*PQ*delta; 

  uij = 2*pf*gamma*delta/(RPQ2**2*zeta*eta**3*sqrt(zeta+eta))*&
            (-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta)+ 3*f0eval);
  uiipart = -2*pf*gamma*delta/(RPQ2*zeta*eta**3*sqrt(zeta+eta))*(-expeval+f0eval);
  do i1=1,3; do i2=1,3;
   uketket(i1,i2) = uij*PQ(i1)*PQ(i2)
   if (i1 == i2) then;      uketket(i1,i2) = uketket(i1,i2) + uiipart;   end if
  end do; end do

  uij = -2*pf/(RPQ2**2*zeta**2*eta**2*sqrt(zeta+eta))*(-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta)+ 3*f0eval);
  uiipart = 2*pf/(RPQ2*zeta**2*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);

  do i1=1,3; do i2=1,3;
   ubraket1(i1,i2) = beta*gamma*uij*PQ(i1)*PQ(i2);    ubraket2(i1,i2) = beta*delta*uij*PQ(i1)*PQ(i2)
   if (i1 == i2) then;
      ubraket1(i1,i2) = ubraket1(i1,i2) + uiipart*beta*gamma
      ubraket2(i1,i2) = ubraket2(i1,i2) + uiipart*beta*delta
   end if
  end do; end do

 uijk = 2*pf*gamma*delta/(RPQ2**3*zeta**2*eta**3*sqrt(zeta+eta))*(-15*f0eval + (expeval/(zeta+eta)**2)*&
         (15*eta**2 + 10*eta*(3 + eta*RPQ2)*zeta + (15 + 2*eta*RPQ2*(5 + 2*eta*RPQ2))*zeta**2))
  uiij = 2*pf*gamma*delta/(RPQ2**2*zeta**2*eta**3*sqrt(zeta+eta))*&
           (3*f0eval-expeval/(zeta+eta)*(3*(eta+zeta)+2*zeta*eta*RPQ2));

 do i1=1,3; do i2=1,3; do i3=1,3;
  ubraketket(i1,i2,i3) = uijk*PQ(i1)*PQ(i2)*PQ(i3)*beta; 
  if (i1 == i2) then;     ubraketket(i1,i2,i3) = ubraketket(i1,i2,i3) + uiij*PQ(i3)*beta; end if
  if (i1 == i3) then;     ubraketket(i1,i2,i3) = ubraketket(i1,i2,i3) + uiij*PQ(i2)*beta; end if
  if (i2 == i3) then;     ubraketket(i1,i2,i3) = ubraketket(i1,i2,i3) + uiij*PQ(i1)*beta; end if
 end do ; end do; end do


 do i1=1,3; do i2=1,3; do i3=1,3;
   tintegral(i1,i2,i3) = tintegral(i1,i2,i3) -AB(i1)*CD(i2)*CD(i3)*u*alpha*gamma*delta/(zeta*eta**2)

   if (i2 == i3) then;
      tintegral(i1,i2,i3) = tintegral(i1,i2,i3) + AB(i1)*u*alpha/(2*zeta*eta)+ubra(i1)/(4*beta*eta);   end if

  tintegral(i1,i2,i3) = tintegral(i1,i2,i3)+ 1d0/(8*beta*gamma*delta)*ubraketket(i1,i2,i3) &
                         + alpha*AB(i1)*(-CD(i2)*uket2(i3)+CD(i3)*uket1(i2))/(2*zeta*eta) & 
                                        - ubra(i1)*gamma*delta*CD(i2)*CD(i3)/(2*beta*eta**2) & 
                                        + AB(i1)*uketket(i2,i3)*alpha/(4*gamma*delta*zeta) & 
                                      + 1d0/(4*beta*eta)*(-CD(i2)*ubraket2(i1,i3) +CD(i3)*ubraket1(i1,i2))   
                                      
  end do; end do; end do;
 end if

 tintegral = tintegral*nc1*nc2
 printopt =0; compare =0;
 do i1=1,3; do i2=1,3; do i3=1,3; 

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then;
  call makenewtwo(newtwo,nobasisfun,bf1,bf2+i1-1,bf3+i2-1,bf4+i3-1, tintegral(i1,i2,i3),2)
end if
#endif

 call digestintsnewcompare(tintegral(i1,i2,i3),bf1,bf2+i1-1,bf3+i2-1,bf4+i3-1,& 
         Galpha,Gbeta,Ptot,Palpha,Pbeta, nobasisfun,printopt,2, normtwoelecmat2,out,compare)
 end do; end do; end do;

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


end subroutine calcggggsppp



!**********************
