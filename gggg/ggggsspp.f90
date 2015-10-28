subroutine calcggggsspp(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggss,ggpp,noggss,noggpp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
 implicit none
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: ggss(1:20,1:noggss),ggpp(1:20,1:noggpp), basis(1:nobasisfun,1:30)
  integer                 :: noggss, noggpp,nobasisfun, i, j, bf1, bf2, bf3, bf4, printopt, braketno, i1, i2
  integer  :: out,compare, lengthtwoelec
  real*8   :: alpha, beta, gamma, delta, zeta, eta, pf, f0eval,RPQ2, nc1, nc2, CD(1:3), PQ(1:3)
  real*8   :: tempint(1:3,1:3), uiipart,uij, u, uketint, uket1(1:3), uket2(1:3), uketket(1:3,1:3),expeval
  real*8   :: normtwoelecmat2(0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  real*8   :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  character(len=30), intent(in)           :: moleculename

#include "DEBUG.h"

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
print *, "in ggggsspp"
call readQchemtwoeecompare(normtwoelecmat2,nobasisfun,basis,moleculename);  compare = 0;
else;  compare = 0;   end if
#endif


  pf = (dacos(-1d0))**(5d0/2d0);

  do i=1,noggss
    alpha = ggss(14,i); beta = ggss(15,i); zeta = ggss(8,i);
    bf1 = int(ggss(1,i)); bf2 = int(ggss(2,i));    nc1 = ggss(9,i);
    do j=1,noggpp
      gamma = ggpp(6,j); delta = ggpp(7,j); eta = ggpp(8,j);
      CD(1:3) = ggpp(11:13,j);
      bf3 = int(ggpp(1,j)); bf4 = int(ggpp(2,j));      nc2 = ggpp(9,j);
      PQ = ggss(4:6,i)-ggpp(3:5,j);       RPQ2 = (PQ(1))**2+(PQ(2))**2+(PQ(3))**2;

      call calcf0(f0eval,zeta*eta/(zeta+eta)*RPQ2);      expeval = exp(-RPQ2*zeta*eta/(zeta+eta));

 tempint = 0d0;

if (RPQ2 .lt. 1d-10) then;
do i1=1,3; tempint(i1,i1) = pf*(2*zeta+3*eta)/(3*zeta*eta**2*(zeta+eta)**(3d0/2d0)); end do;
do i1=1,3; do i2=1,3;
  tempint(i1,i2)=tempint(i1,i2)-2*CD(i1)*CD(i2)*pf*gamma*delta/(zeta*eta**3*sqrt(zeta+eta));
end do; end do;

else
  u = 2*pf/(zeta*eta*sqrt(zeta+eta))*f0eval;

  uketint = 2*pf/(RPQ2*zeta*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);
  uket1 = uketint*PQ*gamma;   uket2 = uketint*PQ*delta; 

  uij = 2*pf*gamma*delta/(RPQ2**2*zeta*eta**3*sqrt(zeta+eta))*&
            (-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta)+ 3*f0eval);
  uiipart = -2*pf*gamma*delta/(RPQ2*zeta*eta**3*sqrt(zeta+eta))*(-expeval+f0eval);
  do i1=1,3; do i2=1,3;
   uketket(i1,i2) = uij*PQ(i1)*PQ(i2)
   if (i1 == i2) then;      uketket(i1,i2) = uketket(i1,i2) + uiipart;   end if
  end do; end do

do i1=1,3; tempint(i1,i1) = u/(2*eta);end do;

do i1=1,3; do i2=1,3;
  tempint(i1,i2)=tempint(i1,i2)-CD(i1)*CD(i2)*gamma*delta*u/(eta**2) + &
                   uketket(i1,i2)/(4*gamma*delta)  - CD(i1)*uket2(i2)/(2*eta) + CD(i2)*uket1(i1)/(2*eta)
end do; end do;

end if

        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then;
            braketno = 1;   if (i .ne. j) then; nc2 = nc2*2d0; end if;
        else; braketno = 2;        end if;

 tempint = tempint*nc1*nc2
printopt = 0;
do i1=1,3; do i2=1,3;
   call digestintsnewcompare(tempint(i1,i2),bf1,bf2,bf3+i1-1,bf4+i2-1,&
                Galpha,Gbeta,Ptot,Palpha,Pbeta, nobasisfun,printopt,braketno, normtwoelecmat2,out,compare)

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3+i1-1,bf4+i2-1,tempint(i1,i2),braketno)
end if
#endif

end do; end do;

   end do
 end do




end subroutine calcggggsspp
