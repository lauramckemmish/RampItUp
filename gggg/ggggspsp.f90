subroutine calcggggspsp(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggsp,noggsp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
 implicit none
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun), ggsp(1:20,1:noggsp)
  real*8, intent(in)      :: basis(1:nobasisfun,1:30)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: noggsp, nobasisfun, i, j, bf1, bf2, bf3, bf4, printopt, braketno, i1, i2
  integer                 :: out,compare, lengthtwoelec
  real*8                  :: alpha, beta, gamma, delta, zeta, eta, pf, f0eval, PQ(1:3), expeval
  real*8                  :: RPQ2, nc1, nc2, tempint(1:3,1:3), AB(1:3), CD(1:3)
  real*8                  :: u, uketint,ubraint,uket2(1:3),ubra2(1:3),ubra2ket2(1:3,1:3), uij, uiipart
  real*8 :: normtwoelecmat2(0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  character(len=30), intent(in)           :: moleculename

#include "DEBUG.h"

#ifdef PSEUDOCHECK
printopt = 0;
if (abs(Palpha(1,1)) .lt. 1d-10) then
call readQchemtwoeecompare(normtwoelecmat2,nobasisfun,basis,moleculename)
print *, "in ggggspsp"
printopt = 0;
compare = 0;
else
compare = 0;
end if
#endif

  pf = (dacos(-1d0))**(5d0/2d0);

  do i=1,noggsp
    alpha = ggsp(6,i); beta = ggsp(7,i); zeta = ggsp(8,i);
    AB = ggsp(11:13,i); 
    bf1 = int(ggsp(1,i)); bf2 = int(ggsp(2,i));    nc1 = ggsp(9,i);
    do j=i,noggsp
      gamma = ggsp(6,j); delta = ggsp(7,j); eta = ggsp(8,j);
      CD = ggsp(11:13,j);      PQ = ggsp(3:5,i) - ggsp(3:5,j);
      bf3 = int(ggsp(1,j)); bf4 = int(ggsp(2,j));      nc2 = ggsp(9,j);
      RPQ2 = (PQ(1))**2+(PQ(2))**2+(PQ(3))**2;

      call calcf0(f0eval,zeta*eta/(zeta+eta)*RPQ2);      expeval = exp(-RPQ2*zeta*eta/(zeta+eta));

 tempint = 0d0;

if (RPQ2 .lt. 1d-10) then; 
 do i1=1,3; do i2=1,3;
   tempint(i1,i2) = 2*AB(i1)*CD(i2)*pf*alpha*gamma/(zeta**2*eta**2*sqrt(zeta+eta));
 end do; end do; 
do i1=1,3; tempint(i1,i1) = tempint(i1,i1) + pf/(3*zeta*eta*(zeta+eta)**(3d0/2d0));end do

else

  u = 2*pf/(zeta*eta*sqrt(zeta+eta))*f0eval;

  ubraint = 2*pf/(RPQ2*eta*zeta**2*sqrt(zeta+eta))*(expeval-f0eval);
  uketint = 2*pf/(RPQ2*zeta*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);
  ubra2 = ubraint*PQ*beta;   uket2 = uketint*PQ*delta; 

  uij = -2*pf/(RPQ2**2*zeta**2*eta**2*sqrt(zeta+eta))*&
            (-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta)+ 3*f0eval);
  uiipart = 2*pf/(RPQ2*zeta**2*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);
  do i1=1,3; do i2=1,3;
   ubra2ket2(i1,i2) = beta*delta*uij*PQ(i1)*PQ(i2)
   if (i1 == i2) then;      ubra2ket2(i1,i2) = ubra2ket2(i1,i2) + uiipart*beta*delta;   end if
  end do; end do


 do i1=1,3; do i2=1,3;
   tempint(i1,i2) =  AB(i1)*CD(i2)*alpha*gamma*u/(zeta*eta) + AB(i1)*alpha*uket2(i2)/(2*delta*zeta) & 
          + CD(i2)*gamma*ubra2(i1)/(2*beta*eta)+ ubra2ket2(i1,i2)/(4*beta*delta); 
 end do; end do;
end if

 if ((bf1 == bf3 .and. bf2 == bf4)) then;
            braketno = 1;            if (i .ne. j) then; nc2 = nc2*2d0; end if;
        else;             braketno = 2;        end if;

 tempint = tempint*nc1*nc2
printopt = 0;

do i1=1,3; do i2=1,3;
 call digestintsnewcompare(tempint(i1,i2),bf1,bf2+i1-1,bf3,bf4+i2-1,Galpha,Gbeta,Ptot,&
    Palpha,Pbeta, nobasisfun,printopt,braketno, normtwoelecmat2,out,compare)

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
  call makenewtwo(newtwo,nobasisfun,bf1,bf2+i1-1,bf3,bf4+i2-1,tempint(i1,i2),braketno)
end if
#endif

end do; end do;

   end do
 end do


#ifdef INTPRINT
if (abs(Palpha(1,1)) .lt. 1d-10) then
!print *, "TESTING DIFFERENCES"
do bf1=1,nobasisfun; do bf2=bf1,nobasisfun; do bf3=bf1,nobasisfun; do bf4=bf2,nobasisfun;
   if (abs(newtwo(bf1,bf2,bf3,bf4,1)) .gt. 1d-10) then;
 if ( abs(newtwo(bf1,bf2,bf3,bf4,1)- normtwoelecmat2(bf1,bf2,bf3,bf4,1)) .gt. 1d-9) then
 !    print *, bf1,bf2,bf3,bf4, newtwo(bf1,bf2,bf3,bf4,1), normtwoelecmat2(bf1,bf2,bf3,bf4,1)
 end if
   end if 
end do; end do; end do; end do;
!print *, "END OF TESTING DIFFERENCES"
end if
#endif

end subroutine calcggggspsp
