subroutine calcggggssss(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggss,noggss,nobasisfun,basis,moleculename, &
               newtwo ,lengthtwoelec)
 implicit none
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun),ggss(1:20,1:noggss),basis(1:nobasisfun,1:30)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: noggss, nobasisfun, i, j, bf1, bf2, bf3, bf4, printopt, m, n, l, s, braketno
  integer :: out,compare, lengthtwoelec
  real*8   :: zeta, eta, pf, f0eval, tempint, P(1:3), Q(1:3), RPQ2, nc1, nc2
  real*8 :: normtwoelecmat2(0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  character(len=30), intent(in)           :: moleculename

#include "DEBUG.h"

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
call readQchemtwoeecompare(normtwoelecmat2,nobasisfun,basis,moleculename)
compare = 0;   print *, "in ggggssss"
else;   compare = 0;   end if
#endif
printopt = 0;

  pf = 2*(dacos(-1d0))**(5d0/2d0);

  do i=1,noggss
    zeta = ggss(8,i);    P(1:3) = ggss(4:6,i); 
    bf1 = int(ggss(1,i)); bf2 = int(ggss(2,i));    nc1 = ggss(9,i);
    do j=i,noggss
      eta = ggss(8,j);
      Q(1:3) = ggss(4:6,j);       RPQ2 = (P(1)-Q(1))**2+(P(2)-Q(2))**2+(P(3)-Q(3))**2;
      bf3 = int(ggss(1,j)); bf4 = int(ggss(2,j));      nc2 = ggss(9,j);

      call calcf0(f0eval,zeta*eta/(zeta+eta)*RPQ2);
      tempint = pf/(zeta*eta*sqrt(zeta+eta))*f0eval*nc1*nc2

        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; 
           braketno = 1;             if (i .ne. j) then; tempint = tempint*2d0; end if;
        else; braketno = 2;        end if;

    call digestintsnewcompare(tempint,bf1,bf2,bf3,bf4,Galpha,Gbeta,Ptot,Palpha,Pbeta,&
       nobasisfun,printopt,braketno, normtwoelecmat2,out,compare)

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,tempint,braketno)
end if
#endif
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


end subroutine calcggggssss





!******
subroutine makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,intt,braketno)
 implicit none
 integer :: bf1,bf2,bf3,bf4,brano,ketno,braketno, nobasisfun
 real*8 :: intt, newtwo(0:nobasisfun,0: nobasisfun,0: nobasisfun,0: nobasisfun,1:1)

     newtwo(bf1,bf2,bf3,bf4,1) = newtwo(bf1,bf2,bf3,bf4,1) + intt*0.5d0*braketno
     newtwo(bf1,bf2,bf4,bf3,1) = newtwo(bf1,bf2,bf4,bf3,1) + intt*0.5d0*braketno
     newtwo(bf2,bf1,bf3,bf4,1) = newtwo(bf2,bf1,bf3,bf4,1) + intt*0.5d0*braketno
     newtwo(bf2,bf1,bf4,bf3,1) = newtwo(bf2,bf1,bf4,bf3,1) + intt*0.5d0*braketno
     newtwo(bf3,bf4,bf1,bf2,1) = newtwo(bf3,bf4,bf1,bf2,1) + intt*0.5d0*braketno
     newtwo(bf3,bf4,bf2,bf1,1) = newtwo(bf3,bf4,bf2,bf1,1) + intt*0.5d0*braketno   
     newtwo(bf4,bf3,bf1,bf2,1) = newtwo(bf4,bf3,bf1,bf2,1) + intt*0.5d0*braketno
     newtwo(bf4,bf3,bf2,bf1,1) = newtwo(bf4,bf3,bf2,bf1,1) + intt*0.5d0*braketno
     
end subroutine makenewtwo
