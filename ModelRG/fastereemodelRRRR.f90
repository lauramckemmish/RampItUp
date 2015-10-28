subroutine fastereemodelRRRR(Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,noatoms,at, &
     rrrrpreftable, maxL, maxrampdegree,newsortRg,lenmodelnewsort,maxrgpa,intcnt,Gtot, newtwo,lengthtwoelec)
  implicit none
  real*8, intent(in)                  :: rrrrpreftable(0:maxrampdegree,0:maxrampdegree,0:maxL),Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)                  :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)               :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer, intent(in)                 :: nobasisfun,noatoms, maxL, maxrampdegree,at
  integer                             :: i, j, k, m, a, n1, n2, bf1, bf2, bf3, bf4, L, st,  K2L(1:(1+maxL)**2), n, s, angL, angm
  integer                             :: braketno, pastbf1, pastbf2,lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, w, v, k2, angL2
  integer                             :: warraylast(0:(1+maxL)**2,1:maxrgpa), intcnt, z, printopt, lengthtwoelec
  real*8                              :: newsortRg(1:4,lenmodelnewsort(3,1)+1:lenmodelnewsort(4,maxrgpa))
  real*8                              :: Gtot(1:nobasisfun,1:nobasisfun), nc1,nc2, eeint
  real*8                              :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This subroutine calculates all (RI|RI) on AAAA for all angular momentum. It uses precomputed values in rrrrpreftable.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "DEBUG.h"
#ifdef DEBUG_PRINT
print *, "Start of concentric (RR|RR)"
printopt = 1;
#endif

do i=1,maxrgpa
     warraylast(0,i) = lenmodelnewsort(3,i);
     warraylast(1:(1+maxL)**2,i) = lenmodelnewsort(3,i)
     do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
        k  = newsortRg(3,w)
        warraylast(k,i) = w
        warraylast(k+1:(1+maxL)**2,i) = w
     end do
end do

   do i=1, maxrgpa
     bf1 = lenmodelnewsort(1,i)
     bf2 = lenmodelnewsort(2,i)

     do j=i,maxrgpa
        bf3 = lenmodelnewsort(1,j)
        bf4 = lenmodelnewsort(2,j)
        eeint = 0d0;
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if;

intcnt = intcnt + 1;
        do k=1,(1+maxL)**2
           do w = warraylast(k-1,i)+1,warraylast(k,i)
              nc1 = newsortRg(1,w)
              n1  = newsortRg(2,w)
              angL = newsortRg(4,w)
              do v= warraylast(k-1,j)+1,warraylast(k,j)
                 nc2 = newsortRg(1,v)
                 n2 = newsortRg(2,v)                                  
                 eeint = eeint + nc1*nc2*rrrrpreftable(n1,n2,angL)
              end do
           end do
        end do

#ifdef DEBUG_PRINT
 call printint(eeint,bf1,bf2,bf3,bf4,printopt)
#endif
 
#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
if (i .ne. j .and. braketno == 1) then;
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,2*eeint,braketno)
else 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,eeint,braketno)
end if
end if
#endif

     Gtot(bf1,bf2) = Gtot(bf1,bf2) + Ptot(bf3,bf4)* eeint*braketno
     Gtot(bf3,bf4) = Gtot(bf3,bf4) + Ptot(bf1,bf2)*eeint*braketno
   
     Galpha(bf1,bf4) = Galpha(bf1,bf4) - 0.5d0*braketno*Palpha(bf2,bf3)* eeint 
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - 0.5d0*braketno*Pbeta(bf2,bf3)* eeint
     Galpha(bf1,bf3) = Galpha(bf1,bf3) - 0.5d0*braketno*Palpha(bf2,bf4)*eeint
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - 0.5d0*braketno*Pbeta(bf2,bf4)*eeint
     Galpha(bf2,bf4) = Galpha(bf2,bf4) - 0.5d0*braketno*Palpha(bf1,bf3)*eeint
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - 0.5d0*braketno* Pbeta(bf1,bf3)*eeint
     Galpha(bf2,bf3) = Galpha(bf2,bf3) - 0.5d0*braketno*Palpha(bf1,bf4)*eeint
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - 0.5d0*braketno*Pbeta(bf1,bf4)*eeint

       end do
  end do

 


  return
end subroutine fastereemodelRRRR

