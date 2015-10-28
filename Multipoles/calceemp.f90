subroutine calceemp(Galpha,Gbeta,Palpha,Pbeta,Ptot,flatdigmp,nobasisfun,noatoms,Tintmat,nLpure,intcnt, Gtot, & 
           newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)                     :: nobasisfun,noatoms, nLpure
  real*8, intent(in)                      :: flatdigmp(1:nobasisfun*(nobasisfun+1)/2,1:nLpure+2)
  real*8, intent(in)                      :: Tintmat(1:nLpure,1:nLpure,0:noatoms,0:noatoms)
  real*8, intent(in)                      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)                      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)                   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8                                  :: Gtot(1:nobasisfun,1:nobasisfun), eempint ,Ptotbf3bf4, x
  integer                                 :: l2, l1, angm2, angm, intcnt, a,b, lengthtwoelec, braketno
  integer                                 :: multipolecenter, i, j,k, m, at1, at2, bf1, bf2, bf3, bf4, j2, k2, n, l, s
  integer, allocatable, dimension(:)      :: index
  integer, allocatable, dimension(:,:,:)  :: listtodo
  real*8, allocatable, dimension(:,:)     :: onlynecmp1, onlynecmp2, eempintmat, midway
  real*8    :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)

 allocate(listtodo(1:3,1:9*nobasisfun,1:noatoms))
 allocate(index(1:noatoms))
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This program calculates electron-electron repulsion integrals via multipole moment interactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "DEBUG.h"

#ifdef DEBUG_PRINT
print *, "start of non-concentric (RI|RI)"
#endif
  index = 0
  j2 = 0;
  do bf3 = 1,nobasisfun; do bf4= bf3,nobasisfun
     j2 = j2+1;
     at1 = int(flatdigmp(j2,nLpure+1));
   if (at1 .gt. 0) then;
      index(at1) = index(at1)+1;
      listtodo(1,index(at1),at1) = bf3;
      listtodo(2,index(at1),at1) = bf4;
      listtodo(3,index(at1),at1) = j2;
   end if
 end do; end do;


do at1 = 1,noatoms

if (index(at1) == 0) then;  cycle;   end if
 allocate( onlynecmp1(1:index(at1),1:nLpure))
 onlynecmp1(1:index(at1),1:nLpure) = flatdigmp(listtodo(3,1:index(at1),at1),1:nLpure)
 allocate( midway(1:nLpure,1:index(at1)))
  do at2 = 1,noatoms

if (index(at2) == 0) then;  cycle;  end if
 allocate( onlynecmp2(1:index(at2),1:nLpure))
 onlynecmp2(1:index(at2),1:nLpure) = flatdigmp(listtodo(3,1:index(at2),at2),1:nLpure)

 call dgemm('N','T',nLpure,index(at1),nLpure,1d0,Tintmat(1:nLpure,1:nLpure,at2,at1),&
     nLpure,onlynecmp1(1:index(at1),1:nLpure),index(at1),0d0,midway,nLpure)

 allocate(eempintmat(1:index(at2),1:index(at1)))

 call dgemm('N','N',index(at2),index(at1),nLpure,1d0,onlynecmp2,index(at2),midway,&
         nLpure,0d0,eempintmat,index(at2))

 do a=1,index(at1)
   bf3 = listtodo(1,a,at1)
   bf4 = listtodo(2,a,at1)
  Ptotbf3bf4 = 2d0*Ptot(bf3,bf4)
 x = 0d0; 
    do b=1,index(at2)
       bf1 = listtodo(1,b,at2)
       bf2 = listtodo(2,b,at2)
        if (bf1 > bf3 .or. (bf1 == bf3 .and. bf2 .gt. bf4)) then;

intcnt = intcnt + 1;
eempint = eempintmat(b,a)

#ifdef DEBUG_PRINT
 call printint(eempint,bf1,bf2,bf3,bf4,1)
#endif

x = x+ 2d0*Ptot(bf1,bf2)*eempint

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;        end if;
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,eempint,braketno)
end if
#endif

     Gtot(bf1,bf2) = Gtot(bf1,bf2) + Ptotbf3bf4* eempint
     
     Galpha(bf1,bf4) = Galpha(bf1,bf4) - Palpha(bf2,bf3)* eempint 
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - Pbeta(bf2,bf3)* eempint
     Galpha(bf1,bf3) = Galpha(bf1,bf3) -Palpha(bf2,bf4)*eempint
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - Pbeta(bf2,bf4)*eempint
     Galpha(bf2,bf4) = Galpha(bf2,bf4) -Palpha(bf1,bf3)*eempint
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - Pbeta(bf1,bf3)*eempint
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -Palpha(bf1,bf4)*eempint
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - Pbeta(bf1,bf4)*eempint

   end if; 
end do; 
Gtot(bf3,bf4) = Gtot(bf3,bf4) +x;  

end do;
deallocate(onlynecmp2,eempintmat)
end do; 

deallocate(midway,onlynecmp1)
end do;
 
deallocate(index,listtodo)

  return
end subroutine calceemp
