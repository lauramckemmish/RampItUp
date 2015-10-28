
subroutine calcLRgenRIssdigestoneRperatom(noggss, ggss,  &
     noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, lenmodelcontnewsort, contnewsortRg, &
     Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,intcnt,Gtot,ang,intt, & 
           newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     :: noatoms, maxrampdegree, noggss,a,nobasisfun
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: atoms(1:5,1:noatoms), Tcutoff, ggss(1:20,1:noggss)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  integer                 :: i,j, ntot, k, atid, l, m, n, swit, warnswit(1:10), bf1, bf2, bf3, bf4
  integer                 :: lenmodelcontnewsort(1:3,1:maxrgpa), maxrgpa, v, st, intcnt
  real*8                  :: nc2, zeta, sqrtpir4, Ax, Ay, Az, QA, QA2, QAx, QAy, QAz,rQA
  real*8                  :: W(0:4), pf1(0:4), ang(1:25), T, totrf
  real*8                  :: contnewsortRg(1: 25*maxrgpa ), x, y, Ptotbf3bf4,  tot, stor(1:4)
  real*8                  :: temp111, intt(1:maxrgpa) 
 real*8, dimension(:), allocatable :: acum2, stor2
 integer, dimension(:), allocatable :: rowbf2
  integer :: lengthtwoelec,braketno,cnt,cnt2
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 
 
allocate(acum2(1:maxrgpa), stor2(1:maxrgpa))
allocate(rowbf2(1:maxrgpa))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LRgenRIss calculates integrals of the form (SS|ss) for cases where only the long-range component is important
  ! i.e. the (RI| and |ss) shell pairs do not overlap significantly (as quantified Qy the T parameter)
  !
  ! All integrals are correct. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "DEBUG.h"
  warnswit = 0;
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!***********
  !!      PRECOMPUTED QUANTITES
  !!   These could be precomputed and given as a table, but easier to compute them here 
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!**********
  do l=0,4;     pf1(l) = 16*(dacos(-1d0)**2)/(2d0*l+1);    end do
  sqrtpir4 = sqrt(dacos(-1d0))*0.25d0;

                 bf1 = lenmodelcontnewsort(1,1)
  
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!***********
  !!      SHORT-RANGE --> LONG-RANGE REPRESENTATION
  !!   This could be done at a shell pair level, but it is done here for convenience - move out later
  !!   if this becomes a time bottleneck
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!**********

 
!!!!!!!***********!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!!!!********* 
  !!   MAIN LOOP OF SHELL QUARTETS
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!**********
rowbf2 = lenmodelcontnewsort(2,1:maxrgpa);
stor2(1:maxrgpa) = Ptot(rowbf2,bf1)
acum2 = 0d0;

        Ax     = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)
  do j=1,noggss
        QAx    = ggss(4,j) - Ax;    QAy     = ggss(5,j) - Ay;          QAz     = ggss(6,j)-Az;
        QA2    = (QAx)**2+(QAy)**2+(QAz)**2 ; 
         zeta     = ggss(8,j)
        T      = zeta*QA2;
        
        !************ Decides whether to do it via long-range methodology; 
        if (T .ge. Tcutoff) then;
           x = 0d0;
     bf3      = ggss(1,j);        bf4      = ggss(2,j);
     Ptotbf3bf4 = 2d0*Ptot(bf3,bf4)
     QA = sqrt(QA2);   rQA = 1/QA;

           call calcfullang(ang,QAx,QAy,QAz);

           W(0) = sqrtpir4*rQA*ggss(14,j)*ggss(9,j);

           do l=1,4;     W(l) = W(l-1)*rQA;       end do

           !********* Loop over Ramps centered on atomic center
ang(1) = ang(1)*pf1(0)*W(0)
ang(2:4) = ang(2:4)*pf1(1)*W(1)
ang(5:9) = ang(5:9)*pf1(2)*W(2)
ang(10:16) = ang(10:16)*pf1(3)*W(3)
ang(17:25) = ang(17:25)*pf1(4)*W(4)

#ifdef INTCNTPRT
intcnt = intcnt + maxrgpa;
#endif

call dgemv('T',25,maxrgpa,1d0,contnewsortRg(lenmodelcontnewsort(3,1)+1: lenmodelcontnewsort(3,maxrgpa)+25), &
       25,ang,1,0d0, intt,1)

stor(1) = Palpha(bf3,bf1)
stor(2) = Pbeta(bf3,bf1)
stor(3) = Palpha(bf4,bf1)
stor(4) = Pbeta(bf4,bf1)

        do i=1,maxrgpa
                 bf2 = rowbf2(i)
    
     Galpha(bf2,bf4) = Galpha(bf2,bf4) -stor(1)*intt(i)
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  -stor(2)*intt(i)
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -stor(3)*intt(i)
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  -stor(4)*intt(i)


#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if;  
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,intt(i),braketno)
end if
#endif

        end do
    acum2 = acum2 + Ptotbf3bf4*intt

     Galpha(bf1,bf4) = Galpha(bf1,bf4) -sum(Palpha(rowbf2,bf3)*intt)
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - sum(Pbeta(rowbf2,bf3)*intt)
     Galpha(bf1,bf3) = Galpha(bf1,bf3) -sum(Palpha(rowbf2,bf4)*intt)
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - sum(Pbeta(rowbf2,bf4)*intt)
     Gtot(bf3,bf4) = Gtot(bf3,bf4) + sum(2d0*stor2*intt);
        end if
  end do
  
do i=1,maxrgpa
  bf2 = rowbf2(i)
 Gtot(bf2,bf1) = Gtot(bf2,bf1) + acum2(i)
end do
#ifdef INTCNTPRT
print *, "Integral cnt in LRgenRIss = ", intcnt
#endif
deallocate(rowbf2,stor2,acum2)

return
end subroutine calcLRgenRIssdigestoneRperatom










!*****************************



subroutine calcLRgenRIssdigest(noggss, ggss,  &
     noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, lenmodelcontnewsort, contnewsortRg, lencont, &
     Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, ang,intcnt,Gtot)
  implicit none
  integer, intent(in)     :: noatoms, maxrampdegree, noggss,a,nobasisfun
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: atoms(1:5,1:noatoms), Tcutoff, ggss(1:20,1:noggss)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  integer                 :: i,j, ntot, k, atid, l, m, n, swit, warnswit(1:10), bf1, bf2, bf3, bf4
  integer                 :: lenmodelcontnewsort(1:3,1:maxrgpa), maxrgpa, lencont, v, st, intcnt
  real*8                  :: nc2, zeta, sqrtzeta, sqrtpir4, Ax, Ay, Az, Qx, Qy, Qz, QA, QA2, QAx, QAy, QAz,rQA
  real*8                  :: W(0:7), pf1(0:7), ang(1:25), T, totrf
  real*8                  :: contnewsortRg(1: lencont ), intt(1:maxrgpa), x, y, Ptotbf3bf4, inttest(1:maxrgpa), tot
 real*8 :: temp(1:25*maxrgpa), tempang(1:25)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! LRgenRIss calculates integrals of the form (SS|ss) for cases where only the long-range component is important
  ! i.e. the (RI| and |ss) shell pairs do not overlap significantly (as quantified Qy the T parameter)
  !
  ! All integrals are correct. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "DEBUG.h"
  warnswit = 0;
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!***********
  !!      PRECOMPUTED QUANTITES
  !!   These could be precomputed and given as a table, but easier to compute them here 
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!**********
  do l=0,4;     pf1(l) = 16*(dacos(-1d0)**2)/(2d0*l+1);    end do
  sqrtpir4 = sqrt(dacos(-1d0))*0.25d0;
  
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!***********
  !!      SHORT-RANGE --> LONG-RANGE REPRESENTATION
  !!   This could be done at a shell pair level, but it is done here for convenience - move out later
  !!   if this becomes a time bottleneck
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!**********

 
!!!!!!!***********!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!!!!********* 
  !!   MAIN LOOP OF SHELL QUARTETS
!!!!!!!!**************!!!!!!!!!!!!!!!!!!!!!!!!*********************!!!!!!!!!!!!!!**********
        Ax     = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)
  do j=1,noggss
     bf3      = ggss(1,j);        bf4      = ggss(2,j);
     nc2      = ggss(9,j);
     Qx       = ggss(4,j);        Qy       = ggss(5,j);       Qz       = ggss(6,j)
     sqrtzeta = ggss(12,j);       zeta     = ggss(8,j)
     Ptotbf3bf4 = 2d0*Ptot(bf3,bf4)
        QAx    = Qx - Ax;    QAy     = Qy - Ay;          QAz     = Qz-Az;
        QA2    = (QAx)**2+(QAy)**2+(QAz)**2 ; QA = sqrt(QA2);   rQA = 1/QA;
        T      = zeta*QA2;
        
        !************ Decides whether to do it via long-range methodology; 
        if (T .ge. Tcutoff) then;
intt = 0d0; ! DONT DELETE THIS LINE!!!
           x = 0d0;
           call calcfullang(ang,QAx,QAy,QAz);

           W(0) = sqrtpir4/(QA*zeta**(3d0/2d0));

           do l=1,4;     W(l) = W(l-1)*rQA;       end do

           !********* Loop over Ramps centered on atomic center
k = 0; do l=0,4; do m=-l,l; k=k+1; ang(k) = ang(k)*W(l)*pf1(l)*nc2      ; end do; end do;

intcnt = intcnt + maxrgpa;

call dgemv('T',25,maxrgpa,1d0,contnewsortRg(lenmodelcontnewsort(3,1)+1: lenmodelcontnewsort(3,maxrgpa)+25), &
       25,ang,1,10d0, intt,1)

        do i=1,maxrgpa
                 bf1 = lenmodelcontnewsort(1,i)
                 bf2 = lenmodelcontnewsort(2,i)
     y =Ptotbf3bf4* intt(i) 
     Gtot(bf1,bf2) = Gtot(bf1,bf2) + y
    
     Galpha(bf1,bf4) = Galpha(bf1,bf4) -Palpha(bf2,bf3)*intt(i)
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - Pbeta(bf2,bf3)*intt(i)
     Galpha(bf1,bf3) = Galpha(bf1,bf3) -Palpha(bf2,bf4)*intt(i)
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - Pbeta(bf2,bf4)*intt(i)
     Galpha(bf2,bf4) = Galpha(bf2,bf4) -Palpha(bf1,bf3)*intt(i)
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - Pbeta(bf1,bf3)*intt(i)
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -Palpha(bf1,bf4)*intt(i)
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - Pbeta(bf1,bf4)*intt(i)

     x = x+2d0*Ptot(bf1,bf2)*intt(i);
        end do

     Gtot(bf3,bf4) = Gtot(bf3,bf4) + x
        end if
  end do
  
#ifdef INTCNTPRT
print *, "Integral cnt in LRgenRIss = ", intcnt
#endif

return
end subroutine calcLRgenRIssdigest

