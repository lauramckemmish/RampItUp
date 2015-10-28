subroutine twoelec(Galpha,Gbeta,Palpha,Pbeta,Ptot ,flatdigmp,Tintmat,&
     totnomodel, nLpure,algchange, noggss, noggsp, noggpp, ggss,ggsp, ggpp, nobasisfun, basis,noatoms, atoms,&
     moleculename, RRRRpreftable, maxL, maxrampdegree, Tcutoff, pts, compareprint,&
     newsortRg,lenmodelnewsort,maxrgpa,norampsperatom, totnoramps ,rampinfo)

  implicit none
  character(len=30), intent(in)       :: moleculename
  integer, intent(in)     :: norampsperatom(1:noatoms), totnoramps, rampinfo(1:2,1:totnoramps)
  integer, intent(in)     :: maxL, maxrampdegree, nLpure, totnomodel
  integer, intent(in)     :: nobasisfun, noatoms,  algchange(1:20), pts
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(out)     :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: RRRRpreftable(0:maxrampdegree,0:maxrampdegree,0:maxL), Tcutoff
  real*8, intent(in)      :: basis(1:nobasisfun,1:30), atoms(1:5,1:noatoms)
  real*8, intent(in)      :: flatdigmp(1:nobasisfun*(nobasisfun+1)/2,1:nLpure+2)
  integer                 :: lenmodelcontnewsort(1:3,1:maxrgpa), stindex, finindex, w, oldk, a,i, j, k, bfR
  integer                 :: noggss, noggsp, noggpp, compareprint, l, n, bf1, bf2, timecnt, arraysize 
  integer                 :: at2, finss, stss, finsp, stsp, finpp, stpp, shortnoggss, shortnoggsp, shortnoggpp
  integer                 :: atid, ntot, lenmodelnewsort(1:4,1:maxrgpa,1:noatoms), maxrgpa, tempmaxrgpa, intcnt(1:20)
  integer                  :: toteeread, lengthtwoelec, bf3, bf4
  real*8                  :: ggss(1:20,1:noggss), ggsp(1:20,1:noggsp), ggpp(1:20,1:noggpp)
  real*8                  :: Tintmat(1:25,1:25,0:noatoms,0:noatoms), ramppreffactor(0:maxrampdegree,0:4)
  real*8                  :: nc1, rf, newsortRg(1:4,1: totnomodel), contnewsortRg(1:maxrgpa*25)
  real*8                  :: absc(1:2*pts+1), weights(1:2*pts+1), Gtot(1:nobasisfun,1:nobasisfun)
 real*8, dimension(:,:), allocatable :: Omrnmatrix, RIssmain, RIspmain, RIppmain,integrand, integrandsp, integrandpp
 real*8, dimension(:), allocatable :: r1, r2, r3, r4, r6, r8,r10,r12,oddpow,evenpow,erfQpr,erfQmr,expQmr,expQpr
 real*8, dimension(:,:,:), allocatable :: totpartpp
  real*8,dimension(:,:,:,:,:), allocatable :: newtwo
  real*8, dimension(:,:), allocatable :: qchemee
  real                    :: t(1:5), t2(1:5), tarray(1:2,1:5), totaltime(1:5)

#include "DEBUG.h"

  Galpha = 0d0;        Gbeta = 0d0;        Gtot = 0d0;   intcnt = 0;

  shortnoggss = ceiling(dble(noggss)/dble(noatoms))
  shortnoggsp = ceiling(dble(noggsp)/dble(noatoms))
  shortnoggpp = ceiling(dble(noggpp)/dble(noatoms))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates two-electron integrals
  ! 
  ! The output for each type is in terms of (unnorm integral, bf1, bf2, normalisation factor)
  !
  ! The normtwoelecmat is in the form (bf1, bf2, bf3, bf4, integral value (normalised)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !****************************************************************************************
!!!************************************ Reads all-Gaussian two-electron integrals from Qchem**************
  !****************************************************************************************
! Must do this first with the multipole calls of twoelecmatdigestion (which is aimed to increase speed of calculation).

  timecnt = 1;
  call etime(tarray(1:2,timecnt),totaltime(timecnt)); t(timecnt) = tarray(1,timecnt); t2(timecnt) = tarray(2,timecnt)


if (abs(Palpha(1,1)) .lt. 1d-10) then;
 lengthtwoelec = nobasisfun;
else
 lengthtwoelec = 1;
end if
 allocate(newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1))
 newtwo = 0d0;
#ifdef KEEPGGGG
call calcggmain(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggss,ggsp,ggpp,noggss,noggsp,noggpp,nobasisfun,&
  basis,moleculename, lengthtwoelec,newtwo,noatoms)
#endif
  timecnt = timecnt +1;
call etime(tarray(1:2,timecnt),totaltime(timecnt)); t(timecnt) = tarray(1,timecnt); t2(timecnt) = tarray(2,timecnt)
  !****************************************************************************************
  !************** Calculating two-electron integrals via multipole-moments **************
  !***************************************************************************************
  call calceemp(Galpha,Gbeta,Palpha,Pbeta,Ptot,flatdigmp,nobasisfun,noatoms,Tintmat,nLpure,intcnt(7),Gtot, newtwo,lengthtwoelec)

allocate(r1(1:2*pts+1), r2(1:2*pts+1), r3(1:2*pts+1), r4(1:2*pts+1),r6(1:2*pts+1), r8(1:2*pts+1), r10(1:2*pts+1),r12(1:2*pts+1))
allocate(expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1),evenpow(1:2*pts+1), oddpow(1:2*pts+1))

allocate(Omrnmatrix(1:2*pts+1,0:maxrampdegree), RIssmain(0:maxrampdegree,0:4), &
      RIspmain(0:maxrampdegree,0:4), RIppmain(0:maxrampdegree,0:4))
allocate(integrand(1:2*pts+1,0:4), integrandsp(1:2*pts+1,0:4), integrandpp(1:2*pts+1,0:4))
allocate(totpartpp(1:9,1:25,0:maxrampdegree))
  do n=0,maxrampdegree; do l=0,maxL;
     ramppreffactor(n,l) = exp(dlgama(3d0+2d0*l)+dlgama(1d0+n)-dlgama(4d0+2d0*l+n))
  end do; end do

  call getquad(weights,absc,pts);
  do ntot = 0,maxrampdegree;    call evalOmrn(Omrnmatrix(1:2*pts+1,ntot),absc,pts,ntot);  end do
call evalallrk(r1,r2,r3,r4, r6, r8,absc,pts);  r10 = r8*r2;  r12 = r8*r4; 

  do a = 1,noatoms

if (norampsperatom(a) > 0) then;

tempmaxrgpa = maxrgpa
do i=maxrgpa,1,-1
   if (lenmodelnewsort(1,i,a) == 0) then
      tempmaxrgpa = i
   end if
end do
  call fastereemodelRRRR(Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,noatoms, a, RRRRpreftable,  maxL, &
       maxrampdegree, newsortRg(1:4,lenmodelnewsort(3,1,a)+1:lenmodelnewsort(4,tempmaxrgpa,a)),&
       lenmodelnewsort(1:4,1: tempmaxrgpa,a:a), tempmaxrgpa,intcnt(8),Gtot, newtwo,lengthtwoelec)

if (norampsperatom(a) == 1) then;

 do w = 1,totnoramps
   if (rampinfo(2,w) == a) then; 
       bfR = rampinfo(1,w)
   end if
 end do

   stindex = 0; finindex = 0;
   contnewsortRg(1:25*tempmaxrgpa) = 0d0;
   do i = 1,tempmaxrgpa
     oldk = 0;
     stindex = finindex 
     do w = lenmodelnewsort(3,i,a)+1, lenmodelnewsort(4,i,a)
        ntot    = int(newsortRg(2,w))
        k       = int(newsortRg(3,w))
        l       = int(newsortRg(4,w))
        rf      = ramppreffactor(ntot,l)       
        nc1     = newsortRg(1,w);
        contnewsortRg(stindex + k) = contnewsortRg(stindex+k) + nc1*rf
     end do
     finindex = stindex + 25;

if (int(bfR) ==  lenmodelnewsort(1,i,a)) then;
   bf1 = bfR
   bf2 = lenmodelnewsort(2,i,a);
else if (bfR ==  lenmodelnewsort(2,i,a)) then;
   bf1 = bfR
   bf2 = lenmodelnewsort(1,i,a);
end if

     lenmodelnewsort(1:2,i,a) = (/bf1,bf2/)

     lenmodelcontnewsort(1:3,i) = (/bf1,bf2,stindex/)
  end do

finss = 0;   finsp = 0;   finpp = 0;
do at2=1,noatoms
 stss = finss;          finss = min(stss+shortnoggss,noggss);
 stsp = finsp;          finsp = min(stsp+shortnoggsp,noggsp);
 stpp = finpp;          finpp = min(stpp+shortnoggpp,noggpp);

 call twoelecRgoneRperatom( Galpha,Gbeta,Palpha,Pbeta,Ptot, &
     finss-stss, finsp-stsp, finpp-stpp, ggss(1:20,stss+1:finss),ggsp(1:20,stsp+1:finsp),&
     ggpp(1:20,stpp+1:finpp), nobasisfun,noatoms, atoms, maxrampdegree, Tcutoff, pts,a, &
     newsortRg(1:4,lenmodelnewsort(3,1,a)+1:lenmodelnewsort(4,tempmaxrgpa,a)),&
     lenmodelnewsort(1:4,1: tempmaxrgpa,a:a), tempmaxrgpa, &
     lenmodelcontnewsort(1:3,1:tempmaxrgpa), contnewsortRg(1: finindex), finindex, &
     Omrnmatrix, RIssmain, RIspmain, RIppmain, integrand,integrandsp,integrandpp, intcnt, Gtot,absc,weights, &
     r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow, totpartpp, & 
           newtwo,lengthtwoelec)
end do

else

   stindex = 0; finindex = 0;
   contnewsortRg(1:25*tempmaxrgpa) = 0d0;
   do i = 1,tempmaxrgpa
     oldk = 0;
     stindex = finindex 
     do w = lenmodelnewsort(3,i,a)+1, lenmodelnewsort(4,i,a)
        ntot    = int(newsortRg(2,w))
        k       = int(newsortRg(3,w))
        l       = int(newsortRg(4,w))
        rf      = ramppreffactor(ntot,l)       
        nc1     = newsortRg(1,w);
        contnewsortRg(stindex + k) = contnewsortRg(stindex+k) + nc1*rf
     end do
     finindex = stindex + 25;
     bf1     = lenmodelnewsort(1,i,a);     bf2     = lenmodelnewsort(2,i,a);
     lenmodelcontnewsort(1,i) = bf1
     lenmodelcontnewsort(2,i) = bf2
     lenmodelcontnewsort(3,i) = stindex
  end do

finss = 0;   finsp = 0;   finpp = 0;
do at2=1,noatoms
 stss = finss;          finss = min(stss+shortnoggss,noggss);
 stsp = finsp;          finsp = min(stsp+shortnoggsp,noggsp);
 stpp = finpp;          finpp = min(stpp+shortnoggpp,noggpp);

 call twoelecRg( Galpha,Gbeta,Palpha,Pbeta,Ptot, &
     finss-stss, finsp-stsp, finpp-stpp, ggss(1:20,stss+1:finss),ggsp(1:20,stsp+1:finsp),&
     ggpp(1:20,stpp+1:finpp), nobasisfun,noatoms, atoms, maxrampdegree, Tcutoff, pts,a, &
     newsortRg(1:4,lenmodelnewsort(3,1,a)+1:lenmodelnewsort(4,tempmaxrgpa,a)),&
     lenmodelnewsort(1:4,1: tempmaxrgpa,a:a), tempmaxrgpa, &
     lenmodelcontnewsort(1:3,1:tempmaxrgpa), contnewsortRg(1: finindex), finindex, &
     Omrnmatrix, RIssmain, integrand,integrandsp,integrandpp, intcnt, Gtot, newtwo,lengthtwoelec)
end do

end if

end if


end do
deallocate(RIssmain,Omrnmatrix,integrand,integrandsp, integrandpp,totpartpp)
deallocate(r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow)

  !****************************************************************************************
!!!*********************** Puts all the various kinds of two electron integrals into a two-electron integral matrix**************
  !****************************************************************************************
  timecnt = timecnt+1;
  call etime(tarray(1:2,timecnt),totaltime(timecnt)); t(timecnt) = tarray(1,timecnt); t2(timecnt) = tarray(2,timecnt)

Galpha = Galpha + Gtot;   Gbeta  = Gbeta  + Gtot;

do i=1,nobasisfun;   Galpha(i,i) = 2d0*Galpha(i,i);  Gbeta(i,i) = 2d0*Gbeta(i,i);end do

do i=1,nobasisfun; do j=i+1,nobasisfun
 Galpha(i,j) = Galpha(i,j) + Galpha(j,i);    Gbeta(i,j) = Gbeta(i,j) + Gbeta(j,i);
 Galpha(j,i) = Galpha(i,j);                  Gbeta(j,i) = Gbeta(i,j);
end do; end do; 

  print *, "***!*@&*(@*&!#*&!#(!#(*@*(@*(@(@*@*&!*^!^&@&^@&^!#*("
  print "(a40,f20.12,a5,f20.12,a5)", &
        "Time for two-elec integral computation is", t(timecnt) -t(1), "secs", t2(timecnt) -t2(1), "secs"
  print "(a40,f20.12,a5)", "Time for gggg two-elec integral computation is", t(2)-t(1), "secs"
  print "(a40,f20.12,a5)", "Time for R-containing two-elec integral comp. is", t(timecnt)-t(2), "secs"
  print *, "***!*@&*(@*&!#*&!#(!#(*@*(@*(@(@*@*&!*^!^&@&^@&^!#*("


#ifdef INTPRINT
if (abs(Palpha(1,1)) .lt. 1d-10) then

#ifdef PSEUDOCHECK
allocate(qchemee(1:6,1:nobasisfun**4))
call readQchempseudotwoee(qchemee,nobasisfun,toteeread,moleculename)
print *, "TESTING AGAINST PSEUDOTWOEE"
print *, "ERRORS ARE"
print *, "               PSEUDO      RAMPITUP         DIF" 
do i=1,nobasisfun**4
bf1 = qchemee(2,i)
bf2 = qchemee(3,i)
bf3 = qchemee(4,i)
bf4 = qchemee(5,i)
if (i .gt. 2 .and. bf1 ==qchemee(2,1) .and. bf2 == qchemee(3,1) &
           .and. bf3 == qchemee(4,1) .and. bf4 == qchemee(5,1)) then;
  go to 200
end if
if (abs(qchemee(1,i)-newtwo(bf1,bf2,bf3,bf4,1)) .gt. 1d-7) then
 print *, bf1,bf2,bf3,bf4,qchemee(1,i),newtwo(bf1,bf2,bf3,bf4,1),qchemee(1,i)-newtwo(bf1,bf2,bf3,bf4,1)
end if
end do
200 print *, "i = ", i, "nobasisfun = ", nobasisfun, "exited out"
deallocate(qchemee)
#endif
print *, "finish TESTING"
print *, ""
#ifdef WRITETWOEE
if (lengthtwoelec .gt. 1) then
 call writetwoee(newtwo,nobasisfun,basis,moleculename)
print *, "written twoelec integrals to Integrals.csv"
end if
#endif
print *, "*@(U$(*@*)$(*!@)($*@!)93*@!)(#*@!)(#*@!)(#*!)@(*#)(@!#"


 do i=1,nobasisfun
  do j=i,nobasisfun
    do k=i,nobasisfun;
     do l=k,nobasisfun;
      if (abs(newtwo(i,j,k,l,1)) .gt. 1d-10) then;
!       print *, i, j, k, l, newtwo(i,j,k,l,1)
      end if
     end do
    end do
  end do
 end do
end if
#endif

deallocate(newtwo)

#ifdef INTCNTPRT 
print "(a40,i10)", "Int cnt for LR (R|ss) = ", intcnt(1)
print "(a40,i10)", "Int cnt for SR (R|ss) = ", intcnt(2)
print "(a40,i10)", "Int cnt for LR (R|sp) = ", 3*intcnt(3) 
print "(a40,i10)", "Int cnt for SR (R|sp) = ", 3*intcnt(4) 
print "(a40,i10)", "Int cnt for LR (R|pp) = ", 9*intcnt(5) 
print "(a40,i10)", "Int cnt for SR (R|pp) = ", 9*intcnt(6)
print "(a40,i10)", "Int cnt for LR (R|R) = ", intcnt(7)
print "(a40,i10)", "Int cnt for SR (R|R) = ", intcnt(8)
print "(a40,i10)", "Int cnt for (Ss|S) = ", intcnt(9)
print "(a40,i10)", "Int cnt for (Sp|P) = ", 9*intcnt(10)

print "(a40,i10)", "Int cnt for (Ss|ss) = ", intcnt(11)
print "(a40,i10)", "Int cnt for (Ss|sp) = ", 3*intcnt(12)
print "(a40,i10)", "Int cnt for (Ss|pp) = ", 9*intcnt(13)
print "(a40,i10)", "Int cnt for (Sp|ss) = ", 3*intcnt(14)
print "(a40,i10)", "Int cnt for (Sp|sp) = ", 9*intcnt(15)
print "(a40,i10)", "Int cnt for (Sp|pp) = ", 27*intcnt(16)


print "(a40,i10)", "Int cnt for (Ss|Ss) = ", intcnt(17)
print "(a40,i10)", "Int cnt for (Sp|Sp) = ", 9*intcnt(18)
#endif
  return
end subroutine twoelec























