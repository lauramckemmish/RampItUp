subroutine twoelecRgoneRperatom(Galpha,Gbeta,Palpha,Pbeta,Ptot, noggss, noggsp, noggpp, ggss,ggsp, ggpp, &
   nobasisfun,noatoms, atoms,  maxrampdegree,  Tcutoff, pts,a,&
   newsortRg,lenmodelnewsort,maxrgpa, lenmodelcontnewsort, contnewsortRg,lencont, &
   Omrnmatrix,RIssmain, RIspmain, RIppmain,integrand, integrandsp,   integrandpp, intcnt, Gtot,absc,weights, &
   r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow, totpartpp, & 
   newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     ::  maxrampdegree, a, nobasisfun, noatoms, pts
  integer, intent(in)     :: noggss, noggsp, noggpp
  real*8, intent(in)      ::  Tcutoff, atoms(1:5,1:noatoms)
  real*8, intent(in)      :: ggss(1:20,1:noggss), ggsp(1:20,1:noggsp), ggpp(1:20,1:noggpp)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: nonzero, lenmodelcontnewsort(1:4,1:maxrgpa), lencont, lengthtwoelec
  integer                 :: lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, intcnt(1:20) 
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun) 
  real*8                  :: newsortRg(1:4,lenmodelnewsort(3,1)+1:lenmodelnewsort(4,maxrgpa))
  real*8                  :: contnewsortRg(1:4,1: lencont), absc(1:2*pts+1), weights(1:2*pts+1)
  real*8                  :: Omrnmatrix(1:2*pts+1,0:maxrampdegree), RIssmain(0:maxrampdegree,0:4)
  real*8                  :: integrand(1:2*pts+1,0:4), integrandsp(1:2*pts+1,0:4), integrandpp(1:2*pts+1,0:4)
  real*8                  :: RIspmain(0:maxrampdegree,0:4), RIppmain(0:maxrampdegree,0:4)
  real*8  :: r1(1:2*pts+1), r2(1:2*pts+1), r3(1:2*pts+1), r4(1:2*pts+1),r6(1:2*pts+1), r8(1:2*pts+1), r10(1:2*pts+1),&
      r12(1:2*pts+1),expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1),&
      evenpow(1:2*pts+1), oddpow(1:2*pts+1)
  real*8   ::  totpartpp(1:9,1:25,0:maxrampdegree), ang(1:25), dangdx(1:25), dangdy(1:25), dangdz(1:25)
  real*8   :: newtwo(0:lengthtwoelec,0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 
  real     :: t1, t2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates two-electron integrals
  ! 
  ! The output for each type is in terms of (unnorm integral, bf1, bf2, normalisation factor)
  !
  ! The normtwoelecmat is in the form (bf1, bf2, bf3, bf4, integral value (normalised)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !****************************************************************************************
  !************** Calculating mixed Ramp-Gaussian integrals **************
  !****************************************************************************************
  
  !****************************************************************************************
  !************** Ramp-only shell pairs and Gaussian-only Shell pairs **************
  !****************************************************************************************

!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  !!                                RIgg INTEGRALS
!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!

!!!********!!!!!!!!*****!!!!!!!!!!!! LR (T > 20) ********!!!!!!!!!!!!************!
  call calcLRgenRIssdigestoneRperatom(noggss, ggss, noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa,lenmodelcontnewsort,&
      contnewsortRg, Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, intcnt(1),Gtot,ang,totpartpp, newtwo,lengthtwoelec)


!!!********!!!!!!!!*****!!!!!!!!!!!! SR (T < 20 & zeta <= betachangealg) ********!!!!!!!!!!!!************!
   call calcSRgenRIssdigestoneRperatom(noggss, ggss, noatoms, atoms, maxrampdegree, Tcutoff, pts,a, newsortRg,&
     lenmodelnewsort,maxrgpa, Omrnmatrix, RIssmain,  Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,integrand, intcnt(2), &
     Gtot,absc,weights, r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow,ang, newtwo,lengthtwoelec)

!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  !!                                RIsp INTEGRALS
!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  
!!!********!!!!!!!!*****!!!!!!!!!!!! LR (T > 20) ********!!!!!!!!!!!!************!
call calcLRgenRIspdigestoneRperatom(noggsp, ggsp, noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, lenmodelcontnewsort,&
      contnewsortRg, Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, intcnt(3),Gtot, totpartpp, newtwo,lengthtwoelec)

!!!********!!!!!!!!*****!!!!!!!!!!!! SR (T < 20 & zeta <= betachangealg) ********!!!!!!!!!!!!************!
 call calcSRgenRIspdigestoneRperatom(noggsp, ggsp, noatoms, atoms, maxrampdegree, Tcutoff,pts,a, newsortRg,& 
     lenmodelnewsort,maxrgpa, Omrnmatrix,Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, integrand, integrandsp, intcnt(4), &
     Gtot,absc,weights, RIssmain, RIspmain,r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow, &
     totpartpp, newtwo,lengthtwoelec)

!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  !!                                RIpp INTEGRALS
!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  
!!!********!!!!!!!!*****!!!!!!!!!!!! LR ********!!!!!!!!!!!!************!
  call calcLRgenRIppdigestoneRperatom(noggpp, ggpp, noatoms, atoms, maxrampdegree, Tcutoff, a,maxrgpa, lenmodelcontnewsort,&
     contnewsortRg,  Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,intcnt(5),Gtot, totpartpp, newtwo,lengthtwoelec)
  
!!!********!!!!!!!!*****!!!!!!!!!!!! SR (T < 20 & zeta <= betachangealg) ********!!!!!!!!!!!!************!
  
 call calcSRgenRIppdigestoneRperatom(noggpp, ggpp,  noatoms, atoms,  maxrampdegree, Tcutoff, pts,a, newsortRg,&
      lenmodelnewsort,maxrgpa, Omrnmatrix,  Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, integrand, integrandsp, &
      integrandpp, intcnt(6), Gtot,absc,weights, RIssmain, RIspmain, RIppmain, &
     r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow,totpartpp, newtwo,lengthtwoelec)


 return
end subroutine twoelecRgoneRperatom











subroutine twoelecRg(Galpha,Gbeta,Palpha,Pbeta,Ptot, noggss, noggsp, noggpp, ggss,ggsp, ggpp, nobasisfun,noatoms,&
   atoms,  maxrampdegree,  Tcutoff, pts,a,newsortRg,lenmodelnewsort,maxrgpa, lenmodelcontnewsort, contnewsortRg,&
   lencont, Omrnmatrix,RIssmain, integrand, integrandsp,   integrandpp, intcnt, Gtot, newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     ::  maxrampdegree, a, nobasisfun, noatoms, pts, noggss, noggsp, noggpp
  real*8, intent(in)      ::  Tcutoff, atoms(1:5,1:noatoms), Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: ggss(1:20,1:noggss), ggsp(1:20,1:noggsp), ggpp(1:20,1:noggpp)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: nonzero, lenmodelcontnewsort(1:4,1:maxrgpa), lencont, lengthtwoelec
  integer                 :: lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, intcnt(1:20)  
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun), contnewsortRg(1:4,1: lencont)
  real*8                  :: newsortRg(1:4,lenmodelnewsort(3,1)+1:lenmodelnewsort(4,maxrgpa))
  real*8                  :: Omrnmatrix(1:2*pts+1,0:maxrampdegree), RIssmain(0:maxrampdegree,0:4)
  real*8                  :: integrand(1:2*pts+1,0:4), integrandsp(1:2*pts+1,0:4), integrandpp(1:2*pts+1,0:4)
  real*8                  :: ang(1:25), dangdx(1:25), dangdy(1:25), dangdz(1:25), d2angdth2(1:25)
  real*8                  :: d2angdx2(1:25), d2angdxy(1:25), d2angdxz(1:25), d2angdph2(1:25),d2angdthph(1:25)
  real*8                  :: d2angdy2(1:25), d2angdyz(1:25), d2angdz2(1:25),dangdth(1:25), dangdph(1:25)
  real*8                  :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates two-electron integrals
  ! 
  ! The output for each type is in terms of (unnorm integral, bf1, bf2, normalisation factor)
  !
  ! The normtwoelecmat is in the form (bf1, bf2, bf3, bf4, integral value (normalised)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !****************************************************************************************
  !************** Calculating mixed Ramp-Gaussian integrals **************
  !****************************************************************************************
  !****************************************************************************************
  !************** Ramp-only shell pairs and Gaussian-only Shell pairs **************
  !****************************************************************************************

!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  !!                                RIgg INTEGRALS
!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!

!!!********!!!!!!!!*****!!!!!!!!!!!! LR (T > 20) ********!!!!!!!!!!!!************!
  call calcLRgenRIssdigest(noggss, ggss, noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, lenmodelcontnewsort, &
    contnewsortRg,lencont, Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,ang, intcnt(1),Gtot, newtwo,lengthtwoelec)

!!!********!!!!!!!!*****!!!!!!!!!!!! SR (T < 20 & zeta <= betachangealg) ********!!!!!!!!!!!!************!
   call calcSRgenRIssdigest(noggss, ggss, noatoms, atoms, maxrampdegree, Tcutoff, pts,a, newsortRg,&
     lenmodelnewsort,maxrgpa, Omrnmatrix, RIssmain, &
     Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,integrand,ang, intcnt(2),Gtot, newtwo,lengthtwoelec)

 !!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  !!                                RIsp INTEGRALS
!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  
!!!********!!!!!!!!*****!!!!!!!!!!!! LR (T > 20) ********!!!!!!!!!!!!************!
call calcLRgenRIspdigest(noggsp, ggsp, noatoms, atoms, maxrampdegree, Tcutoff,&
     a,maxrgpa, lenmodelcontnewsort, contnewsortRg, Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, &
     ang,dangdx, dangdy, dangdz,dangdth,dangdph, intcnt(3),Gtot, newtwo,lengthtwoelec)

!!!********!!!!!!!!*****!!!!!!!!!!!! SR (T < 20 & zeta <= betachangealg) ********!!!!!!!!!!!!************!
 call calcSRgenRIspdigest(noggsp, ggsp, noatoms, atoms, maxrampdegree, Tcutoff,pts,a, newsortRg,&
     lenmodelnewsort,maxrgpa, Omrnmatrix,  Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, integrand, integrandsp,&
     ang,dangdx, dangdy, dangdz,dangdth,dangdph, intcnt(4),Gtot, newtwo,lengthtwoelec)

!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  !!                                RIpp INTEGRALS
!!!!!!!!!!!!************!!!!!!!!!!!!!****************!!!!!!!!!!!!!!***********!!!!!!!!!!!!!!
  
!!!********!!!!!!!!*****!!!!!!!!!!!! LR ********!!!!!!!!!!!!************!
  call calcLRgenRIppdigest(noggpp, ggpp, noatoms, atoms, maxrampdegree, Tcutoff,&
     a,maxrgpa, lenmodelcontnewsort, contnewsortRg,  Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, &
     ang,dangdx, dangdy, dangdz, d2angdx2, d2angdxy,d2angdxz,d2angdy2,d2angdyz,d2angdz2, & 
     dangdth,dangdph,d2angdth2,d2angdph2,d2angdthph,intcnt(5),Gtot, newtwo,lengthtwoelec)
  
!!!********!!!!!!!!*****!!!!!!!!!!!! SR (T < 20 & zeta <= betachangealg) ********!!!!!!!!!!!!************!
  
 call calcSRgenRIppdigest(noggpp, ggpp,  noatoms, atoms, maxrampdegree, Tcutoff, pts,a, newsortRg,&
      lenmodelnewsort,maxrgpa, Omrnmatrix, Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, integrand, integrandsp, &
      integrandpp, ang,dangdx, dangdy, dangdz, d2angdx2, d2angdxy,d2angdxz,d2angdy2,d2angdyz,d2angdz2, & 
      dangdth,dangdph,d2angdth2,d2angdph2,d2angdthph,intcnt(6),Gtot, newtwo,lengthtwoelec)

 return
end subroutine twoelecRg



!!!!!!!!!!!!!!
!!!!!!!!!!!!!!
