subroutine na(combnormnamat,calcmp,algchange,noRR, noRGSs, noRGSp, RR, RGSs, RGSp, nobasisfun,basis,noatoms,atoms,&
     invintatomdist, Tintmat,flatdigmp,L,nLpure,naprint,moleculename,kmax,pts,  &
     ggss, ggsp, ggpp, noggss, noggsp, noggpp)
  implicit none
  integer, intent(in)                     :: nobasisfun, noatoms,kmax, L, nLpure
  integer, intent(in)                     :: naprint, calcmp(1:nobasisfun,1:nobasisfun,1:noatoms)
  integer, intent(in)                     :: noRR, noRGSs, noRGSp, algchange(1:30),pts
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30),invintatomdist(1:noatoms,1:noatoms)
  real*8, intent(in)                      :: Tintmat(1:nLpure,1:nLpure,0:noatoms,0:noatoms), atoms(1:5,1:noatoms)
  real*8, intent(in)                      :: flatdigmp(1:nobasisfun*(nobasisfun+1)/2,1:nLpure+2)
  real*8, intent(in)                      :: RR(1:30,1:noRR), RGSs(1:30,1:noRGSs), RGSp(1:30,1:noRGSp)
  real*8, intent(out)                     :: combnormnamat(0:nobasisfun,0:nobasisfun,1:2)
  integer                                 :: i, j, atid, totnaread, noSP, st, fin, SSSP(1:noRR,1:30), bf1, bf2
  real*8                                  :: nargSsa(1:algchange(1),1:4,1:noatoms)
  real*8                                  :: nargSsg(1:noRGSs-algchange(1),1:4,1:noatoms)
  real*8                                  :: nargSpg(1:3*(noRGSp-algchange(2)),1:4,1:noatoms)
  real*8                                  :: asymSsr(1:30,1:algchange(1)), genSsr(1:30,1:noRGSs-algchange(1))
  real*8                                  :: genSpr(1:30,1:(noRGSp-algchange(2))), naRR(1:noRR,1:4,1:noatoms)
  real*8                                  :: normnamat(0:nobasisfun,0:nobasisfun,1:2,1:noatoms), newpnormalisation
  character(len=30)                       :: moleculename
  real*8,dimension(:,:,:),allocatable     :: namp, namat
  integer, intent(in)                     :: noggss, noggsp, noggpp
  real*8, intent(in)                      :: ggss(1:30,1:noggss), ggsp(1:30,1:noggsp), ggpp(1:30,1:noggpp)
  real*8, allocatable, dimension(:,:,:)   :: naggss, naggsp, naggpp

  newpnormalisation = 1d0/(2d0*sqrt(dacos(-1d0)));
  
  !*******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Master program for calculating nuclear attraction integrals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
#include "DEBUG.h"
allocate(namp(1:nobasisfun,1:nobasisfun,1:noatoms))

  call calcnamp(namp,flatdigmp,nobasisfun,noatoms,Tintmat,nLpure,L)

  do atid = 1,noatoms
     normnamat(1:nobasisfun,1:nobasisfun,1,atid) = -atoms(2,atid)*&
          namp(1:nobasisfun,1:nobasisfun,atid)*calcmp(1:nobasisfun,1:nobasisfun,atid)   
  end do
  normnamat(1:nobasisfun,1:nobasisfun,2,1:noatoms) = calcmp(1:nobasisfun,1:nobasisfun,1:noatoms)
  
deallocate(namp)

 allocate(namat(1:7*nobasisfun+4*nobasisfun**2+noggss+3*noggsp+9*noggpp,1:4,1:noatoms))
namat = 0d0;
  !********* Calculate (S|V|S) ************************
    
  ! ******  Then correct for those that are.
  call r2calcconcSSna(naRR,RR, noRR, noatoms)
  naRR(1:noRR,1,1:noatoms) = naRR(1:noRR,1,1:noatoms)/(4d0*dacos(-1d0));  namat(1:noRR,1:4,1:noatoms) = naRR
  
  !********* Calculate (S|V|s) ************************
  ! *** General case
  noSP = noRGSs;   genSsr = RGSs(1:30,1:noRGSs)
  
  ! Correcting for those that are concentric S & V
  call genrgSsna(nargSsg,noSP,genSsr,noatoms,pts)
  st = noRR+1;  fin = noRR+noSP;
  nargSsg(1:noSP,1,1:noatoms) = nargSsg(1:noSP,1,1:noatoms)/(2d0*sqrt(dacos(-1d0)))
do i=1,noRgSs
 if (int(RgSs(14,i)) == 111) then; 
   print *, "IMPORTANT - na integrals", int(RgSs(1:2,i)), naRgSsg(i,1,1:noatoms)
 end if
end do;
  namat(st:fin,1:4,1:noatoms) = nargSsg(1:noSP,1:4,1:noatoms)
  !*********************
  
  !********* Calculate (S|V|p) where S, V are concentric (S, V, p concentric =0 and are not calculated)************************
  noSP = noRGSp-algchange(2);  genSpr = RGSp(1:30,algchange(2)+1:noRGSp)
  call genrgSpna(nargSpg,noSP,genSpr,noatoms,kmax)
  st = noRR+noRGSs+1;  fin = noRR+noRGSs+3*(noRGSp-algchange(2));
  nargSpg(1:3*noSP,1,1:noatoms) = nargSpg(1:3*noSP,1,1:noatoms)*newpnormalisation;
  namat(st:fin,1:4,1:noatoms) = nargSpg(1:3*noSP,1:4,1:noatoms)

 allocate(naggss(1:noggss,1:4,1:noatoms), naggsp(1:3*noggsp,1:4,1:noatoms), naggpp(1:9*noggpp,1:4,1:noatoms))
  call calcggssna(naggss,noggss,ggss,noatoms,atoms)
  st = fin; fin = st + noggss;    namat(st+1:fin,1:4,1:noatoms) = naggss;
  call calcggspna(naggsp,noggsp,ggsp,noatoms,atoms)
  st = fin; fin = st + 3*noggsp;  namat(st+1:fin,1:4,1:noatoms) = naggsp;
  call calcggppna(naggpp,noggpp,ggpp,noatoms,atoms)
  st = fin; fin = st + 9*noggpp;  namat(st+1:fin,1:4,1:noatoms) = naggpp;

 deallocate(naggss, naggsp,naggpp)

! The digestion is inline now!!
  do atid = 1,noatoms
     do i=1,fin
        bf1 = min(int(namat(i,2,atid)), int(namat(i,3,atid))); 
        bf2 = max(int(namat(i,2,atid)), int(namat(i,3,atid))); 
        normnamat(bf1,bf2,1,atid) = normnamat(bf1,bf2,1,atid)-atoms(2,atid)*namat(i,1,atid)*namat(i,4,atid)    
        normnamat(bf2,bf1,1,atid) = normnamat(bf1,bf2,1,atid)
     end do
  end do
deallocate(namat)

combnormnamat=0d0;
    do atid = 1,noatoms
     combnormnamat(1:nobasisfun,1:nobasisfun,1) = combnormnamat(1:nobasisfun,1:nobasisfun,1) &
          + normnamat(1:nobasisfun,1:nobasisfun,1,atid)
  end do
 

#ifdef WRITEONE
  call  writena(combnormnamat,nobasisfun,basis,moleculename)
#endif 

#ifdef ONE_ELEC_PRT
     
!     print *, ""; print *, "These are digested NA integrals for each atom"
!     do atid = 1,noatoms;   do i=1,nobasisfun;   do j=i,nobasisfun
      !  if (int(normnamat(i,j,2,atid)) == 1) then
 !          print "(a1,i3,a2,i2,a1,i2,a3,f20.15)", "(", i, "|C ", atid,"|", j,") = ", normnamat(i,j,1,atid)
      !  end if;  
!     end do;    end do;   end do
     
!     print *, ""; print *, "These are digested total NA integrals "
!     do i=1,nobasisfun;    do j=i,nobasisfun
  !      if (combnormnamat(i,j,2) == 1) then
 !          print "(a1,i3,a1,i3,a3,f20.15)", "(", i, ",", j,") = ", combnormnamat(i,j,1)
  !      end if;   
!     end do;    end do
     
     print *, ""; print *, "This is the nuclear attraction matrix" 
     do j=1,nobasisfun,7;   do i=1,nobasisfun
        print "(7f20.15)", combnormnamat(i,j:min(nobasisfun,j+6),1)
     end do;   print *,"";    end do;      print *, ""
     
#endif
  
  return
end subroutine na


