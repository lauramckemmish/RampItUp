subroutine kinetic(normkemat,algchange,noRR, noRGSs, noRGSp, RR, RGSsr, RGSpr, nobasisfun,basis,keprint,moleculename, &
    kmax,pts,  ggss, ggsp, ggpp, noggss, noggsp, noggpp)
  implicit none
  integer, intent(in)                 :: nobasisfun,kmax,  noRR, noRGSs, noRGSp, algchange(1:30), pts, keprint
  real*8, intent(in)                  :: RR(1:30,1:noRR), basis(1:nobasisfun,1:30)
  real*8, intent(in)                  :: RGSsr(1:30,1:noRGSs), RGSpr(1:30,1:noRGSp)
  real*8, intent(out)                 :: normkemat(0:nobasisfun,0:nobasisfun,1:2)
  character(len=30)                   :: moleculename
  real*8                              :: ker2SS(1:noRR,1:4), kemat(1:7*nobasisfun+nobasisfun**2+noggss+3*noggsp+9*noggpp,1:4)
  real*8                              :: kergSsg(1:noRGSs,1:4), newpnormalisation
  real*8                              :: kergSpg(1:3*(noRGSp-algchange(2)),1:4), genSpr(1:30,1:noRGSp-algchange(2))
  integer        :: st, fin, no, i, j, totkeread
  integer, intent(in) :: noggss, noggsp, noggpp
  real*8, intent(in)  :: ggss(1:30,1:noggss), ggsp(1:30,1:noggsp), ggpp(1:30,1:noggpp)
  real*8, allocatable, dimension(:,:) :: keggss, keggsp, keggpp


#include "DEBUG.h"
  newpnormalisation = 1d0/(2d0*sqrt(dacos(-1d0)));
  
  kemat = 0d0;
  !*********************************
!!!!!!!!!!!!!!! Master program for calculating kinetic energy integrals 
!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!
  !********* Calculate (S|T|S) ************************
  call r2calcSSke(ker2SS,RR, noRR,nobasisfun)
  ker2SS(1:noRR,1) = ker2SS(1:noRR,1)/(4d0*dacos(-1d0))
  kemat(1:noRR,1:4) = ker2SS(1:noRR,1:4)
   
 !********* Calculate (S|T|s) via quadrature algorithm ************************  
  call genrgSske(kergSsg, noRGSs, RGSsr,pts)
  st = noRR+1;  fin = noRR+noRGSs;
  kergSsg(1: noRGSs,1) = kergSsg(1: noRGSs,1)/(2d0*sqrt(dacos(-1d0)));
  kemat(st:fin,1:4) = kergSsg(1: noRGSs,1:4)
  
  !********* Calculate (S|T|p) via general algorithm (non-concentric cases only) ************************
  no = noRGSp-algchange(2);  genSpr = RGSpr(1:30,algchange(2)+1:noRGSp)
  call genrgSpke(kergSpg,no,genSPr,kmax)
  st = noRR+noRGSs+1;  fin = noRR+noRGSs+3*(noRGSp-algchange(2));
  kergSpg(1:3*no,1) = kergSpg(1:3*no,1)*newpnormalisation;
  kemat(st:fin,1:4) = kergSpg(1:3*no,1:4)
  
 allocate(keggss(1:noggss,1:4), keggsp(1:3*noggsp,1:4), keggpp(1:9*noggpp,1:4))

  call calcggsske(keggss,noggss,ggss)
  st = fin; fin = st + noggss
  kemat(st+1:fin,1:4) = keggss;

  call calcggspke(keggsp,noggsp,ggsp)
  st = fin; fin = st + 3*noggsp
  kemat(st+1:fin,1:4) = keggsp;

  call calcggppke(keggpp,noggpp,ggpp)
  st = fin; fin = st + 9*noggpp
  kemat(st+1:fin,1:4) = keggpp;

deallocate(keggss,keggsp,keggpp)

  call kineticdigestion(normkemat,kemat,nobasisfun,fin,noggss,noggsp,noggpp)

#ifdef WRITEONE
  call  writekinetic(normkemat,nobasisfun,basis,moleculename)
#endif
  
#ifdef ONE_ELEC_PRT
  ! print *, "These are KE integrals"
  !   do i=1,fin
  !     print "(a1,i3,a1,i3,a3,f20.15)", "(", int(kemat(i,2)), ",", int(kemat(i,3)), ") = ", &
  !                kemat(i,1), kemat(i,4)
  !   end do
     
  !  print *, "These are digested KE integrals"
  !   do i=1,nobasisfun
  !      do j=1,nobasisfun
        !   if (normkemat(i,j,2) == 1) then
  !            print "(a1,i3,a1,i3,a3,f20.15)", "(", i, ",", j, ") = ", normkemat(i,j,1)
        !   end if
  !      end do
  !   end do
     
     print *, "This is the kinetic energy matrix" 
     do j=1,nobasisfun,4
        do i=1,nobasisfun
           print "(4f20.15)", normkemat(i,j:min(nobasisfun,j+3),1)
        end do
        print *,""
     end do
#endif


  return
  
end subroutine kinetic


!**********************************************

subroutine kineticdigestion(normkemat,kemat, nobasisfun,nonzero,noggss,noggsp,noggpp)
  implicit none
  real*8, intent(in)    :: kemat(1:7*nobasisfun+nobasisfun**2+noggss+3*noggsp+9*noggpp,1:4)
  real*8, intent(out)   :: normkemat(0:nobasisfun,0:nobasisfun,1:2)
  integer, intent(in)   :: nobasisfun, nonzero,noggss,noggsp,noggpp
  integer               ::  i
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program digests the KE integrals, adding the various components
  ! and normalisation factors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  normkemat = 0d0;
  
  do i=1,nonzero
     normkemat(int(kemat(i,2)),int(kemat(i,3)),1) = normkemat(int(kemat(i,2)),int(kemat(i,3)),1)+kemat(i,1)*kemat(i,4)    
     normkemat(int(kemat(i,3)),int(kemat(i,2)),1) = normkemat(int(kemat(i,2)),int(kemat(i,3)),1)
     normkemat(int(kemat(i,2)),int(kemat(i,3)),2) = 0
     normkemat(int(kemat(i,3)),int(kemat(i,2)),2) = 1  
  end do
  
  return
end subroutine kineticdigestion
