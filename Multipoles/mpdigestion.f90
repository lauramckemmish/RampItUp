subroutine mpdigestion(flatdigmp,calcmp,mpRRSS, mpRGSs,mpRGSp,nobasisfun, noRR,noSs,noSp,nLpure,noatoms)
  implicit none
  integer, intent(in)                    :: nLpure,noRR, noSs, noSp, nobasisfun, noatoms
  real*8, intent(in)                     :: mpRRSS(1:noRR,1:nLpure+4), mpRGSs(1:noSs,1:nLpure+4), mpRGSp(1:3*noSp,1:nLpure+4)
  real*8, intent(out)                    :: flatdigmp(1:nobasisfun*(1+nobasisfun)/2,1:nLpure+2)
  integer, intent(out)                   :: calcmp(1:nobasisfun,1:nobasisfun,1:noatoms)
  real*8                                 :: digmp(0:nobasisfun,0:nobasisfun,1:nLpure+2)
  integer                                :: i, j, cnt, bf1, bf2, bf3, bf4, at1, at2, k
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program digests the multipole moments of (SS|, (Ss| and (Sp| shell pairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  digmp= 0d0;  flatdigmp =0d0;  calcmp = 0;
  
  do i=1,noRR
     digmp(int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+2)),1:nLpure) = &
          digmp(int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+2)),1:nLpure) + mpRRSS(i,1:nLpure)*mpRRSS(i,nLpure+3)
     digmp(int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+2)),nLpure+1) = mpRRSS(i,nLpure+4);
     digmp(int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+2)),nLpure+2) = 1;
     
     digmp(int(mpRRSS(i,nLpure+2)),int(mpRRSS(i,nLpure+1)),1:nLpure+2) = &
          digmp(int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+2)),1:nLpure+2);
     
     ! This says that all na integrals involving this shell pair are done via mp interactions
     calcmp(int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+2)),1:noatoms) = 1
     calcmp(int(mpRRSS(i,nLpure+2)),int(mpRRSS(i,nLpure+1)),1:noatoms) = 1
     
     ! This says that the na integral where the ramp and nuclei are concentric cannot be done like this;
     calcmp(int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+2)),int(mpRRSS(i,nLpure+4))) = 0
     calcmp(int(mpRRSS(i,nLpure+2)),int(mpRRSS(i,nLpure+1)),int(mpRRSS(i,nLpure+4))) = 0
  end do
  
  do i=1,noSs
     digmp(int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+2)),1:nLpure) = &
          digmp(int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+2)),1:nLpure) + mpRGSs(i,1:nLpure)*mpRGSs(i,nLpure+3)
     digmp(int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+2)),nLpure+1) = mpRGSs(i,nLpure+4);
     digmp(int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+2)),nLpure+2) = 1;
     
     digmp(int(mpRGSs(i,nLpure+2)),int(mpRGSs(i,nLpure+1)),1:nLpure+2) = &
          digmp(int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+2)),1:nLpure+2);
     
     ! This says that all na integrals involving this shell pair are done via mp interactions
     calcmp(int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+2)),1:noatoms) = 1
     calcmp(int(mpRGSs(i,nLpure+2)),int(mpRGSs(i,nLpure+1)),1:noatoms) = 1
     
     ! This says that the na integral where the ramp and nuclei are concentric cannot be done like this;
     calcmp(int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+2)),int(mpRGSs(i,nLpure+4))) = 0
     calcmp(int(mpRGSs(i,nLpure+2)),int(mpRGSs(i,nLpure+1)),int(mpRGSs(i,nLpure+4))) = 0
  end do
  
  do i=1,3*noSp
     digmp(int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+2)),1:nLpure) = &
          digmp(int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+2)),1:nLpure) + mpRGSp(i,1:nLpure)*mpRGSp(i,nLpure+3)
     digmp(int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+2)),nLpure+1) = mpRGSp(i,nLpure+4);
     digmp(int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+2)),nLpure+2) = 1;
     
     digmp(int(mpRGSp(i,nLpure+2)),int(mpRGSp(i,nLpure+1)),1:nLpure+2) = &
          digmp(int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+2)),1:nLpure+2) 
     
     ! This says that all na integrals involving this shell pair are done via mp interactions
     calcmp(int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+2)),1:noatoms) = 1
     calcmp(int(mpRGSp(i,nLpure+2)),int(mpRGSp(i,nLpure+1)),1:noatoms) = 1
     
     ! This says that the na integral where the ramp and nuclei are concentric cannot be done like this;
     calcmp(int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+2)),int(mpRGSp(i,nLpure+4))) = 0
     calcmp(int(mpRGSp(i,nLpure+2)),int(mpRGSp(i,nLpure+1)),int(mpRGSp(i,nLpure+4))) = 0
  end do
  
  cnt = 0;
  do i=1,nobasisfun;    do j=i,nobasisfun
     cnt = cnt + 1;
      flatdigmp(cnt,1:nLpure+2) = digmp(i,j,1:nLpure+2);   
   end do;  end do
   
   return
 end subroutine mpdigestion
