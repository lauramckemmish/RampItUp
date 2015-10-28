subroutine modelRg(totnomodel, &
     noRR,noRgSs,noRgSp, RR,RgSs,RgSp,nobasisfun,noatoms,maxmodellen,modelRgprint, algchange, &
     maxrampdegree,maxL,kmax,lenmdlconc, concmodeloption,rrrrpreftable, Lnoncon,mdlchk,&
     newsortRg, lenmodelnewsort,maxrgpa)
  
  implicit none
  integer, intent(in)  :: maxrampdegree, maxL,kmax, Lnoncon,  algchange(1:30), concmodeloption
  integer, intent(in)  :: modelRgprint,noRgSs, noRgSp, noRR, mdlchk, nobasisfun, noatoms, maxmodellen, lenmdlconc
  real*8, intent(in)   :: rrrrpreftable(0:maxrampdegree,0:maxrampdegree,0:maxL)
  real*8, intent(in)   :: RR(1:30,1:noRR), RgSs(1:30,1:noRgSs), RgSp(1:30,1:noRgSp)
  integer, intent(out) :: totnomodel
  real*8, dimension(:,:,:), allocatable :: shellmodelRgmat
  real*8, dimension(:,:), allocatable   :: shellmodelids
  integer, dimension(:), allocatable    :: lenshellmodel
  real*8, dimension(:,:), allocatable   :: RRshellmodelids
  real*8, dimension(:,:), allocatable   :: bfmodelids, concRgSs, concRgSp, concSsshellmodelids, concSpshellmodelids
  real*8, dimension(:,:), allocatable   :: nonconcRgSs, nonconcRgSp, nonconcSsshellmodelids, nonconcSpshellmodelids
  real*8, dimension(:,:,:), allocatable :: basisfunmodelRgmat, concSsshellmodelRgmat
  real*8, dimension(:,:,:), allocatable :: RRshellmodelRgmat
  real*8, dimension(:,:,:), allocatable :: concSpshellmodelRgmat, nonconcSsshellmodelRgmat, nonconcSpshellmodelRgmat
  real*8,dimension(:,:,:,:),allocatable :: normmodelRgmat
  integer, dimension(:), allocatable    :: lenbfmodel, concSslenshellmodel, nonconcSslenshellmodel 
  integer, dimension(:), allocatable    :: concSplenshellmodel, nonconcSplenshellmodel, RRlenshellmodel
  integer, dimension(:), allocatable    :: RRSPlenshellmodel, RRPPlenshellmodel
  integer              :: i, j, K, Lcnt, Lx, Ly, Lz, bf1, bf2, totno, nohere, st, fin, cnt, a
  integer              :: lennormmodel(1:nobasisfun,1:nobasisfun,1:2), timecnt, maxlenbf, maxrgpa
  real*8               :: t(1:10)
  real*8, intent(out)  :: newsortRg(1:4,1:nobasisfun*noatoms*maxmodellen)
  integer, intent(out) :: lenmodelnewsort(1:4,1:10*nobasisfun,1:noatoms)
  integer, dimension(:), allocatable     :: RgPslenshellmodel, RgPplenshellmodel
  real*8, dimension(:,:), allocatable    :: RgPsshellmodelids, RgPpshellmodelids
  real*8, dimension(:,:,:), allocatable :: RgPsshellmodelRgmat, RgPpshellmodelRgmat



  allocate(shellmodelRgmat(1:noRgSs&
                            +3*(noRgSp-algchange(2)) + algchange(2) +noRR,1:3,1:maxmodellen))
  allocate(shellmodelids(1:noRgSs+3*(noRgSp-algchange(2)) &
                                + algchange(2)+noRR,1:4))
  allocate(lenshellmodel(1:noRgSs+3*(noRgSp-algchange(2))+ algchange(2)+noRR))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program oversees the creation of model ramps for mixed Gaussian-ramp shell pairs 
  ! (also adds ramp-ramp shell pair data to be processed consecutively) 
  ! Ramp models are defined by the basis functions they come, their atomic centre and their angular momentum. 
  ! There is also a coefficient associated with them which takes into account all contraction coefficients
  ! and normalisation constants
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  timecnt = 0;  timecnt = timecnt+1;  call cpu_time(t(timecnt))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!***** Model concentric (Ss| where s is not too large by sum of S-ramps******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(concRgSs(1:30,1:algchange(4)),concSsshellmodelRgmat(1:algchange(4),1:3,1:lenmdlconc))
  allocate(concSsshellmodelids(1:algchange(4),1:4),concSslenshellmodel(1:algchange(4)))

  concRgSs(1:30,1:algchange(4))     = RgSs(1:30,1:algchange(4))
  nohere                                          = algchange(4)
  if (mdlchk == 1) then;  print *, ""; print *, "Modelling Concentric Ss shell pairs"; end if
  call concmdlctrl(concSsshellmodelRGmat,concSslenshellmodel,concSsshellmodelids, nohere,&
       concRgSs, 13, mdlchk)
  st = 0;                                         fin = algchange(4)
  
  shellmodelRgmat(st+1:fin,1:3,1:lenmdlconc)      = concSsshellmodelRgmat
  lenshellmodel(st+1:fin)                          = concSslenshellmodel
  shellmodelids(st+1:fin,1:4)                      = concSsshellmodelids
  
  timecnt = timecnt+1;  call cpu_time(t(timecnt))
  print *, "Time for conc Ss modelling", t(timecnt) - t(timecnt-1)

deallocate(concRgSs , concSsshellmodelids, concSsshellmodelRgmat, concSslenshellmodel)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***** Model non-concentric (Ss| where s is not too large by sum of Ramps of varying angular momentum ******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(nonconcSsshellmodelids(1:noRgSs-algchange(4),1:4))
  allocate(nonconcRgSs(1:30,1:noRgSs-algchange(4)))
  allocate(nonconcSsshellmodelRgmat(1:noRgSs-algchange(4),1:3,1:maxmodellen))
  allocate(nonconcSslenshellmodel(1:noRgSs-algchange(4)))

  nonconcRgSs(1:30,1:noRgSs-algchange(4))     = RgSs(1:30,algchange(4)+1:noRgSs)
  nohere                                       = noRgSs-algchange(4)

  if (nohere .gt. 0) then;
  call Ssmdlctrl(nonconcSsshellmodelRGmat,nonconcSslenshellmodel,nonconcSsshellmodelids, nohere,&
       nonconcRgSs, maxmodellen,kmax, lenmdlconc,Lnoncon,  concmodeloption,rrrrpreftable,maxrampdegree, maxL, mdlchk)

  st = fin; fin = fin+noRgSs-algchange(4)
  shellmodelRgmat(st+1:fin,1:3,1:maxmodellen)      = nonconcSsshellmodelRgmat
  lenshellmodel(st+1:fin)                          = nonconcSslenshellmodel
  shellmodelids(st+1:fin,1:4)                      = nonconcSsshellmodelids
  end if

deallocate(nonconcRgSs,nonconcSsshellmodelRgmat, nonconcSslenshellmodel, nonconcSsshellmodelids)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***** Model concentric (Sp| where s is not too large by sum of P-ramps ******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(concRgSp(1:30,1:algchange(2)),concSpshellmodelRgmat(1:algchange(2),1:3,1:lenmdlconc))
  allocate(concSpshellmodelids(1:algchange(2),1:4),concSplenshellmodel(1:algchange(2)))
  st = fin;                                        fin = st + algchange(2)
  concRgSp(1:30,1:algchange(2))      = RgSp(1:30,1:algchange(2))
  nohere                                           = algchange(2)
  
  call concmdlctrl(concSpshellmodelRGmat,concSplenshellmodel,concSpshellmodelids, nohere, concRgSp, 13, mdlchk)
  
  shellmodelRgmat(st+1:fin,1:3,1:lenmdlconc)      = concSpshellmodelRgmat
  lenshellmodel(st+1:fin)                          = concSplenshellmodel
  shellmodelids(st+1:fin,1:4)                      = concSpshellmodelids

deallocate(concRgSp , concSpshellmodelids , concSpshellmodelRgmat, concSplenshellmodel)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***** Model non-concentric (Sp| where p is not too large by sum of Ramps of varying angular momentum ******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(nonconcRgSp(1:30,1:noRgSp-algchange(2)))
  allocate(nonconcSpshellmodelRgmat(1:3*(noRgSp-algchange(2)),1:3,1:maxmodellen))
  allocate(nonconcSpshellmodelids(1:3*(noRgSp-algchange(2)),1:4))
  allocate(nonconcSplenshellmodel(1:3*(noRgSp-algchange(2))))

  st = fin; fin = fin+3*(noRgSp-algchange(2))
  nonconcRgSp(1:30,1:(fin-st)/3)              = RgSp(1:30,algchange(2)+1:noRgSp)
  nohere                                      = (noRgSp - algchange(2))
  if (mdlchk == 1) then;  print *,""; print *, "Modelling Non-concentric Sp shell pairs";   end if;

  if (nohere .gt. 0) then;
  call Spmdlctrl(nonconcSpshellmodelRGmat,nonconcSplenshellmodel,nonconcSpshellmodelids, nohere,&
       nonconcRgSp, maxmodellen,kmax,  lenmdlconc,Lnoncon,concmodeloption,rrrrpreftable,maxrampdegree, maxL,mdlchk)
  shellmodelRgmat(st+1:fin,1:3,1:maxmodellen)      = nonconcSpshellmodelRgmat
  lenshellmodel(st+1:fin)                          = nonconcSplenshellmodel
  shellmodelids(st+1:fin,1:4)                      = nonconcSpshellmodelids
  end if

  timecnt = timecnt+1;  call cpu_time(t(timecnt))
  print *, "Time for non-conc Sp modelling", t(timecnt) - t(timecnt-1)

deallocate(nonconcSplenshellmodel,  nonconcRgSp,  nonconcSpshellmodelRgmat, nonconcSpshellmodelids)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***** Add the (SS| terms to our array of models ******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate(RRshellmodelids(1:noRR,1:4), RRlenshellmodel(1:noRR),RRshellmodelRgmat(1:noRR,1:3,1:1))

  call insertRRtomodel(RRshellmodelRgmat,RRlenshellmodel,RRshellmodelids,noRR,RR)
  
  st = fin; fin = st + noRR
  shellmodelRgmat(st+1:fin,1:3,1:1)           = RRshellmodelRgmat;
  shellmodelids(st+1:fin,1:4)                 = RRshellmodelids;
  lenshellmodel(st+1:fin)                     = RRlenshellmodel;
  
deallocate(RRshellmodelRgmat,RRshellmodelids, RRlenshellmodel)

  if (modelRgprint ==1) then; print *, "Shell pair models"
     do i=1,(noRgSs)+(algchange(2))+3*(noRgSp-algchange(2))+noRR
        print *, "bf1 = ", int(shellmodelids(i,1)), " bf2 = ", int(shellmodelids(i,2))
        do j=1,lenshellmodel(i) 
           print "(i3,f20.15,2f4.0,f30.15)", j, shellmodelRgmat(i,1:3,j), shellmodelids(i,4)
        end do;  
     end do; 
  end if;  print *, ""
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***** Convert from shell pair models to basis function pairs models (mostly applicable for concentric Sp shell pairs - the rest is already done) ******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(basisfunmodelRgmat(1:(noRgSs)+3*(noRgSp)+noRR,1:3,1:maxmodellen))
  allocate(bfmodelids(1:(noRgSs)+3*(noRgSp)+noRR,1:4))
  allocate(lenbfmodel(1:(noRgSs)+3*(noRgSp)+noRR))

  call converttobasisfunmodel(basisfunmodelRgmat, lenbfmodel, bfmodelids,shellmodelRgmat,lenshellmodel, shellmodelids,&
       noRR,noRgSs,algchange(2), 3*(noRgSp-algchange(2)),maxmodellen, maxL)

  if (modelRgprint ==1) then;  print *, "Basis Function Pair Models";
     do i=1,(noRgSs)+3*(noRgSp)+noRR
        print *, "bf1 = ", int(bfmodelids(i,1)), " bf2 = ", int(bfmodelids(i,2))
        do j=1,lenbfmodel(i) 
           print "(i3,f25.15,2f4.0,f30.15)", j, basisfunmodelRgmat(i,1:3,j), bfmodelids(i,4)
        end do
     end do;
  end if; print *, ""
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***** Digesting all models for itive quantities to a single model for each basis function pair ******************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

maxlenbf = 10;
do i=1,(noRgSs)+3*(noRgSp)+noRR
   if (lenbfmodel(i) > maxlenbf) then; 
      maxlenbf = lenbfmodel(i)
   end if
end do
maxlenbf = 9*maxlenbf;

allocate(normmodelRgmat(1:nobasisfun,1:nobasisfun,1:3,1:maxlenbf));
! The digestion takes all the models from RS to RI (i.e. remove the 1/(2*sqrt(pi)) factor)
  call modelRgdigestion(normmodelRgmat, lennormmodel ,basisfunmodelRgmat, lenbfmodel, &
       bfmodelids, maxlenbf,nobasisfun,(noRgSs)+3*(noRgSp)+noRR)

  if (modelRgprint ==1) then; print *, "Digested models"
     do bf1=1,nobasisfun;   do bf2=bf1,nobasisfun
        if (lennormmodel(bf1,bf2,1) > 0) then
           print *, "";     print *, "bf1 = ", bf1, "  bf2 = ", bf2; 
        end if
        do j=1,lennormmodel(bf1,bf2,1) 
           print "(i3,f25.15,2f4.0)", j, normmodelRgmat(bf1,bf2,1:3,j)
        end do;   
     end do; end do;  
  end if; print *, ""
  
  call secondsortmodelRG(newsortRg, lenmodelnewsort, totnomodel, normmodelRgmat, lennormmodel, &
     nobasisfun, noatoms, maxmodellen, maxL,maxlenbf, maxrgpa)

print *, "";
  if (modelRgprint ==1) then
     do a = 1,noatoms
        do k=1,maxrgpa;
print *, "";
print *, "*********************************"
           print *, "BC This is bf ids", lenmodelnewsort(1:4,k,a), "atid = ", a
print *, "*********************************"
           do i= lenmodelnewsort(3,k,a)+1, lenmodelnewsort(4,k,a)
              print "(a2,f25.15,3f5.0)", "AA ",newsortRg(1:4,i)
           end do;
 if (lenmodelnewsort(4,k,a) - lenmodelnewsort(3,k,a) .gt. 10) then
           do i= lenmodelnewsort(3,k,a)+1, lenmodelnewsort(4,k,a)
              print "(a2,f25.15,3f5.0)", "BC ",newsortRg(1:4,i)
           end do;
 end if
        end do;  
     end do;
  end if

deallocate(shellmodelRgmat, shellmodelids, lenshellmodel,basisfunmodelRgmat, bfmodelids,lenbfmodel,normmodelRgmat)

  return
end subroutine modelRg
