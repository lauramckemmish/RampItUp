subroutine  converttobasisfunmodel(basisfunmodelRgmat, lenbfmodel, bfmodelids,shellmodelRgmat,lenshellmodel, shellmodelids,&
     noRRall,noSs,noconcSp,noOTHERRg,maxmodellen, maxL)
  implicit none
  integer, intent(in)                   :: noRRall, noSs, noconcSp, noOTHERRg, maxmodellen, maxL
  integer, intent(in)                   :: lenshellmodel(1:noSs+noconcSp+ noOTHERRg +noRRall)
  integer, intent(out)                  :: lenbfmodel(1:noSs+3*noconcSp+ noOTHERRg +noRRall)
  real*8, intent(in)                    :: shellmodelids(1:noSs+noconcSp+noOTHERRg+noRRall,1:4)
  real*8, intent(out)                   :: bfmodelids(1:noSs+3*noconcSp+noOTHERRg+noRRall,1:4)
  real*8, intent(in)                    :: shellmodelRGmat(1:noSs+noconcSp+noOTHERRg+noRRall,1:3,1:maxmodellen)
  real*8, intent(out)                   :: basisfunmodelRGmat(1:noSs+3*noconcSp+ noOTHERRg +noRRall,1:3,1:maxmodellen)
  integer                               :: i, stbf, finbf, stshell, finshell, swit, k, j
  integer                               :: lenpx, lenpy, lenpz,  nKpx, nKpy, nKpz, Ks, m
  integer                               :: Kpx(1:2*maxL+1), Kpy(1:2*maxL+1), Kpz(1:2*maxL+1)
  real*8                                :: cpx(1:2*maxL+1), cpy(1:2*maxL+1), cpz(1:2*maxL+1), betapf
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This function converts the shell pair models to basis function pair models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  basisfunmodelRGmat = 0d0;  lenbfmodel = 0;  bfmodelids = 0;
  
!*******************************!*******************************!*******************************
!!!!!!!!!!!! Converting shell pair to basis function pair for (Ss| models!!!!!!!!!!!!!!!!
!*******************************!*******************************!*******************************
  basisfunmodelRGmat(1:noSs,1:3,1:maxmodellen) = shellmodelRGmat(1:noSs,1:3,1:maxmodellen)
  lenbfmodel(1:noSs)       = lenshellmodel(1:noSs);
  bfmodelids(1:noSs,1:4)   = shellmodelids(1:noSs,1:4);

!*******************************!*******************************!*******************************
!!!!!!!!!!!! Converting shell pair to basis function pair for concentric (Sp|  models!!!!!!!!!!!!!!!!
!*******************************!*******************************!*******************************
  do i=1,noconcSp
     ! This is the s basis function
     bfmodelids(noSs+3*i-2: noSs+3*i,1)   = shellmodelids(noSs+i,1);
     
     ! This is the p basis function 
     bfmodelids(noSs+3*i-2,2)   = shellmodelids(noSs+i,2)
     bfmodelids(noSs+3*i-1,2)   = shellmodelids(noSs+i,2)+1
     bfmodelids(noSs+3*i  ,2)   = shellmodelids(noSs+i,2)+2
     
     ! This is the ramp center
     bfmodelids(noSs+3*i-2:noSs+3*i,3) = shellmodelids(noSs+i,3)
     
     ! This is the normalisation constant
     bfmodelids(noSs+3*i-2:noSs+3*i,4) = shellmodelids(noSs+i,4)
     
     lenbfmodel(noSs+3*i-2:noSs+3*i) = lenshellmodel(noSs+i)
     basisfunmodelRgmat(noSs+3*i-2,  1,1:maxmodellen) = -shellmodelRGmat(noSs+i,1,1:maxmodellen)
     basisfunmodelRgmat(noSs+3*i-2,  2,1:maxmodellen) =  shellmodelRGmat(noSs+i,2,1:maxmodellen)
     !*** The minus takes into account that PI = -x*SI
     basisfunmodelRgmat(noSs+3*i-2,  3,1:maxmodellen) = 4;
     
     basisfunmodelRgmat(noSs+3*i-1,  1,1:maxmodellen) = -shellmodelRGmat(noSs+i,1,1:maxmodellen)
     basisfunmodelRgmat(noSs+3*i-1,  2,1:maxmodellen) =  shellmodelRGmat(noSs+i,2,1:maxmodellen)
     !*** The minus takes into account that PI = -y*SI
     basisfunmodelRgmat(noSs+3*i-1,  3,1:maxmodellen) = 2;
     
     basisfunmodelRgmat(noSs+3*i  ,1:2,1:maxmodellen) =  shellmodelRGmat(noSs+i,1:2,1:maxmodellen)
     !*** No minus here coz                 PI =  z*SI
     basisfunmodelRgmat(noSs+3*i  ,  3,1:maxmodellen) = 3;
     
     ! This is altering the normalisation constant
     bfmodelids(noSs+3*i-2,4) = shellmodelids(noSs+i,4)/sqrt(3d0)
     bfmodelids(noSs+3*i-1,4) = shellmodelids(noSs+i,4)/sqrt(3d0)
     bfmodelids(noSs+3*i  ,4) = shellmodelids(noSs+i,4)/sqrt(3d0)
     
  end do
  
  
  !*******************************!*******************************!*******************************
!!!!!!!!!!!! All others !!!!!!!!!!!!!!!!
  !*******************************!*******************************!*******************************
  stbf    = noSs+3*noconcSp;    finbf    = stbf +noOTHERRg+ noRRall
  stshell = noSs+noconcSp;      finshell = stshell +noOTHERRg+ noRRall
  
  lenbfmodel(stbf+1:finbf) = lenshellmodel(stshell+1:finshell)
  basisfunmodelRGmat(stbf+1:finbf,1:3,1:maxmodellen) = shellmodelRGmat(stshell+1:finshell,1:3,1:maxmodellen)
  bfmodelids(stbf+1:finbf,1:4) = shellmodelids(stshell+1:finshell,1:4)
  
  return
end subroutine converttobasisfunmodel


!***********************************************************************************************
