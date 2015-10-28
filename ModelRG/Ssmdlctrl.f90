subroutine Ssmdlctrl(shellmodelRGmat,lenshellmodel,shellmodelids, no,&
     RGin, maxmodellen,kmax, lenmdl,Ltrial,  concmodeloption,rrrrpreftable,maxrampdegree, Lmax, mdlchk)
  
  implicit none
  real*8, intent(out)            :: shellmodelRGmat(1:no,1:3,1:maxmodellen), shellmodelids(1:no,1:4)
  integer, intent(out)           :: lenshellmodel(1:no)
  integer, intent(in)            :: no,  maxmodellen, kmax, concmodeloption, mdlchk
  integer,intent(in)             :: lenmdl, Ltrial, maxrampdegree, Lmax
  real*8, intent(in)             :: RGin(1:30,1:no), rrrrpreftable(0:maxrampdegree,0:maxrampdegree,0:Lmax)
  real*8                         :: beta, BA2
  integer                        :: K, rn, m, index, L, mdlbasis(1:3,1:lenmdl*(1+Ltrial)**2), fin, st, len, i
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function chooses the ramp basis functions that are to be used in the modelling.
  ! This should probably be modified in the future to suggest the optimal ramp basis functions for a specific
  ! class of Rg shell pairs.
  !
  ! Then, we call the function that finds the coefficients
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
if (Ltrial > Lmax) then; print *, "Ltrial is greater than Lmax"; print *, "Exiting"; call exit; end if;
  lenshellmodel = 0d0; shellmodelids = 0d0; shellmodelRGmat = 0d0;   
   
st = 0; fin =0;
do i=1,no
   st =i-1;
   fin = i;

   beta = RGin(8,i)
   BA2 = Rgin(4,i)**2+Rgin(5,i)**2+Rgin(6,i)**2

   if (beta .lt. 0.1d0) then; 
            len = floor(lenmdl/3d0)
   else if (beta .lt. 0.5) then; 
             len = floor(lenmdl/2d0);
   else if (beta .lt. 1d0) then
            len = floor(lenmdl/1.5d0);
   else
          len = lenmdl
   end if

   if (BA2 .gt. 20d0) then; 
         len = floor(len/1.2d0)
   end if

  if (BA2 .gt. 30d0) then;
         len = floor(len/2d0)
   end if

  if (BA2 .gt. 40d0) then;
         len = floor(len/2d0)
   end if





   !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
   !     MODEL BASIS FUNCTION SELECTION
   !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
   mdlbasis = 0d0; index = 0;
   K = 0; do L=0,Ltrial; do m=-L,L; K = K + 1;  
   do rn=0,ceiling((len-1d0)/((L/2d0)+1d0))
      index = index + 1;
      mdlbasis(1:3,index) = (/rn, K, L/)
   end do;end do;end do


if (beta .gt. 5d0) then
  mdlbasis = 0d0; index = 0;
   K = 0; do L=0,Ltrial; do m=-L,L; K = K + 1;
   do rn=0,ceiling((len-1d0)/((L/3d0)+1d0))
      index = index + 1;
      mdlbasis(1:3,index) = (/rn, K, L/)
   end do;end do;end do
end if

if (beta .lt. 0.5d0) then
  mdlbasis = 0d0; index = 0;
   K = 0; do L=0,2; do m=-L,L; K = K + 1;
   do rn=0,ceiling((len-1d0)/((L/1d0)+1d0))
      index = index + 1;
      mdlbasis(1:3,index) = (/rn, K, L/)
   end do;end do;end do
end if

if (beta .lt. 1.d0) then
  mdlbasis = 0d0; index = 0;
   K = 0; do L=0,3; do m=-L,L; K = K + 1;
   do rn=0,ceiling((len-1d0)/((L/1.5d0)+1d0))
      index = index + 1;
      mdlbasis(1:3,index) = (/rn, K, L/)
   end do;end do;end do
end if



!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     FINDING COEFFICIENTS OF MODEL BASIS FUNCTIONS
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************

   call calcSsmodels(shellmodelRGmat(st+1:fin,1:3,1:maxmodellen),lenshellmodel(st+1:fin),&
        shellmodelids(st+1:fin,1:4), fin-st, RGin(1:30,st+1:fin), maxmodellen,kmax, lenmdl,Ltrial,&
        concmodeloption,rrrrpreftable,maxrampdegree, Lmax,mdlbasis(1:3,1:index),index, mdlchk)
   
end do
end subroutine Ssmdlctrl


!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************



subroutine Spmdlctrl(shellmodelRGmat,lenshellmodel,shellmodelids, no,&
     RGin, maxmodellen,kmax, lenmdl,Ltrial,  concmodeloption,rrrrpreftable,maxrampdegree, Lmax,mdlchk)
  
  implicit none
  real*8, intent(out)            :: shellmodelRGmat(1:3*no,1:3,1:maxmodellen), shellmodelids(1:3*no,1:4)
  integer, intent(out)           :: lenshellmodel(1:3*no)
  integer, intent(in)            :: no,  maxmodellen, kmax, concmodeloption,mdlchk
  integer,intent(in)             :: lenmdl, Ltrial, maxrampdegree, Lmax
  real*8, intent(in)             :: RGin(1:30,1:no), rrrrpreftable(0:maxrampdegree,0:maxrampdegree,0:Lmax)
  real*8                         :: beta, BA2
  integer                        :: mdlbasis(1:3,1:lenmdl*(1+Ltrial)**2), fin, st
  integer                        :: K, i, j, info, n1, K1, L1, n2, K2, L2, rn, m, index, L, len
  
  lenshellmodel = 0d0; shellmodelids = 0d0; shellmodelRGmat = 0d0;
if (Ltrial > Lmax) then; print *, "Ltrial is greater than Lmax"; print *, "Exiting"; call exit; end if;
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function chooses the ramp basis functions that are to be used in the modelling.
!
! Then, we call the function that finds the coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i=1,no

   !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
   !     MODEL BASIS FUNCTION SELECTION
   !********!!!!!!!!!!!!************!!!!!!!!!!!!**************

   st =3*i-3;
   fin = 3*i;

   beta = RGin(8,i)
   BA2 = Rgin(4,i)**2+Rgin(5,i)**2+Rgin(6,i)**2
 
  len = lenmdl
   if (beta .lt. 0.1d0) then; 
            len = floor(lenmdl/3d0)
   else if (beta .lt. 0.5d0) then; 
             len = floor(lenmdl/2d0);
   else if (beta .lt. 1d0) then
            len = floor(lenmdl/1.5d0);
   else
             len = lenmdl
   end if

   if (BA2 .gt. 20d0) then; 
         len = floor(len/1.2d0)
   end if 

  if (BA2 .gt. 40d0) then;
         len = floor(len/1.2d0)
   end if




 mdlbasis = 0d0;  index = 0;
   K = 0; do L=0,Ltrial; do m=-L,L; K = K + 1;  
   do rn=0,ceiling((len-1d0)/((L/2d0)+1d0))
      index = index + 1;
      mdlbasis(1,index) = rn;
      mdlbasis(2,index) = K;
      mdlbasis(3,index) = L;
   end do;end do;end do

if (beta .gt. 5d0) then
  mdlbasis = 0d0; index = 0;
   K = 0; do L=0,Ltrial; do m=-L,L; K = K + 1;
   do rn=0,ceiling((len-1d0)/((L/3d0)+1d0))
      index = index + 1;
      mdlbasis(1:3,index) = (/rn, K, L/)
   end do;end do;end do
end if

if (beta .lt. 0.5d0) then
  mdlbasis = 0d0; index = 0;
   K = 0; do L=0,2; do m=-L,L; K = K + 1;
   do rn=0,ceiling((len-1d0)/((L/1d0)+1d0))
      index = index + 1;
      mdlbasis(1:3,index) = (/rn, K, L/)
   end do;end do;end do
end if

if (beta .lt. 1d0) then
  mdlbasis = 0d0; index = 0;
   K = 0; do L=0,3; do m=-L,L; K = K + 1;
   do rn=0,ceiling((len-1d0)/((L/1.5d0)+1d0))
      index = index + 1;
      mdlbasis(1:3,index) = (/rn, K, L/)
   end do;end do;end do
end if



!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     FINDING COEFFICIENTS OF MODEL BASIS FUNCTIONS
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
   
   call calcSpmodels(shellmodelRGmat(st+1:fin,1:3,1:maxmodellen),lenshellmodel(st+1:fin),&
        shellmodelids(st+1:fin,1:4),1, RGin(1:30,i), maxmodellen,kmax, lenmdl,Ltrial,&
        concmodeloption,rrrrpreftable,maxrampdegree, Lmax,mdlbasis(1:3,1:index),index,mdlchk)
   
   end do

 end subroutine Spmdlctrl
