subroutine calcSpmodels(Spmdl,lenshmdl1,shmdl1ids, no, RGinr, maxmodellen,kmax, lenmdl,Ltrial,&
  concmodeloption,rrrrpreftable,maxrampdegree, Lmax, mdlbasis, index,mdlchk)
 implicit none
 real*8, intent(out)            :: Spmdl(1:3,1:3,1:maxmodellen), shmdl1ids(1:3,1:4)
 integer, intent(out)           :: lenshmdl1(1:3)
 integer, intent(in)            :: no,  maxmodellen, kmax, concmodeloption,mdlchk, index
 integer,intent(in)             :: lenmdl, Ltrial, maxrampdegree, Lmax,mdlbasis(1:3,1:index)
 real*8, intent(in)             :: RGinr(1:30,1:1)
 real*8, intent(in)             :: rrrrpreftable(0:maxrampdegree,0:maxrampdegree,0:Lmax)
 integer                        :: K,ntot(1:1), i, j, info, n1, K1, L1, n2, K2, L2, rn, m, L,ai
 integer                        :: nLCart, nLPure, cntx, cnty, cntz, kcurrent, st, fin, swit, rank
 real*8                        :: RGr(1:30,1:1), OlrgRs1(1:1,1:4), n(1:1)
 real*8                        :: pi, Mat(1:index,1:index), Matfixed(1:index,1:index)
 real*8                        :: Vecx(1:index,1:1), DummyArray(1:index),Matnew(1:index,1:index)
 real*8                        :: Vecy(1:index,1:1),Vecz(1:index,1:1), Vecnew(1:index,1:1)
 real*8                        :: OLRpx(1:1),OLRpy(1:1), OLRpz(1:1), x
 real*8                        :: WORK(1:10*index),qpDummyMat(1:index,1:index)

pi = dacos(-1d0); lenshmdl1 = 0d0; shmdl1ids = 0d0; Spmdl = 0d0;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the coefficients of the model for Sp in terms of a sum of ramps of different degrees
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n(1:1) = RGinr(7,1:1);    RGr = real(RGinr,16);

!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     MATRIX CONSTRUCTION
! Calculate overlap of model basis functions with each other 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************

do i=1, index; 
   n1 = mdlbasis(1,i);   K1 = mdlbasis(2,i);   L1 = mdlbasis(3,i);
do j=1,index; 
   n2 = mdlbasis(1,j);   K2 = mdlbasis(2,j);   L2 = mdlbasis(3,j);
 if (K1 == K2) then;
!*** This could be precomputed (as all are constants) **** 
      Matfixed(i,j) =  exp(dlgama(3d0+2d0*l1)+dlgama(1d0+n1+n2)-dlgama(4d0+2d0*l1+n1+n2));
 end if
end do; end do; 


!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     VECTOR CONSTRUCTION
! Calculate overlap of model basis functions with gaussian 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
nLCart = (Lmax+1)*(Lmax+2)*(Lmax+3)/6;     nLPure = (1+Lmax)**2; 
do i = 1, index; 
   n1 = mdlbasis(1,i);   K1 = mdlbasis(2,i);   L1 = mdlbasis(3,i);
   RGr(7,1:1) = n1;
   call calcSpol(OLRpx, OLRpy, OLRpz,Rgr,no,K1,L1, kmax)
   Vecx(i,1) = OLRpx(1);         Vecy(i,1) = OLRpy(1);         Vecz(i,1) = OLRpz(1);  
end do


!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     LINEAR SOLVE
! Calculate coefficients by solving Mc = v 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************

  cntx = 0; cnty = 0; cntz = 0;   st = 0; fin =0;
do kcurrent = 1,1;
  st = fin; swit = 0;
  do k=st+1,index; if (swit == 0 .and. mdlbasis(2,k) >  kcurrent) then;  swit = 1; fin = k-1; end if;  end do

     Mat(st+1:fin,st+1:fin) = Matfixed(st+1:fin,st+1:fin);
     do ai=st+1,fin
         n1 = mdlbasis(1,ai)+n(1)
         Mat(ai,fin+1) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
         Mat(fin+1,ai) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
     end do
    Vecnew = Vecx
   RGr(7,1:1) = n(1);
   call calcSpol(OLRpx, OLRpy, OLRpz,Rgr,1,1,0, kmax)

    Vecnew(fin+1,1) = OLRpx(1)

     call dgesv(fin-st+1,1,Mat(st+1:fin+1,st+1:fin+1),fin-st+1,DummyArray,Vecnew(st+1:fin+1,1), fin-st+1, info)
     !call dgesv(fin-st,1,Mat(st+1:fin,st+1:fin),fin-st,DummyArray,Vecnew(st+1:fin,1), fin-st, info)

print *, "LAMBDA = ", Vecnew(fin+1,1)


  do k=st+1,fin
    if (Abs(Vecnew(k,1)) > 1d-10) then
       cntx = cntx + 1;
       Spmdl(1,1,cntx) = dble(Vecnew(k,1));          ! Coefficient of ramp
       Spmdl(1,2,cntx) = mdlbasis(1,k)+n(1);       ! Degree of model ramp
       Spmdl(1,3,cntx) = mdlbasis(2,k);            ! Angular momentum of model ramp in K form
    end if
  end do

     Mat(st+1:fin,st+1:fin) = Matfixed(st+1:fin,st+1:fin);
     do ai=st+1,fin
         n1 = mdlbasis(1,ai)+n(1)
         Mat(ai,fin+1) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
         Mat(fin+1,ai) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
     end do
    Vecnew = Vecy
    Vecnew(fin+1,1) = OLRpy(1)

     call dgesv(fin-st+1,1,Mat(st+1:fin+1,st+1:fin+1),fin-st+1,DummyArray,Vecnew(st+1:fin+1,1), fin-st+1, info)

!     call dgesv(fin-st,1,Mat(st+1:fin,st+1:fin),fin-st,DummyArray,Vecnew(st+1:fin,1), fin-st, info)

print *, "LAMBDA = ", Vecnew(fin+1,1)

  do k=st+1,fin
    if (Abs(Vecnew(k,1)) > 1d-10) then
       cnty = cnty + 1;
       Spmdl(2,1,cnty) = dble(Vecnew(k,1));          ! Coefficient of ramp
       Spmdl(2,2,cnty) = mdlbasis(1,k)+n(1);       ! Degree of model ramp
       Spmdl(2,3,cnty) = mdlbasis(2,k);            ! Angular momentum of model ramp in K form
    end if
  end do

     Mat(st+1:fin,st+1:fin) = Matfixed(st+1:fin,st+1:fin);
     do ai=st+1,fin
         n1 = mdlbasis(1,ai)+n(1)
         Mat(ai,fin+1) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
         Mat(fin+1,ai) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
     end do
    Vecnew = Vecz 
    Vecnew(fin+1,1) = OLRpz(1)


     call dgesv(fin-st+1,1,Mat(st+1:fin+1,st+1:fin+1),fin-st+1,DummyArray,Vecnew(st+1:fin+1,1), fin-st+1, info)

    ! call dgesv(fin-st,1,Mat(st+1:fin,st+1:fin),fin-st,DummyArray,Vecnew(st+1:fin,1), fin-st, info)

print *, "LAMBDA = ", Vecnew(fin+1,1)

x = 0;
  do k=st+1,fin
    if (Abs(Vecnew(k,1)) > 1d-10) then
       cntz = cntz + 1;
       Spmdl(3,1,cntz) = dble(Vecnew(k,1));            ! Coefficient of ramp
 x = x + Vecnew(k,1)*Mat(k,fin+1) 
      Spmdl(3,2,cntz) = mdlbasis(1,k)+n(1);         ! Degree of model ramp
       Spmdl(3,3,cntz) = mdlbasis(2,k);              ! Angular momentum of model ramp in K form
    end if
!print *, "AAA IMPORTANT = ", x, OLRpz(1)
  end do
end do



do kcurrent = 2,(1+Lmax)**2;
  st = fin; swit = 0;
  do k=st+1,index; if (swit == 0 .and. mdlbasis(2,k) >  kcurrent) then;  swit = 1; fin = k-1; end if;  end do
  if (swit == 0) then;      fin = index;  end if
  if (fin-st .ne. 0) then;
     Mat(st+1:fin,st+1:fin) = Matfixed(st+1:fin,st+1:fin);  
     call dgesv(fin-st,1,Mat(st+1:fin,st+1:fin),fin-st,DummyArray,Vecx(st+1:fin,1), fin-st, info)
     Mat(st+1:fin,st+1:fin) = Matfixed(st+1:fin,st+1:fin);
     call dgesv(fin-st,1,Mat(st+1:fin,st+1:fin),fin-st,DummyArray,Vecy(st+1:fin,1), fin-st, info)
     Mat(st+1:fin,st+1:fin) = Matfixed(st+1:fin,st+1:fin); 
     call dgesv(fin-st,1,Mat(st+1:fin,st+1:fin),fin-st,DummyArray,Vecz(st+1:fin,1), fin-st, info)
  end if

  do k=st+1,fin
    if (Abs(Vecx(k,1)) > 1d-10) then
       cntx = cntx + 1;
       Spmdl(1,1,cntx) = dble(Vecx(k,1));          ! Coefficient of ramp    
       Spmdl(1,2,cntx) = mdlbasis(1,k)+n(1);       ! Degree of model ramp
       Spmdl(1,3,cntx) = mdlbasis(2,k);            ! Angular momentum of model ramp in K form
    end if
    if (Abs(Vecy(k,1)) > 1d-10) then
       cnty = cnty + 1;
       Spmdl(2,1,cnty) = dble(Vecy(k,1));          ! Coefficient of ramp    
       Spmdl(2,2,cnty) = mdlbasis(1,k)+n(1);       ! Degree of model ramp
       Spmdl(2,3,cnty) = mdlbasis(2,k);            ! Angular momentum of model ramp in K form
    end if
    if (Abs(Vecz(k,1)) > 1d-10) then
       cntz = cntz + 1;
       Spmdl(3,1,cntz) = dble(Vecz(k,1));            ! Coefficient of ramp    
       Spmdl(3,2,cntz) = mdlbasis(1,k)+n(1);         ! Degree of model ramp
       Spmdl(3,3,cntz) = mdlbasis(2,k);              ! Angular momentum of model ramp in K form
    end if
  end do
end do
lenshmdl1(1:3) = (/cntx,cnty,cntz/);

  shmdl1ids(1:3,1) = RGr(1,1);    shmdl1ids(1  ,2) = RGr(2,1);  shmdl1ids(2  ,2) = RGr(2,1)+1;
  shmdl1ids(3  ,2) = RGr(2,1)+2;  shmdl1ids(1:3,3) = RGr(3,1);  shmdl1ids(1:3,4) = RGr(9,1);

return
end subroutine calcSpmodels

!************************************************************************************


subroutine calcSpol(OLRpx,OLRpy,OLRpz,RGSp,no,K,L, kmax) 
 implicit none
 integer, intent(in)      :: no, K, L, kmax
 real*8, intent(in)      :: RGSp(1:30,1:1)
 real*8, intent(out)     :: OLRpx(1:1), OLRpy(1:1), OLRpz(1:1)
 real*8                  :: Rs(1:1,(2+L)*(3+L)*(4+L)/6,0:L+1)
 integer                  :: nLCart, nLPure, i, n(1:1)
 real*8                  :: CartRpx(1:1,1:(1+L)*(2+L)*(3+L)/6,0:L)
 real*8                  :: CartRpy(1:1,1:(1+L)*(2+L)*(3+L)/6,0:L)
 real*8                  :: CartRpz(1:1,1:(1+L)*(2+L)*(3+L)/6,0:L)
 real*8                  :: PureDef1Rpx(1:1,1:(1+L)**2), PureDef1Rpy(1:1,1:(1+L)**2)
 real*8                  :: PureDef1Rpz(1:1,1:(1+L)**2)
 real*8                  :: Fund(1:1,0:L+1), beta(1:1), BAx(1:1), BAy(1:1), BAz(1:1)
 real*8                  :: T(1:1), h(1:1,0:kmax,0:L+1)

 nLCart = (1+L)*(2+L)*(3+L)/6;    nLPure = (1+L)**2

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!  This subroutine calculates pure overlap [R|p] integral over all space
 !!  int( (1-r)^n r^{l} RY_{l,m}[theta,ph] e^{-\beta |r-P|^2} ) r^2 sin[theta] dr dtheta dphi
 !!
 !!  i.e. overlap between ramp and gaussian 
 !! 
 !! where RY_{l,m} =         Y_{l,m}, m ==0
 !!                = 1/sqrt(2) * (Y_{l, m} + (-1)^m Y_{l,-m}) , m > 0
 !!                = I/sqrt(2) * (Y_{l,-m} + (-1)^m Y_{l, m}) , m < 0
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! The methodology here is to fiRs1t calculate integral over all space
 !!  int( (1-r)^n x^{i} y^{j} z^{k} e^{-\beta |r-P|^2} ) r^2 sin[theta] dr dtheta dphi
 !! 
 !! Then transform these to the desired integrals. 

 !! In all cases, the beta is small or P = 0 (i.e. not concentric and large beta), so we can use general method (truncating infinite series)

   beta =RGSp(8,1:1) ;   BAx = RGSp(4,1:1) ;   BAy = RGSp(5,1:1) ;   BAz = RGSp(6,1:1) ;
   T = beta * (BAx*BAx + BAy*BAy + BAz*BAz);     n   = RgSp(7,1:1) ;

   call NEWcalcmphk(h,kmax,L+1,T,no)                                   !  Compute h_k(T) and their derivatives
   call NEWcalcintmpSs(Fund,h,kmax,L+1,n,beta,no)                      !  Compute [S|s]^(m) integrals
   call NEWcalcmpCartRs(Rs,Fund,(2+L)*(3+L)*(4+L)/6,L+1,beta,BAx,BAy,BAz,no)        !  Compute d([R|s])/dB_i integrals

   call NEWCartRs2CartRp(CartRpx,CartRpy,CartRpz,Rs,L,nLCart,BAx,BAy,BAz,no);
   
   call NEWCart2PureDef1(PureDef1Rpx,CartRpx,no,nLCart,nLPure, L);
   OLRpx = PureDef1Rpx(1:1,K);

   call NEWCart2PureDef1(PureDef1Rpy,CartRpy,no,nLCart,nLPure, L);
   OLRpy = PureDef1Rpy(1:1,K);

   call NEWCart2PureDef1(PureDef1Rpz,CartRpz,no,nLCart,nLPure, L);
   OLRpz = PureDef1Rpz(1:1,K);

return
end subroutine calcSpol



!**************************************************************************************************

subroutine NEWCartRs2CartRp(CartRpx,CartRpy,CartRpz,Rs,L,nLCart,BAx,BAy,BAz,no)
 implicit none
 integer, intent(in) :: nLCart, L, no
 real*8, intent(in)  :: BAx(1:1), BAy(1:1), BAz(1:1), Rs(1:1,(2+L)*(3+L)*(4+L)/6,0:L+1)
 real*8, intent(out) :: CartRpx(1:1,1:nLCart), CartRpy(1:1,1:nLCart), CartRpz(1:1,1:nLCart)
 integer i

!************** Application of the Horizontal Recurrence relation to move angular momentum from ramp to Gaussian *****************

  CartRpx(1:1,1) = -BAx*Rs(1:1,1,0)+ Rs(1:1,2,0)      ! [S|px]
  CartRpy(1:1,1) = -BAy*Rs(1:1,1,0)+ Rs(1:1,3,0)      ! [S|py]
  CartRpz(1:1,1) = -BAz*Rs(1:1,1,0)+ Rs(1:1,4,0)      ! [S|pz]

if (L < 1) then; return; end if
  CartRpx(1:1,2) = -BAx*Rs(1:1,2,0)+ Rs(1:1,5,0)             
  CartRpy(1:1,2) = -BAy*Rs(1:1,2,0)+ Rs(1:1,6,0)             
  CartRpz(1:1,2) = -BAz*Rs(1:1,2,0)+ Rs(1:1,8,0)                

  CartRpx(1:1,3) = -BAx*Rs(1:1,3,0)+ Rs(1:1,6,0)             
  CartRpy(1:1,3) = -BAy*Rs(1:1,3,0)+ Rs(1:1,7,0)             
  CartRpz(1:1,3) = -BAz*Rs(1:1,3,0)+ Rs(1:1,9,0)                

  CartRpx(1:1,4) = -BAx*Rs(1:1,4,0)+ Rs(1:1,8,0)             
  CartRpy(1:1,4) = -BAy*Rs(1:1,4,0)+ Rs(1:1,9,0)             
  CartRpz(1:1,4) = -BAz*Rs(1:1,4,0)+ Rs(1:1,10,0)     


if (L < 2) then; return; end if
  CartRpx(1:1,5) = -BAx*Rs(1:1,5,0)+ Rs(1:1,11,0)             
  CartRpy(1:1,5) = -BAy*Rs(1:1,5,0)+ Rs(1:1,12,0)             
  CartRpz(1:1,5) = -BAz*Rs(1:1,5,0)+ Rs(1:1,15,0)                

  CartRpx(1:1,6) = -BAx*Rs(1:1,6,0)+ Rs(1:1,12,0)             
  CartRpy(1:1,6) = -BAy*Rs(1:1,6,0)+ Rs(1:1,13,0)             
  CartRpz(1:1,6) = -BAz*Rs(1:1,6,0)+ Rs(1:1,16,0)                

  CartRpx(1:1,7) = -BAx*Rs(1:1,7,0)+ Rs(1:1,13,0)             
  CartRpy(1:1,7) = -BAy*Rs(1:1,7,0)+ Rs(1:1,14,0)             
  CartRpz(1:1,7) = -BAz*Rs(1:1,7,0)+ Rs(1:1,17,0)                

  CartRpx(1:1,8) = -BAx*Rs(1:1,8,0)+ Rs(1:1,15,0)             
  CartRpy(1:1,8) = -BAy*Rs(1:1,8,0)+ Rs(1:1,16,0)             
  CartRpz(1:1,8) = -BAz*Rs(1:1,8,0)+ Rs(1:1,18,0)  

  CartRpx(1:1,9) = -BAx*Rs(1:1,9,0)+ Rs(1:1,16,0)             
  CartRpy(1:1,9) = -BAy*Rs(1:1,9,0)+ Rs(1:1,17,0)             
  CartRpz(1:1,9) = -BAz*Rs(1:1,9,0)+ Rs(1:1,19,0)                
              
  CartRpx(1:1,10) = -BAx*Rs(1:1,10,0)+ Rs(1:1,18,0)             
  CartRpy(1:1,10) = -BAy*Rs(1:1,10,0)+ Rs(1:1,19,0)             
  CartRpz(1:1,10) = -BAz*Rs(1:1,10,0)+ Rs(1:1,20,0)                
              


! Octupoles of Sp basis function pairs

if (L < 3) then; return; end if
  CartRpx(1:1,11) = -BAx*Rs(1:1,11,0)+ Rs(1:1,21,0)             
  CartRpy(1:1,11) = -BAy*Rs(1:1,11,0)+ Rs(1:1,22,0)             
  CartRpz(1:1,11) = -BAz*Rs(1:1,11,0)+ Rs(1:1,26,0) 

  CartRpx(1:1,12) = -BAx*Rs(1:1,12,0)+ Rs(1:1,22,0)             
  CartRpy(1:1,12) = -BAy*Rs(1:1,12,0)+ Rs(1:1,23,0)             
  CartRpz(1:1,12) = -BAz*Rs(1:1,12,0)+ Rs(1:1,27,0)                
              
  CartRpx(1:1,13) = -BAx*Rs(1:1,13,0)+ Rs(1:1,23,0)             
  CartRpy(1:1,13) = -BAy*Rs(1:1,13,0)+ Rs(1:1,24,0)             
  CartRpz(1:1,13) = -BAz*Rs(1:1,13,0)+ Rs(1:1,28,0) 

  CartRpx(1:1,14) = -BAx*Rs(1:1,14,0)+ Rs(1:1,24,0)             
  CartRpy(1:1,14) = -BAy*Rs(1:1,14,0)+ Rs(1:1,25,0)             
  CartRpz(1:1,14) = -BAz*Rs(1:1,14,0)+ Rs(1:1,29,0) 

  CartRpx(1:1,15) = -BAx*Rs(1:1,15,0)+ Rs(1:1,26,0)             
  CartRpy(1:1,15) = -BAy*Rs(1:1,15,0)+ Rs(1:1,27,0)             
  CartRpz(1:1,15) = -BAz*Rs(1:1,15,0)+ Rs(1:1,30,0) 

  CartRpx(1:1,16) = -BAx*Rs(1:1,16,0)+ Rs(1:1,27,0)             
  CartRpy(1:1,16) = -BAy*Rs(1:1,16,0)+ Rs(1:1,28,0)             
  CartRpz(1:1,16) = -BAz*Rs(1:1,16,0)+ Rs(1:1,31,0) 

  CartRpx(1:1,17) = -BAx*Rs(1:1,17,0)+ Rs(1:1,28,0)             
  CartRpy(1:1,17) = -BAy*Rs(1:1,17,0)+ Rs(1:1,29,0)             
  CartRpz(1:1,17) = -BAz*Rs(1:1,17,0)+ Rs(1:1,32,0) 

  CartRpx(1:1,18) = -BAx*Rs(1:1,18,0)+ Rs(1:1,30,0)             
  CartRpy(1:1,18) = -BAy*Rs(1:1,18,0)+ Rs(1:1,31,0)             
  CartRpz(1:1,18) = -BAz*Rs(1:1,18,0)+ Rs(1:1,33,0) 

  CartRpx(1:1,19) = -BAx*Rs(1:1,19,0)+ Rs(1:1,31,0)             
  CartRpy(1:1,19) = -BAy*Rs(1:1,19,0)+ Rs(1:1,32,0)             
  CartRpz(1:1,19) = -BAz*Rs(1:1,19,0)+ Rs(1:1,34,0)

  CartRpx(1:1,20) = -BAx*Rs(1:1,20,0)+ Rs(1:1,33,0)             
  CartRpy(1:1,20) = -BAy*Rs(1:1,20,0)+ Rs(1:1,34,0)             
  CartRpz(1:1,20) = -BAz*Rs(1:1,20,0)+ Rs(1:1,35,0) 


! ********** L = 4 ***********

if (L < 4) then; return; end if
  CartRpx(1:1,21) = -BAx*Rs(1:1,21,0)+ Rs(1:1,36,0)             
  CartRpy(1:1,21) = -BAy*Rs(1:1,21,0)+ Rs(1:1,37,0)             
  CartRpz(1:1,21) = -BAz*Rs(1:1,21,0)+ Rs(1:1,42,0) 

  CartRpx(1:1,22) = -BAx*Rs(1:1,22,0)+ Rs(1:1,37,0)             
  CartRpy(1:1,22) = -BAy*Rs(1:1,22,0)+ Rs(1:1,38,0)             
  CartRpz(1:1,22) = -BAz*Rs(1:1,22,0)+ Rs(1:1,43,0) 

  CartRpx(1:1,23) = -BAx*Rs(1:1,23,0)+ Rs(1:1,38,0)             
  CartRpy(1:1,23) = -BAy*Rs(1:1,23,0)+ Rs(1:1,39,0)             
  CartRpz(1:1,23) = -BAz*Rs(1:1,23,0)+ Rs(1:1,44,0) 

  CartRpx(1:1,24) = -BAx*Rs(1:1,24,0)+ Rs(1:1,39,0)             
  CartRpy(1:1,24) = -BAy*Rs(1:1,24,0)+ Rs(1:1,40,0)             
  CartRpz(1:1,24) = -BAz*Rs(1:1,24,0)+ Rs(1:1,45,0) 

  CartRpx(1:1,25) = -BAx*Rs(1:1,25,0)+ Rs(1:1,40,0)             
  CartRpy(1:1,25) = -BAy*Rs(1:1,25,0)+ Rs(1:1,41,0)             
  CartRpz(1:1,25) = -BAz*Rs(1:1,25,0)+ Rs(1:1,46,0) 

  CartRpx(1:1,26) = -BAx*Rs(1:1,26,0)+ Rs(1:1,42,0)             
  CartRpy(1:1,26) = -BAy*Rs(1:1,26,0)+ Rs(1:1,43,0)             
  CartRpz(1:1,26) = -BAz*Rs(1:1,26,0)+ Rs(1:1,47,0) 

  CartRpx(1:1,27) = -BAx*Rs(1:1,27,0)+ Rs(1:1,43,0)             
  CartRpy(1:1,27) = -BAy*Rs(1:1,27,0)+ Rs(1:1,44,0)             
  CartRpz(1:1,27) = -BAz*Rs(1:1,27,0)+ Rs(1:1,48,0) 

  CartRpx(1:1,28) = -BAx*Rs(1:1,28,0)+ Rs(1:1,44,0)             
  CartRpy(1:1,28) = -BAy*Rs(1:1,28,0)+ Rs(1:1,45,0)             
  CartRpz(1:1,28) = -BAz*Rs(1:1,28,0)+ Rs(1:1,49,0) 

  CartRpx(1:1,29) = -BAx*Rs(1:1,29,0)+ Rs(1:1,45,0)             
  CartRpy(1:1,29) = -BAy*Rs(1:1,29,0)+ Rs(1:1,46,0)             
  CartRpz(1:1,29) = -BAz*Rs(1:1,29,0)+ Rs(1:1,50,0) 

  CartRpx(1:1,30) = -BAx*Rs(1:1,30,0)+ Rs(1:1,47,0)             
  CartRpy(1:1,30) = -BAy*Rs(1:1,30,0)+ Rs(1:1,48,0)             
  CartRpz(1:1,30) = -BAz*Rs(1:1,30,0)+ Rs(1:1,51,0) 

  CartRpx(1:1,31) = -BAx*Rs(1:1,31,0)+ Rs(1:1,48,0)             
  CartRpy(1:1,31) = -BAy*Rs(1:1,31,0)+ Rs(1:1,49,0)             
  CartRpz(1:1,31) = -BAz*Rs(1:1,31,0)+ Rs(1:1,52,0) 

  CartRpx(1:1,32) = -BAx*Rs(1:1,32,0)+ Rs(1:1,49,0)             
  CartRpy(1:1,32) = -BAy*Rs(1:1,32,0)+ Rs(1:1,50,0)             
  CartRpz(1:1,32) = -BAz*Rs(1:1,32,0)+ Rs(1:1,53,0) 

  CartRpx(1:1,33) = -BAx*Rs(1:1,33,0)+ Rs(1:1,51,0)             
  CartRpy(1:1,33) = -BAy*Rs(1:1,33,0)+ Rs(1:1,52,0)             
  CartRpz(1:1,33) = -BAz*Rs(1:1,33,0)+ Rs(1:1,54,0) 

  CartRpx(1:1,34) = -BAx*Rs(1:1,34,0)+ Rs(1:1,52,0)             
  CartRpy(1:1,34) = -BAy*Rs(1:1,34,0)+ Rs(1:1,53,0)             
  CartRpz(1:1,34) = -BAz*Rs(1:1,34,0)+ Rs(1:1,55,0) 

  CartRpx(1:1,35) = -BAx*Rs(1:1,35,0)+ Rs(1:1,54,0)             
  CartRpy(1:1,35) = -BAy*Rs(1:1,35,0)+ Rs(1:1,55,0)             
  CartRpz(1:1,35) = -BAz*Rs(1:1,35,0)+ Rs(1:1,56,0) 



! ********** L = 5 ***********

if (L < 5) then; return; end if
  CartRpx(1:1,36) = -BAx*Rs(1:1,36,0)+ Rs(1:1,57,0)             
  CartRpy(1:1,36) = -BAy*Rs(1:1,36,0)+ Rs(1:1,58,0)             
  CartRpz(1:1,36) = -BAz*Rs(1:1,36,0)+ Rs(1:1,64,0) 

  CartRpx(1:1,37) = -BAx*Rs(1:1,37,0)+ Rs(1:1,58,0)             
  CartRpy(1:1,37) = -BAy*Rs(1:1,37,0)+ Rs(1:1,59,0)             
  CartRpz(1:1,37) = -BAz*Rs(1:1,37,0)+ Rs(1:1,65,0)

  CartRpx(1:1,38) = -BAx*Rs(1:1,38,0)+ Rs(1:1,59,0)             
  CartRpy(1:1,38) = -BAy*Rs(1:1,38,0)+ Rs(1:1,60,0)             
  CartRpz(1:1,38) = -BAz*Rs(1:1,38,0)+ Rs(1:1,66,0)

  CartRpx(1:1,39) = -BAx*Rs(1:1,39,0)+ Rs(1:1,60,0)             
  CartRpy(1:1,39) = -BAy*Rs(1:1,39,0)+ Rs(1:1,61,0)             
  CartRpz(1:1,39) = -BAz*Rs(1:1,39,0)+ Rs(1:1,67,0)

  CartRpx(1:1,40) = -BAx*Rs(1:1,40,0)+ Rs(1:1,61,0)             
  CartRpy(1:1,40) = -BAy*Rs(1:1,40,0)+ Rs(1:1,62,0)             
  CartRpz(1:1,40) = -BAz*Rs(1:1,40,0)+ Rs(1:1,68,0)

  CartRpx(1:1,41) = -BAx*Rs(1:1,41,0)+ Rs(1:1,62,0)             
  CartRpy(1:1,41) = -BAy*Rs(1:1,41,0)+ Rs(1:1,63,0)             
  CartRpz(1:1,41) = -BAz*Rs(1:1,41,0)+ Rs(1:1,69,0)

  CartRpx(1:1,42) = -BAx*Rs(1:1,42,0)+ Rs(1:1,64,0)             
  CartRpy(1:1,42) = -BAy*Rs(1:1,42,0)+ Rs(1:1,65,0)             
  CartRpz(1:1,42) = -BAz*Rs(1:1,42,0)+ Rs(1:1,70,0)

  CartRpx(1:1,43) = -BAx*Rs(1:1,43,0)+ Rs(1:1,65,0)             
  CartRpy(1:1,43) = -BAy*Rs(1:1,43,0)+ Rs(1:1,66,0)             
  CartRpz(1:1,43) = -BAz*Rs(1:1,43,0)+ Rs(1:1,71,0)

  CartRpx(1:1,44) = -BAx*Rs(1:1,44,0)+ Rs(1:1,66,0)             
  CartRpy(1:1,44) = -BAy*Rs(1:1,44,0)+ Rs(1:1,67,0)             
  CartRpz(1:1,44) = -BAz*Rs(1:1,44,0)+ Rs(1:1,72,0)

  CartRpx(1:1,45) = -BAx*Rs(1:1,45,0)+ Rs(1:1,67,0)             
  CartRpy(1:1,45) = -BAy*Rs(1:1,45,0)+ Rs(1:1,68,0)             
  CartRpz(1:1,45) = -BAz*Rs(1:1,45,0)+ Rs(1:1,73,0)

  CartRpx(1:1,46) = -BAx*Rs(1:1,46,0)+ Rs(1:1,68,0)             
  CartRpy(1:1,46) = -BAy*Rs(1:1,46,0)+ Rs(1:1,69,0)             
  CartRpz(1:1,46) = -BAz*Rs(1:1,46,0)+ Rs(1:1,74,0)

  CartRpx(1:1,47) = -BAx*Rs(1:1,47,0)+ Rs(1:1,70,0)             
  CartRpy(1:1,47) = -BAy*Rs(1:1,47,0)+ Rs(1:1,71,0)             
  CartRpz(1:1,47) = -BAz*Rs(1:1,47,0)+ Rs(1:1,75,0)

  CartRpx(1:1,48) = -BAx*Rs(1:1,48,0)+ Rs(1:1,71,0)             
  CartRpy(1:1,48) = -BAy*Rs(1:1,48,0)+ Rs(1:1,72,0)             
  CartRpz(1:1,48) = -BAz*Rs(1:1,48,0)+ Rs(1:1,76,0)

  CartRpx(1:1,49) = -BAx*Rs(1:1,49,0)+ Rs(1:1,72,0)             
  CartRpy(1:1,49) = -BAy*Rs(1:1,49,0)+ Rs(1:1,73,0)             
  CartRpz(1:1,49) = -BAz*Rs(1:1,49,0)+ Rs(1:1,77,0)

  CartRpx(1:1,50) = -BAx*Rs(1:1,50,0)+ Rs(1:1,73,0)             
  CartRpy(1:1,50) = -BAy*Rs(1:1,50,0)+ Rs(1:1,74,0)             
  CartRpz(1:1,50) = -BAz*Rs(1:1,50,0)+ Rs(1:1,78,0)

  CartRpx(1:1,51) = -BAx*Rs(1:1,51,0)+ Rs(1:1,75,0)             
  CartRpy(1:1,51) = -BAy*Rs(1:1,51,0)+ Rs(1:1,76,0)             
  CartRpz(1:1,51) = -BAz*Rs(1:1,51,0)+ Rs(1:1,79,0)

  CartRpx(1:1,52) = -BAx*Rs(1:1,52,0)+ Rs(1:1,76,0)             
  CartRpy(1:1,52) = -BAy*Rs(1:1,52,0)+ Rs(1:1,77,0)             
  CartRpz(1:1,52) = -BAz*Rs(1:1,52,0)+ Rs(1:1,80,0)

  CartRpx(1:1,53) = -BAx*Rs(1:1,53,0)+ Rs(1:1,77,0)             
  CartRpy(1:1,53) = -BAy*Rs(1:1,53,0)+ Rs(1:1,78,0)             
  CartRpz(1:1,53) = -BAz*Rs(1:1,53,0)+ Rs(1:1,81,0)

  CartRpx(1:1,54) = -BAx*Rs(1:1,54,0)+ Rs(1:1,79,0)             
  CartRpy(1:1,54) = -BAy*Rs(1:1,54,0)+ Rs(1:1,80,0)             
  CartRpz(1:1,54) = -BAz*Rs(1:1,54,0)+ Rs(1:1,82,0)

  CartRpx(1:1,55) = -BAx*Rs(1:1,55,0)+ Rs(1:1,80,0)             
  CartRpy(1:1,55) = -BAy*Rs(1:1,55,0)+ Rs(1:1,81,0)             
  CartRpz(1:1,55) = -BAz*Rs(1:1,55,0)+ Rs(1:1,83,0)

  CartRpx(1:1,56) = -BAx*Rs(1:1,56,0)+ Rs(1:1,82,0)             
  CartRpy(1:1,56) = -BAy*Rs(1:1,56,0)+ Rs(1:1,83,0)             
  CartRpz(1:1,56) = -BAz*Rs(1:1,56,0)+ Rs(1:1,84,0)



! ********** L = 6 ***********

if (L < 6) then; return; end if

  CartRpx(1:1,57) = -BAx*Rs(1:1,57,0)+ Rs(1:1,85,0)             
  CartRpy(1:1,57) = -BAy*Rs(1:1,57,0)+ Rs(1:1,86,0)             
  CartRpz(1:1,57) = -BAz*Rs(1:1,57,0)+ Rs(1:1,93,0)

  CartRpx(1:1,58) = -BAx*Rs(1:1,58,0)+ Rs(1:1,86,0)             
  CartRpy(1:1,58) = -BAy*Rs(1:1,58,0)+ Rs(1:1,87,0)             
  CartRpz(1:1,58) = -BAz*Rs(1:1,58,0)+ Rs(1:1,94,0)

  CartRpx(1:1,59) = -BAx*Rs(1:1,59,0)+ Rs(1:1,87,0)             
  CartRpy(1:1,59) = -BAy*Rs(1:1,59,0)+ Rs(1:1,88,0)             
  CartRpz(1:1,59) = -BAz*Rs(1:1,59,0)+ Rs(1:1,95,0)

  CartRpx(1:1,60) = -BAx*Rs(1:1,60,0)+ Rs(1:1,88,0)             
  CartRpy(1:1,60) = -BAy*Rs(1:1,60,0)+ Rs(1:1,89,0)             
  CartRpz(1:1,60) = -BAz*Rs(1:1,60,0)+ Rs(1:1,96,0)

  CartRpx(1:1,61) = -BAx*Rs(1:1,61,0)+ Rs(1:1,89,0)             
  CartRpy(1:1,61) = -BAy*Rs(1:1,61,0)+ Rs(1:1,90,0)             
  CartRpz(1:1,61) = -BAz*Rs(1:1,61,0)+ Rs(1:1,97,0)

  CartRpx(1:1,62) = -BAx*Rs(1:1,62,0)+ Rs(1:1,90,0)             
  CartRpy(1:1,62) = -BAy*Rs(1:1,62,0)+ Rs(1:1,91,0)             
  CartRpz(1:1,62) = -BAz*Rs(1:1,62,0)+ Rs(1:1,98,0)

  CartRpx(1:1,63) = -BAx*Rs(1:1,63,0)+ Rs(1:1,91,0)             
  CartRpy(1:1,63) = -BAy*Rs(1:1,63,0)+ Rs(1:1,92,0)             
  CartRpz(1:1,63) = -BAz*Rs(1:1,63,0)+ Rs(1:1,99,0)

  CartRpx(1:1,64) = -BAx*Rs(1:1,64,0)+ Rs(1:1,93,0)             
  CartRpy(1:1,64) = -BAy*Rs(1:1,64,0)+ Rs(1:1,94,0)             
  CartRpz(1:1,64) = -BAz*Rs(1:1,64,0)+ Rs(1:1,100,0)

  CartRpx(1:1,65) = -BAx*Rs(1:1,65,0)+ Rs(1:1,94,0)             
  CartRpy(1:1,65) = -BAy*Rs(1:1,65,0)+ Rs(1:1,95,0)             
  CartRpz(1:1,65) = -BAz*Rs(1:1,65,0)+ Rs(1:1,101,0)

  CartRpx(1:1,66) = -BAx*Rs(1:1,66,0)+ Rs(1:1,95,0)             
  CartRpy(1:1,66) = -BAy*Rs(1:1,66,0)+ Rs(1:1,96,0)             
  CartRpz(1:1,66) = -BAz*Rs(1:1,66,0)+ Rs(1:1,102,0)

  CartRpx(1:1,67) = -BAx*Rs(1:1,67,0)+ Rs(1:1,96,0)             
  CartRpy(1:1,67) = -BAy*Rs(1:1,67,0)+ Rs(1:1,97,0)             
  CartRpz(1:1,67) = -BAz*Rs(1:1,67,0)+ Rs(1:1,103,0)

  CartRpx(1:1,68) = -BAx*Rs(1:1,68,0)+ Rs(1:1,97,0)             
  CartRpy(1:1,68) = -BAy*Rs(1:1,68,0)+ Rs(1:1,98,0)             
  CartRpz(1:1,68) = -BAz*Rs(1:1,68,0)+ Rs(1:1,104,0)

  CartRpx(1:1,69) = -BAx*Rs(1:1,69,0)+ Rs(1:1,98,0)             
  CartRpy(1:1,69) = -BAy*Rs(1:1,69,0)+ Rs(1:1,99,0)             
  CartRpz(1:1,69) = -BAz*Rs(1:1,69,0)+ Rs(1:1,105,0)

  CartRpx(1:1,70) = -BAx*Rs(1:1,70,0)+ Rs(1:1,100,0)             
  CartRpy(1:1,70) = -BAy*Rs(1:1,70,0)+ Rs(1:1,101,0)             
  CartRpz(1:1,70) = -BAz*Rs(1:1,70,0)+ Rs(1:1,106,0)

  CartRpx(1:1,71) = -BAx*Rs(1:1,71,0)+ Rs(1:1,101,0)             
  CartRpy(1:1,71) = -BAy*Rs(1:1,71,0)+ Rs(1:1,102,0)             
  CartRpz(1:1,71) = -BAz*Rs(1:1,71,0)+ Rs(1:1,107,0)

  CartRpx(1:1,72) = -BAx*Rs(1:1,72,0)+ Rs(1:1,102,0)             
  CartRpy(1:1,72) = -BAy*Rs(1:1,72,0)+ Rs(1:1,103,0)             
  CartRpz(1:1,72) = -BAz*Rs(1:1,72,0)+ Rs(1:1,108,0)

  CartRpx(1:1,73) = -BAx*Rs(1:1,73,0)+ Rs(1:1,103,0)             
  CartRpy(1:1,73) = -BAy*Rs(1:1,73,0)+ Rs(1:1,104,0)             
  CartRpz(1:1,73) = -BAz*Rs(1:1,73,0)+ Rs(1:1,109,0)

  CartRpx(1:1,74) = -BAx*Rs(1:1,74,0)+ Rs(1:1,104,0)             
  CartRpy(1:1,74) = -BAy*Rs(1:1,74,0)+ Rs(1:1,105,0)             
  CartRpz(1:1,74) = -BAz*Rs(1:1,74,0)+ Rs(1:1,110,0)

  CartRpx(1:1,75) = -BAx*Rs(1:1,75,0)+ Rs(1:1,106,0)             
  CartRpy(1:1,75) = -BAy*Rs(1:1,75,0)+ Rs(1:1,107,0)             
  CartRpz(1:1,75) = -BAz*Rs(1:1,75,0)+ Rs(1:1,111,0)

  CartRpx(1:1,76) = -BAx*Rs(1:1,76,0)+ Rs(1:1,107,0)             
  CartRpy(1:1,76) = -BAy*Rs(1:1,76,0)+ Rs(1:1,108,0)             
  CartRpz(1:1,76) = -BAz*Rs(1:1,76,0)+ Rs(1:1,112,0)

  CartRpx(1:1,77) = -BAx*Rs(1:1,77,0)+ Rs(1:1,108,0)             
  CartRpy(1:1,77) = -BAy*Rs(1:1,77,0)+ Rs(1:1,109,0)             
  CartRpz(1:1,77) = -BAz*Rs(1:1,77,0)+ Rs(1:1,113,0)

  CartRpx(1:1,78) = -BAx*Rs(1:1,78,0)+ Rs(1:1,109,0)             
  CartRpy(1:1,78) = -BAy*Rs(1:1,78,0)+ Rs(1:1,110,0)             
  CartRpz(1:1,78) = -BAz*Rs(1:1,78,0)+ Rs(1:1,114,0)

  CartRpx(1:1,79) = -BAx*Rs(1:1,79,0)+ Rs(1:1,111,0)             
  CartRpy(1:1,79) = -BAy*Rs(1:1,79,0)+ Rs(1:1,112,0)             
  CartRpz(1:1,79) = -BAz*Rs(1:1,79,0)+ Rs(1:1,115,0)

  CartRpx(1:1,80) = -BAx*Rs(1:1,80,0)+ Rs(1:1,112,0)             
  CartRpy(1:1,80) = -BAy*Rs(1:1,80,0)+ Rs(1:1,113,0)             
  CartRpz(1:1,80) = -BAz*Rs(1:1,80,0)+ Rs(1:1,116,0)

  CartRpx(1:1,81) = -BAx*Rs(1:1,81,0)+ Rs(1:1,113,0)             
  CartRpy(1:1,81) = -BAy*Rs(1:1,81,0)+ Rs(1:1,114,0)             
  CartRpz(1:1,81) = -BAz*Rs(1:1,81,0)+ Rs(1:1,117,0)

  CartRpx(1:1,82) = -BAx*Rs(1:1,82,0)+ Rs(1:1,115,0)             
  CartRpy(1:1,82) = -BAy*Rs(1:1,82,0)+ Rs(1:1,116,0)             
  CartRpz(1:1,82) = -BAz*Rs(1:1,82,0)+ Rs(1:1,118,0)

  CartRpx(1:1,83) = -BAx*Rs(1:1,83,0)+ Rs(1:1,116,0)             
  CartRpy(1:1,83) = -BAy*Rs(1:1,83,0)+ Rs(1:1,117,0)             
  CartRpz(1:1,83) = -BAz*Rs(1:1,83,0)+ Rs(1:1,119,0)

  CartRpx(1:1,84) = -BAx*Rs(1:1,84,0)+ Rs(1:1,118,0)             
  CartRpy(1:1,84) = -BAy*Rs(1:1,84,0)+ Rs(1:1,119,0)             
  CartRpz(1:1,84) = -BAz*Rs(1:1,84,0)+ Rs(1:1,120,0)



return
end subroutine NEWCartRs2CartRp


!************************************************************************************

