subroutine calcSsmodels(shellmodelRGmat,lenshellmodel,shellmodelids, no, RGinr, maxmodellen,kmax, lenmdl,Ltrial,&
     concmodeloption,rrrrpreftable,maxrampdegree, Lmax, mdlbasis, index,mdlchk)
  implicit none
  real*8, intent(out)      :: shellmodelRGmat(1:1,1:3,1:maxmodellen), shellmodelids(1:1,1:4)
  integer, intent(out)     :: lenshellmodel(1:1)
  integer, intent(in)      :: no,  maxmodellen, kmax, concmodeloption, mdlchk
  integer,intent(in)       :: lenmdl, Ltrial, maxrampdegree, Lmax,mdlbasis(1:3,1:index), index
  real*8, intent(in)       :: RGinr(1:30,1:1), rrrrpreftable(0:maxrampdegree,0:maxrampdegree,0:Lmax)
  real*8                   :: Matnew(1:index,1:index),Vec(1:index,1:1), DummyArray(1:index), Rs(1:1), RGr(1:30,1:1)
  real*8                   :: Mat(1:index,1:index), condition,  t(1:10),  n(1:1), Vecnew(1:index,1:1), pi, x
  integer                  :: rank,nLCart, nLPure, cnt, st, fin, swit, timecnt, ai
  integer                  :: K,ntot, i, j, info, n1, K1, L1, n2, K2, L2, rn, m, L, kcurrent
  character(len=1)         :: outputchar1
  
  lenshellmodel = 0d0; shellmodelids = 0d0; shellmodelRGmat = 0d0;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This calculates the coefficients of the model for non-concentric Ss in terms of a sum of ramps of different angular momentum and degrees
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  pi = dacos(-1d0); n(1) = RGinr(7,1);    RGr = RGinr;

     !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     MATRIX CONSTRUCTION
! Calculate overlap of model basis functions with each other 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  Mat = 0d0;
  do i=1, index; 
     n1 = mdlbasis(1,i);   K1 = mdlbasis(2,i);   L1 = mdlbasis(3,i);
     do j=1,index; 
        n2 = mdlbasis(1,j);   K2 = mdlbasis(2,j);   L2 = mdlbasis(3,j);
        if (K1 == K2) then;
           Mat(i,j) = exp(dlgama(3d0+2d0*l1)+dlgama(1d0+n1+n2)-dlgama(4d0+2d0*l1+n1+n2));
        end if
     end do; 
  end do; 

!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     VECTOR CONSTRUCTION
! Calculate overlap of model basis functions with gaussian 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************

  nLCart = (Lmax+1)*(Lmax+2)*(Lmax+3)/6;     nLPure = (1+Lmax)**2; 
  do i = 1, index; 
     n1 = mdlbasis(1,i);   K1 = mdlbasis(2,i);   L1 = mdlbasis(3,i);     RGr(7,1:1) = n1;
     call calcpureRs(Rs,Rgr,no,K1,L1, kmax)
     if (Rgr(8,1) > 10) then; call calcpureRs(Rs(1),Rgr(1:30,1),1,K1,L1,10*kmax); end if
     Vec(i,1) = Rs(1);  
  end do

!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     LINEAR SOLVE
! Calculate coefficients by solving Mc = v 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
     st = 0; fin =0;  cnt = 0;


     do kcurrent = 1,1;
        st = fin; swit = 0;
        do k=st+1,index; if (swit == 0 .and. mdlbasis(2,k) >  kcurrent) then;  swit = 1; fin = k-1; end if;
        end do
        if (swit == 0) then;      fin = index;          end if
        if (fin-st .ne. 0) then;
           Matnew(st+1:fin,st+1:fin) = Mat(st+1:fin,st+1:fin); 
           do ai=st+1,fin;
              n1 = mdlbasis(1,ai)+n(1)
              Matnew(ai,fin+1) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
              Matnew(fin+1,ai) = 2d0/((n1+1d0)*(n1+2)*(n1+3))
        end do
           Vecnew = Vec;
           RGr(7,1:1) = n(1); 
           call calcpureRs(Rs,Rgr,1,1,0, kmax)
           Vecnew(fin+1,1) = Rs(1);
           call dgesv(fin-st+1,1, Matnew(st+1:fin+1,st+1:fin+1),fin-st+1, DummyArray, Vecnew(st+1:fin+1,1), fin-st+1, info)
    
      end if
     do k=st+1,fin ;   
           if (Abs(Vecnew(k,1)) > 1d-10) then
              cnt = cnt + 1;
              shellmodelRgmat(1,1,cnt) = dble(Vecnew(k,1));  
              shellmodelRGmat(1,2,cnt) = mdlbasis(1,k)+n(1); ! Degree of model ramp
              shellmodelRGmat(1,3,cnt) = mdlbasis(2,k);      ! Angular momentum of model ramp in K form
           end if;
        end do

     end do


     do kcurrent = 2,(1+Lmax)**2;
        st = fin; swit = 0;
        do k=st+1,index; if (swit == 0 .and. mdlbasis(2,k) >  kcurrent) then;  swit = 1; fin = k-1; end if;
        end do
        if (swit == 0) then;      fin = index;          end if
        if (fin-st .ne. 0) then;
           Matnew(st+1:fin,st+1:fin) = Mat(st+1:fin,st+1:fin); 
           call dgesv(fin-st,1, Matnew(st+1:fin,st+1:fin),fin-st, DummyArray, Vec(st+1:fin,1), fin-st, info)
        end if
        do k=st+1,fin ;   
           if (Abs(Vec(k,1)) > 1d-10) then
              cnt = cnt + 1;
              shellmodelRgmat(1,1,cnt) = dble(Vec(k,1));  
              shellmodelRGmat(1,2,cnt) = mdlbasis(1,k)+n(1); ! Degree of model ramp
              shellmodelRGmat(1,3,cnt) = mdlbasis(2,k);      ! Angular momentum of model ramp in K form
           end if;
        end do
     end do

  lenshellmodel(1) = cnt;
  shellmodelids(1,1:3) = RGr(1:3,1); shellmodelids(1,4) = RGr(9,1);

  return
end subroutine calcSsmodels


!********************************************************

subroutine calcpureRs(Rs,RgSs,no,K,L, kmax) 
  implicit none
  integer, intent(in)     :: no,  K, L, kmax
  real*8, intent(in)     :: RGSs(1:30,1:1)
  real*8, intent(out)    :: Rs(1:1)
  integer                :: nLCart, nLPure, i, n(1:1)
  real*8                 :: CartRs(1:1,1:(1+L)*(2+L)*(3+L)/6,0:L), PureDef1Rs(1:1,1:(1+L)**2), Fund(1:1,0:L)
  real*8                 :: beta(1:1), BAx(1:1), BAy(1:1), BAz(1:1), T(1:1), h(1:1,0:kmax,0:L)
  
  nLCart = (1+L)*(2+L)*(3+L)/6;     nLPure = (1+L)**2
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!  This subroutine calculates pure [R|s] overlap integral over all space
 !!  int( (1-r)^n r^{l} RY_{l,m}[theta,ph] e^{-\beta |r-P|^2} ) r^2 sin[theta] dr dtheta dphi
 !!
 !! where RY_{l,m} =         Y_{l,m}, m ==0
 !!                = 1/sqrt(2) * (Y_{l, m} + (-1)^m Y_{l,-m}) , m > 0
 !!                = 1/sqrt(2) * (Y_{l,-m} + (-1)^m Y_{l, m}) , m < 0
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 !! The methodology here is to first calculate integral over all space
 !!  int( (1-r)^n x^{i} y^{j} z^{k} e^{-\beta |r-P|^2} ) r^2 sin[theta] dr dtheta dphi
 !! 
 !! Then transform these to the desired integrals. 

 !! In all cases, the beta is small or P = 0 (i.e. not concentric and large beta), so we can use general method (truncating infinite series)

  beta =RGSs(8,1:1) ;   BAx = RGSs(4,1:1) ;   BAy = RGSs(5,1:1) ;   BAz = RGSs(6,1:1) ;
  T = beta * (BAx*BAx + BAy*BAy + BAz*BAz);  n = RgSs(7,1:1)
  
  call NEWcalcmphk(h,kmax,L,T,no)                                   !  Compute h_k(T) and their derivatives
  call NEWcalcintmpSs(Fund,h,kmax,L,n,beta,no)                      !  Compute [S|s]^(m) integrals
  call NEWcalcmpCartRs(CartRs,Fund,nLCart,L,beta,BAx,BAy,BAz,no)        !  Compute [R|s] integrals
  call NEWCart2PureDef1(PureDef1Rs,CartRs,no,nLCart,nLPure, L)
  Rs = PureDef1Rs(1:1,K)
  
  return
end subroutine calcpureRs



!************************************************************************************

!****************************************************************

subroutine NEWcalcmphk(h,kmax,L,T,no)
  implicit none
  integer,intent(in)            :: kmax, L, no
  real*8 ,intent(in)           :: T(1:1)
  real*8 ,intent(out)          :: h(1:1,0:kmax,0:L)
  integer                       :: k,m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the h_{i}^{m} functions that form the basis of [R|s] overlap integral evaluation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  h = 0d0;
  h(1:1,0,0) = exp(-T) * 2
  h(1:1,1,0) = (4*T-6) * h(1:1,0,0)
  do k = 2,kmax                                      !  Compute h_k(T) using forward recuCartRsion
     h(1:1,k,0) = (4*T-8*k+2) * h(1:1,k-1,0) - 8*(k-1)*(2*k-1) * h(1:1,k-2,0)
  end do
  do m = 1,L                                         !  Compute the mth derivatives of h_k(T)
     h(1:1,0,m) = - h(1:1,0,m-1)
     do k = 1,kmax
        h(1:1,k,m) = - h(1:1,k,m-1) - 4*k * h(1:1,k-1,m)
     end do
  end do
  return
end subroutine NEWcalcmphk

!************************************************************************************************

subroutine NEWcalcintmpSs(Fund,h,kmax,L,n,beta,no)
  implicit none
  integer,intent(in)            :: kmax,L,n(1:1), no
  real*8,intent(in)            :: h(1:1,0:kmax,0:L), beta(1:1)
  real*8,intent(out)           :: Fund(1:1,0:L)
  integer                       :: k,m
  real*8                       :: coeff(1:1)
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculates the [0]^m integrals for [S|s]^{m} 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Fund = 0d0
  coeff = 12.56637061435917d0 / (n+1d0) / beta
  do k = 0,kmax                                      !  Compute [S|s]^(m) from h_k^(m)(T)
     coeff = coeff * beta / ((n+2*k+2d0)*(n+2d0*k+3d0))
     do m = 0,L
        Fund(1:1,m) = Fund(1:1,m) + (k+1)*coeff * h(1:1,k,m)
     end do
  end do
  return
end subroutine NEWcalcintmpSs

!************************************************************************************************

subroutine NEWcalcmpCartRs(CartRs,Fund,nLCart,L,beta,BAx,BAy,BAz,no)
  implicit none
  integer,intent(in)  :: nLCart, L
  real*8,intent(in)  :: beta(1:1), BAx(1:1), BAy(1:1), BAz(1:1), Fund(1:1,0:L)
  real*8,intent(out) :: CartRs(1:1,1:nLCart,0:L)
  integer             :: m, i, no
  real*8             :: r2beta(1:1),t000(1:1),t100(1:1),t010(1:1)
  real*8             :: t001(1:1),t200(1:1),t110(1:1),t020(1:1),t002(1:1)
  real*8             :: t300(1:1),t120(1:1),t030(1:1),t201(1:1),t012(1:1),t003(1:1)
  real*8             :: t210(1:1),t102(1:1),t021(1:1),t400(1:1),t040(1:1),t004(1:1)
  real*8             :: t220(1:1),t022(1:1),t202(1:1),t031(1:1),t013(1:1),t310(1:1)
  real*8             :: t500(1:1),t050(1:1),t005(1:1),t320(1:1),t230(1:1),t410(1:1)
  real*8             :: t302(1:1),t130(1:1),t203(1:1),t122(1:1),t140(1:1),t103(1:1)
  real*8             :: t014(1:1),t023(1:1),t032(1:1)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This converts the fundamental [0]^m integrals to [R|s] overlap integrals via "VRR" type approach 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  CartRs = 0d0;
  
  r2beta = 1d0 / (2*beta)
  
  CartRs(1:1,1,0:L) = Fund(1:1,0:L)                                   !  [S|s]^(m)
  
  do m = L-1,0,-1                                         !  [P|s]^(m)  4L flops
     t000 =  CartRs(1:1, 1,m) + CartRs(1:1, 1,m+1)
     CartRs(1:1, 2,m) = BAx * t000                      !  100
     CartRs(1:1, 3,m) = BAy * t000                      !  010
     CartRs(1:1, 4,m) = BAz * t000                      !  001
  end do
  
  do m = L-2,0,-1                                         !  [D|s]^(m)  14(L-1) flops
     t000 = (CartRs(1:1, 1,m) + CartRs(1:1, 1,m+1)) * r2beta
     t100 =  CartRs(1:1, 2,m) + CartRs(1:1, 2,m+1)
     t010 =  CartRs(1:1, 3,m) + CartRs(1:1, 3,m+1)
     t001 =  CartRs(1:1, 4,m) + CartRs(1:1, 4,m+1)
     CartRs(1:1, 5,m) = BAx * t100 + t000               !  200
     CartRs(1:1, 6,m) = BAy * t100                      !  110
     CartRs(1:1, 7,m) = BAy * t010 + t000               !  020
     CartRs(1:1, 8,m) = BAx * t001                      !  101
     CartRs(1:1, 9,m) = BAy * t001                      !  011
     CartRs(1:1,10,m) = BAz * t001 + t000               !  002
  end do
  
  do m = L-3,0,-1                                         !  [F|s]^(m)  26(L-2) flops
     t100 = (CartRs(1:1, 2,m) + CartRs(1:1, 2,m+1)) * r2beta
     t010 = (CartRs(1:1, 3,m) + CartRs(1:1, 3,m+1)) * r2beta
     t001 = (CartRs(1:1, 4,m) + CartRs(1:1, 4,m+1)) * r2beta
     t200 =  CartRs(1:1, 5,m) + CartRs(1:1, 5,m+1)
     t110 =  CartRs(1:1, 6,m) + CartRs(1:1, 6,m+1)
     t020 =  CartRs(1:1, 7,m) + CartRs(1:1, 7,m+1)
     t002 =  CartRs(1:1,10,m) + CartRs(1:1,10,m+1)
     CartRs(1:1,11,m) = BAx * t200 + t100 * 2           !  300
     CartRs(1:1,12,m) = BAy * t200                      !  210
     CartRs(1:1,13,m) = BAx * t020                      !  120
     CartRs(1:1,14,m) = BAy * t020 + t010 * 2           !  030
     CartRs(1:1,15,m) = BAz * t200                      !  201
     CartRs(1:1,16,m) = BAz * t110                      !  111
     CartRs(1:1,17,m) = BAz * t020                      !  021
     CartRs(1:1,18,m) = BAx * t002                      !  102
     CartRs(1:1,19,m) = BAy * t002                      !  012
     CartRs(1:1,20,m) = BAz * t002 + t001 * 2           !  003
  end do
  
  do m = L-4,0,-1                                         !  [G|s]^(m)  36(L-3) flops
     t200 = (CartRs(1:1, 5,m) + CartRs(1:1, 5,m+1)) * r2beta
     t020 = (CartRs(1:1, 7,m) + CartRs(1:1, 7,m+1)) * r2beta
     t002 = (CartRs(1:1,10,m) + CartRs(1:1,10,m+1)) * r2beta
     t300 =  CartRs(1:1,11,m) + CartRs(1:1,11,m+1)
     t120 =  CartRs(1:1,13,m) + CartRs(1:1,13,m+1)
     t030 =  CartRs(1:1,14,m) + CartRs(1:1,14,m+1)
     t201 =  CartRs(1:1,15,m) + CartRs(1:1,15,m+1)
     t012 =  CartRs(1:1,19,m) + CartRs(1:1,19,m+1)
     t003 =  CartRs(1:1,20,m) + CartRs(1:1,20,m+1)
     CartRs(1:1,21,m) = BAx * t300 + t200 * 3           !  400
     CartRs(1:1,22,m) = BAy * t300                      !  310
     CartRs(1:1,23,m) = BAx * t120 + t020               !  220
     CartRs(1:1,24,m) = BAx * t030                      !  130
     CartRs(1:1,25,m) = BAy * t030 + t020 * 3           !  040
     CartRs(1:1,26,m) = BAz * t300                      !  301
     CartRs(1:1,27,m) = BAy * t201                      !  211
     CartRs(1:1,28,m) = BAz * t120                      !  121
     CartRs(1:1,29,m) = BAz * t030                      !  031
     CartRs(1:1,30,m) = BAz * t201 + t200               !  202
     CartRs(1:1,31,m) = BAx * t012                      !  112
     CartRs(1:1,32,m) = BAy * t012 + t002               !  022
     CartRs(1:1,33,m) = BAx * t003                      !  103
     CartRs(1:1,34,m) = BAy * t003                      !  013
     CartRs(1:1,35,m) = BAz * t003 + t002 * 3           !  004
   end do


   do m = L-5,0,-1                                         !  [H|s]^(m)  
      t300 = (CartRs(1:1,11,m) + CartRs(1:1,11,m+1)) * r2beta
      t030 = (CartRs(1:1,14,m) + CartRs(1:1,14,m+1)) * r2beta
      t003 = (CartRs(1:1,20,m) + CartRs(1:1,20,m+1)) * r2beta
      t120 = (CartRs(1:1,13,m) + CartRs(1:1,13,m+1)) * r2beta
      t210 = (CartRs(1:1,12,m) + CartRs(1:1,12,m+1)) * r2beta
      t102 = (CartRs(1:1,18,m) + CartRs(1:1,18,m+1)) * r2beta
      t201 = (CartRs(1:1,15,m) + CartRs(1:1,15,m+1)) * r2beta
      t012 = (CartRs(1:1,19,m) + CartRs(1:1,19,m+1)) * r2beta
      t021 = (CartRs(1:1,17,m) + CartRs(1:1,17,m+1)) * r2beta

      t400 =  CartRs(1:1,21,m) + CartRs(1:1,21,m+1)
      t040 =  CartRs(1:1,25,m) + CartRs(1:1,25,m+1)
      t004 =  CartRs(1:1,35,m) + CartRs(1:1,35,m+1)
      t220 =  CartRs(1:1,23,m) + CartRs(1:1,23,m+1)
      t202 =  CartRs(1:1,30,m) + CartRs(1:1,30,m+1)
      t022 =  CartRs(1:1,32,m) + CartRs(1:1,32,m+1)
      t013 =  CartRs(1:1,34,m) + CartRs(1:1,34,m+1)
      t031 =  CartRs(1:1,29,m) + CartRs(1:1,29,m+1)
      t310 =  CartRs(1:1,22,m) + CartRs(1:1,22,m+1)

      CartRs(1:1,36,m) = BAx * t400 + t300 * 4          !  500
      CartRs(1:1,37,m) = BAy * t400                     !  410
      CartRs(1:1,38,m) = BAx * t220 + t120 * 2          !  320 
      CartRs(1:1,39,m) = BAy * t220 + t210 * 2          !  230
      CartRs(1:1,40,m) = BAx * t040                     !  140
      CartRs(1:1,41,m) = BAy * t040 + t030 * 4          !  050
      CartRs(1:1,42,m) = BAz * t400                     !  401
      CartRs(1:1,43,m) = BAz * t310                     !  311
      CartRs(1:1,44,m) = BAz * t220                     !  221
      CartRs(1:1,45,m) = BAx * t031                     !  131
      CartRs(1:1,46,m) = BAz * t040                     !  041
      CartRs(1:1,47,m) = BAx * t202 + t102 * 2          !  302
      CartRs(1:1,48,m) = BAy * t202                     !  212
      CartRs(1:1,49,m) = BAx * t022                     !  122
      CartRs(1:1,50,m) = BAy * t022 + t012 * 2          !  032
      CartRs(1:1,51,m) = BAz * t202 + t201 * 2          !  203
      CartRs(1:1,52,m) = BAx * t013                     !  113
      CartRs(1:1,53,m) = BAz * t022 + t021 * 2          !  023
      CartRs(1:1,54,m) = BAx * t004                     !  104
      CartRs(1:1,55,m) = BAy * t004                     !  014
      CartRs(1:1,56,m) = BAz * t004 + t003 * 4          !  005



   end do 

   do m = L-6,0,-1

      t400 = (CartRs(1:1,21,m) + CartRs(1:1,21,m+1)) * r2beta
      t040 = (CartRs(1:1,25,m) + CartRs(1:1,25,m+1)) * r2beta
      t004 = (CartRs(1:1,35,m) + CartRs(1:1,35,m+1)) * r2beta
      t130 = (CartRs(1:1,24,m) + CartRs(1:1,24,m+1)) * r2beta
      t103 = (CartRs(1:1,33,m) + CartRs(1:1,33,m+1)) * r2beta
      t013 = (CartRs(1:1,34,m) + CartRs(1:1,34,m+1)) * r2beta
      t220 = (CartRs(1:1,23,m) + CartRs(1:1,23,m+1)) * r2beta
      t202 = (CartRs(1:1,30,m) + CartRs(1:1,30,m+1)) * r2beta
      t022 = (CartRs(1:1,32,m) + CartRs(1:1,32,m+1)) * r2beta
 
      t500 =  CartRs(1:1,36,m) + CartRs(1:1,36,m+1)
      t050 =  CartRs(1:1,41,m) + CartRs(1:1,41,m+1)
      t005 =  CartRs(1:1,56,m) + CartRs(1:1,56,m+1)
      t410 =  CartRs(1:1,37,m) + CartRs(1:1,37,m+1)
      t140 =  CartRs(1:1,40,m) + CartRs(1:1,40,m+1)
      t230 =  CartRs(1:1,39,m) + CartRs(1:1,39,m+1)
      t320 =  CartRs(1:1,38,m) + CartRs(1:1,38,m+1)
      t302 =  CartRs(1:1,47,m) + CartRs(1:1,47,m+1)
      t014 =  CartRs(1:1,55,m) + CartRs(1:1,55,m+1)
      t023 =  CartRs(1:1,53,m) + CartRs(1:1,53,m+1)
      t203 =  CartRs(1:1,51,m) + CartRs(1:1,51,m+1)
      t122 =  CartRs(1:1,49,m) + CartRs(1:1,49,m+1)
      t032 =  CartRs(1:1,50,m) + CartRs(1:1,50,m+1)


      CartRs(1:1,57,m) = BAz * t500 + t400 * 5       ! 600
      CartRs(1:1,58,m) = BAy * t500                  ! 510
      CartRs(1:1,59,m) = BAx * t320 + t220 * 3       ! 420
      CartRs(1:1,60,m) = BAx * t230 + t130 * 2       ! 330
      CartRs(1:1,61,m) = BAy * t230 + t220 * 3       ! 240
      CartRs(1:1,62,m) = BAx * t050                  ! 150
      CartRs(1:1,63,m) = BAy * t050 + t040 * 5       ! 060
      CartRs(1:1,64,m) = BAz * t500                  ! 501
      CartRs(1:1,65,m) = BAz * t410                  ! 411
      CartRs(1:1,66,m) = BAz * t320                  ! 321
      CartRs(1:1,67,m) = BAz * t230                  ! 231
      CartRs(1:1,68,m) = BAz * t140                  ! 141
      CartRs(1:1,69,m) = BAz * t050                  ! 051
      CartRs(1:1,70,m) = BAx * t302 + t202 * 3       ! 402
      CartRs(1:1,71,m) = BAy * t302                  ! 312
      CartRs(1:1,72,m) = BAx * t122 + t022           ! 222
      CartRs(1:1,73,m) = BAx * t032                  ! 132
      CartRs(1:1,74,m) = BAy * t032 + t022 * 3       ! 042
      CartRs(1:1,75,m) = BAx * t203 + t103 * 2       ! 303
      CartRs(1:1,76,m) = BAy * t203                  ! 213
      CartRs(1:1,77,m) = BAx * t023                  ! 123
      CartRs(1:1,78,m) = BAy * t023 + t013 * 2       ! 033
      CartRs(1:1,79,m) = BAz * t203 + t202 * 3       ! 204
      CartRs(1:1,80,m) = BAx * t014                  ! 114
      CartRs(1:1,81,m) = BAz * t023 + t022 * 3       ! 024
      CartRs(1:1,82,m) = BAx * t005                  ! 105
      CartRs(1:1,83,m) = BAy * t005                  ! 015
      CartRs(1:1,84,m) = BAz * t005 + t004 * 5       ! 006
   end do
 


      return
end subroutine NEWcalcmpCartRs

!**********************************************************************


!****************************************************
subroutine NEWCart2PureDef1(PureDef1Rs, CartRs,no,nLCart,nLPure, L)
   implicit none
   integer, intent(in)    :: L, no, nLCart, nLPure
   real*8, intent(in)    :: CartRs(1:1,1:nLCart,0:L)
   real*8, intent(out)   :: PureDef1Rs(1:1,1:nLPure)
   real*8                :: t(1:1), invfact(0:2*L)
   integer                :: lindex, mindex, cnt, i, x, fact(0:2*L)
   real*8                :: CY500(1:1), CY410(1:1), CY320(1:1), CY230(1:1), CY140(1:1)
   real*8                :: CY050(1:1), CY401(1:1), CY311(1:1), CY221(1:1), CY131(1:1)
   real*8                :: CY041(1:1), CY302(1:1), CY212(1:1), CY122(1:1), CY032(1:1)
   real*8                :: CY203(1:1), CY113(1:1), CY023(1:1), CY104(1:1), CY014(1:1), CY005(1:1)
   real*8 :: CY220(1:1), CY202(1:1), CY022(1:1), CY400(1:1), CY040(1:1), CY004(1:1)
   real*8 :: CY600(1:1),CY510(1:1),CY420(1:1),CY330(1:1),CY240(1:1),CY150(1:1),CY060(1:1)
   real*8 :: CY501(1:1),CY411(1:1),CY321(1:1),CY231(1:1),CY141(1:1),CY051(1:1)
   real*8 :: CY402(1:1),CY312(1:1),CY222(1:1),CY132(1:1),CY042(1:1)
   real*8 :: CY303(1:1),CY213(1:1),CY123(1:1),CY033(1:1)
   real*8 :: CY204(1:1),CY114(1:1),CY024(1:1)
   real*8 :: CY105(1:1),CY015(1:1),CY006(1:1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This converts Cartesian [R|s] overlap integrals to Pure [R|s] overlap integrals.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   PureDef1Rs = 0d0;

   PureDef1Rs(1:1,1) =  0.28209479177387814347d0*CartRs(1:1,1,0)

if (nLpure == 1) then; return; end if
   PureDef1Rs(1:1,2) = -0.48860251190291992159d0*CartRs(1:1,3,0)
   PureDef1Rs(1:1,3) =  0.48860251190291992159d0*CartRs(1:1,4,0)
   PureDef1Rs(1:1,4) = -0.48860251190291992159d0*CartRs(1:1,2,0)

if (nLpure == 4) then; return; end if
   PureDef1Rs(1:1,5) = 1.0925484305920790705d0*CartRs(1:1,6,0) !l=2, m=-2
   PureDef1Rs(1:1,6) = -1.0925484305920790705d0*CartRs(1:1,9,0) !l=2, m=-1
   PureDef1Rs(1:1,7) = 0.31539156525252000603d0*(CartRs(1:1,10,0)*2-CartRs(1:1,5,0)-CartRs(1:1,7,0))   !l=2, m=0
   PureDef1Rs(1:1,8) = -1.0925484305920790705d0*CartRs(1:1,8,0) !l=2, m=1
   PureDef1Rs(1:1,9) = 0.54627421529603953527d0*(CartRs(1:1,5,0) - CartRs(1:1,7,0))  !l=2, m=2

if (nLpure == 9) then; return; end if

   PureDef1Rs(1:1,10) = -0.59004358992664351035d0*(3*CartRs(1:1,12,0) - CartRs(1:1,14,0)) !l=3, m=-3
   PureDef1Rs(1:1,11) = 2.8906114426405540554d0*CartRs(1:1,16,0) !l=3, m=-2
   PureDef1Rs(1:1,12) = 0.45704579946446573616d0*(CartRs(1:1,12,0)+CartRs(1:1,14,0) - 4d0*CartRs(1:1,19,0)) !l=3, m=-1
   PureDef1Rs(1:1,13) = -1.1195289977703461742d0*(CartRs(1:1,15,0)+CartRs(1:1,17,0)) + 0.74635266518023078283d0*CartRs(1:1,20,0) !l=3,m=0
   PureDef1Rs(1:1,14) = 0.45704579946446573616d0*(CartRs(1:1,11,0)+CartRs(1:1,13,0)-4d0*CartRs(1:1,18,0))  !l=3, m=1
   PureDef1Rs(1:1,15) = 1.4453057213202770277d0*(CartRs(1:1,15,0)-CartRs(1:1,17,0)) !l=3, m=2
   PureDef1Rs(1:1,16) =  -0.59004358992664351035d0*(CartRs(1:1,11,0) - 3d0*CartRs(1:1,13,0)) !l=3, m=3


if (nLpure == 16) then; return; end if

   PureDef1Rs(1:1,17) = 2.5033429417967045383d0*(CartRs(1:1,22,0) - CartRs(1:1,24,0)) !l=4, m=-4
   PureDef1Rs(1:1,18) = -1.7701307697799305310d0*(3*CartRs(1:1,27,0) - CartRs(1:1,29,0)) !l=4, m=-3
   PureDef1Rs(1:1,19) = -0.94617469575756001809d0*(CartRs(1:1,22,0)+CartRs(1:1,24,0)-6d0*CartRs(1:1,31,0))  !l=4, m=-2
   PureDef1Rs(1:1,20) = 0.66904654355728916795d0*(3*CartRs(1:1,27,0)+3d0*CartRs(1:1,29,0) - 4*CartRs(1:1,34,0)) !l=4,m=-1


 CY220 = CartRs(1:1,23,0)
 CY202 = CartRs(1:1,30,0)
 CY022 = CartRs(1:1,32,0)
 CY400 = CartRs(1:1,21,0)
 CY040 = CartRs(1:1,25,0)
 CY004 = CartRs(1:1,35,0)

   PureDef1Rs(1:1,21) = 0.10578554691520430380d0*(8*CY004 - 24*CY022 + 3*CY040 -24*CY202 + 6*CY220 + 3*CY400) !l=4, m=0
   PureDef1Rs(1:1,22) = 0.66904654355728916795d0*(3*CartRs(1:1,26,0)+3d0*CartRs(1:1,28,0) -4*CartRs(1:1,33,0)) !l=4, m=1
   PureDef1Rs(1:1,23) = -0.55901699437494742410d0*(CartRs(1:1,21,0) -CartRs(1:1,25,0)&
          - 6*CartRs(1:1,30,0) + 6*CartRs(1:1,32,0)) !l=4, m=2  
   PureDef1Rs(1:1,24) = -1.7701307697799305310d0*(CartRs(1:1,26,0)-3d0*CartRs(1:1,28,0))!l=4, m=3
   PureDef1Rs(1:1,25) = 0.62583573544917613459d0*(CartRs(1:1,21,0)-6d0*CartRs(1:1,23,0)+CartRs(1:1,25,0)) !l=4, m=4 


if (nLpure == 25) then; return; end if

CY500 = CartRs(1:1,36,0)
CY410 = CartRs(1:1,37,0)
CY320 = CartRs(1:1,38,0)
CY230 = CartRs(1:1,39,0)
CY140 = CartRs(1:1,40,0)
CY050 = CartRs(1:1,41,0)
CY401 = CartRs(1:1,42,0)
CY311 = CartRs(1:1,43,0)
CY221 = CartRs(1:1,44,0)
CY131 = CartRs(1:1,45,0)
CY041 = CartRs(1:1,46,0)
CY302 = CartRs(1:1,47,0)
CY212 = CartRs(1:1,48,0)
CY122 = CartRs(1:1,49,0)
CY032 = CartRs(1:1,50,0)
CY203 = CartRs(1:1,51,0)
CY113 = CartRs(1:1,52,0)
CY023 = CartRs(1:1,53,0)
CY104 = CartRs(1:1,54,0)
CY014 = CartRs(1:1,55,0)
CY005 = CartRs(1:1,56,0)

   PureDef1Rs(1:1,26) = -0.65638205684017010281d0*(5*CY410 - 10*CY230 +   CY050 )           
   PureDef1Rs(1:1,27) =  8.3026492595241651160d0 *(  CY311 -    CY131           )             
   PureDef1Rs(1:1,28) =  0.48923829943525038768d0*(8*CY032 -    CY050 - 24*CY212 + 2*CY230 + 3*CY410)            
   PureDef1Rs(1:1,29) =   4.7935367849733237550d0*(2*CY113 -    CY131 -    CY311)            
   PureDef1Rs(1:1,30) = -0.45294665119569692130d0*(8*CY014 - 12*CY032 +    CY050 -12*CY212 + 2*CY230 +   CY410)            
   PureDef1Rs(1:1,31) =  0.11695032245342359644d0*(8*CY005 - 40*CY023 + 15*CY041 -40*CY203 +30*CY221 +15*CY401)            
   PureDef1Rs(1:1,32) = -0.45294665119569692130d0*(8*CY104 - 12*CY122 +    CY140 -12*CY302 + 2*CY320 +   CY500)              
   PureDef1Rs(1:1,33) =  -2.3967683924866618775d0*(2*CY023 -    CY041 -  2*CY203 +   CY401)           
   PureDef1Rs(1:1,34) =  0.48923829943525038768d0*(  CY500 + 24*CY122 -  3*CY140 - 8*CY302 - 2*CY320)           
   PureDef1Rs(1:1,35) =   2.0756623148810412790d0*(  CY401 -  6*CY221 +   CY041  )              
   PureDef1Rs(1:1,36) = -0.65638205684017010281d0*(5*CY140 - 10*CY320 +   CY500  )  

if (nLPure == 36) then; return; end if

 CY600 = CartRs(1:1,57,0)
 CY510 = CartRs(1:1,58,0)
 CY420 = CartRs(1:1,59,0)
 CY330 = CartRs(1:1,60,0)
 CY240 = CartRs(1:1,61,0)
 CY150 = CartRs(1:1,62,0)
 CY060 = CartRs(1:1,63,0)
 CY501 = CartRs(1:1,64,0)
 CY411 = CartRs(1:1,65,0)
 CY321 = CartRs(1:1,66,0)
 CY231 = CartRs(1:1,67,0)
 CY141 = CartRs(1:1,68,0)
 CY051 = CartRs(1:1,69,0)
 CY402 = CartRs(1:1,70,0)
 CY312 = CartRs(1:1,71,0)
 CY222 = CartRs(1:1,72,0)
 CY132 = CartRs(1:1,73,0)
 CY042 = CartRs(1:1,74,0)
 CY303 = CartRs(1:1,75,0)
 CY213 = CartRs(1:1,76,0)
 CY123 = CartRs(1:1,77,0)
 CY033 = CartRs(1:1,78,0)
 CY204 = CartRs(1:1,79,0)
 CY114 = CartRs(1:1,80,0)
 CY024 = CartRs(1:1,81,0)
 CY105 = CartRs(1:1,82,0)
 CY015 = CartRs(1:1,83,0)
 CY006 = CartRs(1:1,84,0)           

   PureDef1Rs(1:1,37) =    1.3663682103838286439d0*( 3*CY150 - 10*CY330 +  3*CY510    )     
   PureDef1Rs(1:1,38) =   -2.3666191622317520320d0*(   CY051 - 10*CY231 +  5*CY411     )           
   PureDef1Rs(1:1,39) =   -2.0182596029148966370d0*(10*CY132 -    CY150 - 10*CY312 +    CY510     )            

   PureDef1Rs(1:1,40) =  -0.92120525951492349914d0*(24*CY213 -  8*CY033 +  3*CY051 -  6*CY231 -  9*CY411)           

   PureDef1Rs(1:1,41) =  -0.92120525951492349914d0*(16*CY132 -    CY150 + 16*CY312 -  2*CY330 -    CY510 - 16*CY114)            
   PureDef1Rs(1:1,42) =  -0.58262136251873138884d0*( 8*CY015 - 20*CY033 +  5*CY051 - 20*CY213 + 10*CY231 +  5*CY411)    
           
   PureDef1Rs(1:1,43) =  0.063569202267628425933d0*(16*CY006 -120*CY024 + 90*CY042 -  5*CY060 -120*CY204 +180*CY222 &
                               -15* CY240 + 90*CY402 -15*CY420 - 5*CY600)            
   PureDef1Rs(1:1,44) =  -0.58262136251873138884d0*( 8*CY105 - 20*CY123 +  5*CY141 - 20*CY303 + 10*CY321 +  5*CY501)            
   PureDef1Rs(1:1,45) =  -0.46060262975746174957d0*(16*CY024 - 16*CY042 +    CY060 - 16*CY204 +    CY240 + 16*CY402 &
                               - CY420 - CY600)             
   PureDef1Rs(1:1,46) =   0.92120525951492349914d0*(24*CY123 -  9*CY141 -  8*CY303 -  6*CY321 +  3*CY501)           
   PureDef1Rs(1:1,47) =   0.50456490072872415925d0*(10*CY042 -    CY060 - 60*CY222 +  5*CY240 + 10*CY402 &
                               + 5*CY420 - CY600)           
   PureDef1Rs(1:1,48) =   -2.3666191622317520320d0*(5*CY141 - 10*CY321 +    CY501        )             
   PureDef1Rs(1:1,49) =  -0.68318410519191432197d0*(  CY060 - 15*CY240 + 15*CY420 -  CY600)                 


return
end subroutine NEWCart2PureDef1

!*********************************************************************************
