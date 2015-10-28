!************************************************************************************

subroutine genrgSpke(kergSp,no,RGSp,kmax)
  implicit none
  integer,parameter                :: L = 3, nL = (L+1)*(L+2)*(L+3)/6
  integer, intent(in)              :: no, kmax
  real*8, intent(in)               :: RGSp(1:30,1:no)
  real*8, intent(out)              :: kergSp(1:3*no,1:4)
  integer                          :: n(1:no), i, lindex,mindex
  real*8                           :: keSP(1:3*no), h(1:no,0:kmax,0:L)
  real*8                           :: beta(1:no),  BAx(1:no), BAy(1:no), BAz(1:no), T(1:no)
  real*8                           :: Rs(1:no,1:nL,0:L), Fund(1:no,0:L)
  real*8                           :: Sp(1:no,1,2:4), Pp(1:no,2:4,2:4), Pd(1:no,2:4,5:10)
  real*8                           :: Dp(1:no,5:10,2:4), Sd(1:no,1,5:10), Sf(1:no,1,11:20)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program calculates (S|s) using the general Cold-like algorithm (truncation of expression)

! kmax = number of terms before truncation
! Sp = (S|T|p)^(0)
! Pp = (P|T|p)^(0)
! Pd = (P|T|d)^(0)
! Dp = (D|T|p)^(0)
! Sd = (S|T|d)^(0)
! Sf = (S|T|f)^(0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  kergSP = 0d0;
  do i=1,no
     kergSP(3*i-2:3*i,2) = RGSp(1,i)
     kergSP(3*i-2    ,3) = RGSp(2,i)
     kergSP(3*i-1    ,3) = RGSp(2,i)+1
     kergSP(3*i      ,3) = RGSp(2,i)+2
     kergSP(3*i-2:3*i,4) = RGSp(9,i)
  end do
  
  beta = RGSp(8,1:no) 
  BAx = RGSp(4,1:no) 
  BAy = RGSp(5,1:no) 
  BAz = RGSp(6,1:no)
  T = beta * (BAx*BAx + BAy*BAy + BAz*BAz)
  n = int(RGSp(7,1:no))
  
  call calckehkforSp(h,kmax,L,T,no)      
  call calcFundkeSp(Fund,h,kmax,L,n,beta,no)   
  call calcRskeSp(Rs,Fund,nL,L,beta,BAx,BAy,BAz,no)
  
  call calcSp0(Sp,Rs,no,nL,L,BAx,BAy,BAz) 
  call calcPp0(Pp,Rs,no,nL,L,BAx,BAy,BAz) 
  call calcDp0(Dp,Rs,no,nL,L,BAx,BAy,BAz)
  call calcPd0(Pd,Dp,Pp,no,nL,L,BAx,BAy,BAz)
  call calcSd0(Sd,Sp,Pp,no,nL,L,BAx,BAy,BAz)
  call calcSf0(Sf,Sd,Pd,no,nL,L,BAx,BAy,BAz)
  call calckeSp(keSp,Sp,Sf,no,beta)
  
  kergSP(1:3*no,1) = keSp
  
  return
end subroutine genrgSpke

!************************************************************************************************

subroutine calckehkforSp(h,kmax,L,T,no)
  implicit none
  integer,intent(in)            :: kmax, L, no
  real*8 ,intent(in)            :: T(1:no)
  real*8 ,intent(out)           :: h(1:no,0:kmax,0:L)
  integer                       :: k,m
  
  h = 0d0;
  h(1:no,0,0) = exp(-T) * 2
  h(1:no,1,0) = (4*T-6) * h(1:no,0,0)
  do k = 2,kmax                                      !  Cokeute h_k(T) using forward recursion
     h(1:no,k,0) = (4*T-8*k+2) * h(1:no,k-1,0) - 8*(k-1)*(2*k-1) * h(1:no,k-2,0)
  end do
  do m = 1,L                                         !  Cokeute the mth derivatives of h_k(T)
     h(1:no,0,m) = - h(1:no,0,m-1)
     do k = 1,kmax
        h(1:no,k,m) = - h(1:no,k,m-1) - 4*k * h(1:no,k-1,m)
     end do
  end do
  return
end subroutine calckehkforSp

!************************************************************************************************


subroutine calcFundkeSp(Fund,h,kmax,L,n,beta,no)
  implicit none
  integer,intent(in)            :: kmax,L,n(1:no), no
  real*8 ,intent(in)            :: h(1:no,0:kmax,0:L), beta(1:no)
  real*8 ,intent(out)           :: Fund(1:no,0:L)
  integer                       :: k,m
  real*8                        :: coeff(1:no)
  
  Fund = 0d0
  coeff =  4 * dacos(-1d0) / (n+1) / beta
  do k = 0,kmax                                      !  Calculate [S|s]^(m) from h_k^(m)(T)
     coeff = coeff * beta / ((n+2*k+2)*(n+2*k+3))
     do m = 0,L
        Fund(1:no,m) = Fund(1:no,m) + (k+1)*coeff * h(1:no,k,m)
     end do
  end do
  
  return
end subroutine calcFundkeSp


!**********************************************************

subroutine calcRskeSp(Rs,Fund,nL,L,beta,BAx,BAy,BAz,no)
  implicit none
  integer,intent(in)  :: nL, L, no
  real*8 ,intent(in)  :: beta(1:no), BAx(1:no), BAy(1:no), BAz(1:no),  Fund(1:no,0:L)
  real*8 ,intent(out) :: Rs(1:no,1:nL,0:L)
  integer             :: m, i
  real*8              :: r2beta(1:no),t000(1:no),t100(1:no),t010(1:no)
  real*8              :: t001(1:no),t200(1:no),t110(1:no),t020(1:no),t002(1:no)
  real*8              :: t300(1:no),t120(1:no),t030(1:no),t201(1:no),t012(1:no),t003(1:no)
  
  Rs = 0d0;
  r2beta = 1d0 / (2*beta)
  Rs(1:no,1,0:L) = Fund(1:no,0:L)                                   !  [S|s]^(m)
  
  do m = L-1,0,-1                                         !  [P|s]^(m)  4L flops
     t000 =  Rs(1:no, 1,m) + Rs(1:no, 1,m+1)
     Rs(1:no, 2,m) = BAx * t000                      !  100
     Rs(1:no, 3,m) = BAy * t000                      !  010
     Rs(1:no, 4,m) = BAz * t000                      !  001
  end do
  
  do m = L-2,0,-1                                         !  [D|s]^(m)  14(L-1) flops
     t000 = (Rs(1:no, 1,m) + Rs(1:no, 1,m+1)) * r2beta
     t100 =  Rs(1:no, 2,m) + Rs(1:no, 2,m+1)
     t010 =  Rs(1:no, 3,m) + Rs(1:no, 3,m+1)
     t001 =  Rs(1:no, 4,m) + Rs(1:no, 4,m+1)
     Rs(1:no, 5,m) = BAx * t100 + t000               !  200
     Rs(1:no, 6,m) = BAy * t100                      !  110
     Rs(1:no, 7,m) = BAy * t010 + t000               !  020
     Rs(1:no, 8,m) = BAx * t001                      !  101
     Rs(1:no, 9,m) = BAy * t001                      !  011
     Rs(1:no,10,m) = BAz * t001 + t000               !  002
  end do
  
  do m = L-3,0,-1                                         !  [F|s]^(m)  26(L-2) flops
     t100 = (Rs(1:no, 2,m) + Rs(1:no, 2,m+1)) * r2beta
     t010 = (Rs(1:no, 3,m) + Rs(1:no, 3,m+1)) * r2beta
     t001 = (Rs(1:no, 4,m) + Rs(1:no, 4,m+1)) * r2beta
     t200 =  Rs(1:no, 5,m) + Rs(1:no, 5,m+1)
     t110 =  Rs(1:no, 6,m) + Rs(1:no, 6,m+1)
     t020 =  Rs(1:no, 7,m) + Rs(1:no, 7,m+1)
     t002 =  Rs(1:no,10,m) + Rs(1:no,10,m+1)
     Rs(1:no,11,m) = BAx * t200 + t100 * 2           !  300
     Rs(1:no,12,m) = BAy * t200                      !  210
     Rs(1:no,13,m) = BAx * t020                      !  120
     Rs(1:no,14,m) = BAy * t020 + t010 * 2           !  030
     Rs(1:no,15,m) = BAz * t200                      !  201
     Rs(1:no,16,m) = BAz * t110                      !  111
     Rs(1:no,17,m) = BAz * t020                      !  021
     Rs(1:no,18,m) = BAx * t002                      !  102
     Rs(1:no,19,m) = BAy * t002                      !  012
     Rs(1:no,20,m) = BAz * t002 + t001 * 2           !  003
  end do
  
end subroutine calcRskeSp

!***************************************************
subroutine calcSp0(Sp,Rs,no,nL,L,BAx,BAy,BAz)
  implicit none
  integer, intent(in)   :: no, nL, L
  real*8, intent(out)   :: Sp(1:no,1,2:4)
  real*8, intent(in)    :: Rs(1:no,1:nL,0:L), BAx(1:no), BAy(1:no), BAz(1:no)
  
  Sp(1:no,1,2) = -BAx*Rs(1:no,1,0) + Rs(1:no,2,0)
  Sp(1:no,1,3) = -BAy*Rs(1:no,1,0) + Rs(1:no,3,0)
  Sp(1:no,1,4) = -BAz*Rs(1:no,1,0) + Rs(1:no,4,0)
  
  return
end subroutine calcSp0

!***************************************************

subroutine calcDp0(Dp,Rs,no,nL,L,BAx,BAy,BAz)
  implicit none
  integer, intent(in)   :: no, nL, L
  real*8, intent(out)   :: Dp(1:no,5:10,2:4)
  real*8, intent(in)    :: Rs(1:no,1:nL,0:L),  BAx(1:no), BAy(1:no), BAz(1:no)
  
  Dp(1:no, 5,2) = -BAx*Rs(1:no, 5,0) + Rs(1:no,11,0)
  Dp(1:no, 6,2) = -BAx*Rs(1:no, 6,0) + Rs(1:no,12,0)
  Dp(1:no, 7,2) = -BAx*Rs(1:no, 7,0) + Rs(1:no,13,0)
  Dp(1:no, 8,2) = -BAx*Rs(1:no, 8,0) + Rs(1:no,15,0)
  Dp(1:no, 9,2) = -BAx*Rs(1:no, 9,0) + Rs(1:no,16,0)
  Dp(1:no,10,2) = -BAx*Rs(1:no,10,0) + Rs(1:no,18,0)
  
  Dp(1:no, 5,3) = -BAy*Rs(1:no, 5,0) + Rs(1:no,12,0)
  Dp(1:no, 6,3) = -BAy*Rs(1:no, 6,0) + Rs(1:no,13,0)
  Dp(1:no, 7,3) = -BAy*Rs(1:no, 7,0) + Rs(1:no,14,0)
  Dp(1:no, 8,3) = -BAy*Rs(1:no, 8,0) + Rs(1:no,16,0)
  Dp(1:no, 9,3) = -BAy*Rs(1:no, 9,0) + Rs(1:no,17,0)
  Dp(1:no,10,3) = -BAy*Rs(1:no,10,0) + Rs(1:no,19,0)
  
  Dp(1:no, 5,4) = -BAz*Rs(1:no, 5,0) + Rs(1:no,15,0)
  Dp(1:no, 6,4) = -BAz*Rs(1:no, 6,0) + Rs(1:no,16,0)
  Dp(1:no, 7,4) = -BAz*Rs(1:no, 7,0) + Rs(1:no,17,0)
  Dp(1:no, 8,4) = -BAz*Rs(1:no, 8,0) + Rs(1:no,18,0)
  Dp(1:no, 9,4) = -BAz*Rs(1:no, 9,0) + Rs(1:no,19,0)
  Dp(1:no,10,4) = -BAz*Rs(1:no,10,0) + Rs(1:no,20,0)
  
  return
end subroutine calcDp0

!***************************************************
subroutine calcPp0(Pp,Rs,no,nL,L,BAx,BAy,BAz)
  implicit none
  integer, intent(in)   :: no, nL, L
  real*8, intent(out)   :: Pp(1:no,2:4,2:4)
  real*8, intent(in)    :: Rs(1:no,1:nL,0:L),  BAx(1:no), BAy(1:no), BAz(1:no)
  
  Pp(1:no,2,2) = -BAx*Rs(1:no,2,0) + Rs(1:no,5,0)
  Pp(1:no,3,2) = -BAy*Rs(1:no,2,0) + Rs(1:no,6,0)
  Pp(1:no,4,2) = -BAz*Rs(1:no,2,0) + Rs(1:no,8,0)
  
  Pp(1:no,2,3) = -BAx*Rs(1:no,3,0) + Rs(1:no,6,0)
  Pp(1:no,3,3) = -BAy*Rs(1:no,3,0) + Rs(1:no,7,0)
  Pp(1:no,4,3) = -BAz*Rs(1:no,3,0) + Rs(1:no,9,0)
  
  Pp(1:no,2,4) = -BAx*Rs(1:no,4,0) + Rs(1:no,8,0)
  Pp(1:no,3,4) = -BAy*Rs(1:no,4,0) + Rs(1:no,9,0)
  Pp(1:no,4,4) = -BAz*Rs(1:no,4,0) + Rs(1:no,10,0)
  
  return
end subroutine calcPp0

!**********************************************************


subroutine calcPd0(Pd,Dp,Pp,no,nL,L,BAx,BAy,BAz)
  implicit none
  integer, intent(in)   :: no, nL, L
  real*8, intent(out)   :: Pd(1:no,2:4,5:10)
  real*8, intent(in)    :: Pp(1:no,2:4,2:4), Dp(1:no,5:10,2:4), BAx(1:no), BAy(1:no), BAz(1:no)
  
  Pd(1:no, 2,5) = -BAx*Pp(1:no,2,2) + Dp(1:no,5,2)
  Pd(1:no, 3,5) = -BAx*Pp(1:no,3,2) + Dp(1:no,6,2)
  Pd(1:no, 4,5) = -BAx*Pp(1:no,4,2) + Dp(1:no,8,2)
  
  Pd(1:no, 2,6) = -BAx*Pp(1:no,2, 3) + Dp(1:no,5,3)
  Pd(1:no, 3,6) = -BAx*Pp(1:no,3, 3) + Dp(1:no,6,3)
  Pd(1:no, 4,6) = -BAx*Pp(1:no,4, 3) + Dp(1:no,8,3)
  
  Pd(1:no, 2,7) = -BAy*Pp(1:no,2,3) + Dp(1:no,6,3)
  Pd(1:no, 3,7) = -BAy*Pp(1:no,3,3) + Dp(1:no,7,3)
  Pd(1:no, 4,7) = -BAy*Pp(1:no,4,3) + Dp(1:no,9,3)
  
  Pd(1:no, 2,8) = -BAx*Pp(1:no,2,4) + Dp(1:no,5,4)
  Pd(1:no, 3,8) = -BAx*Pp(1:no,3,4) + Dp(1:no,6,4)
  Pd(1:no, 4,8) = -BAx*Pp(1:no,4,4) + Dp(1:no,8,4)
  
  Pd(1:no, 2,9) = -BAy*Pp(1:no,2,4) + Dp(1:no,6,4)
  Pd(1:no, 3,9) = -BAy*Pp(1:no,3,4) + Dp(1:no,7,4)
  Pd(1:no, 4,9) = -BAy*Pp(1:no,4,4) + Dp(1:no,9,4)
  
  Pd(1:no, 2,10) = -BAz*Pp(1:no,2,4) + Dp(1:no,8,4)
  Pd(1:no, 3,10) = -BAz*Pp(1:no,3,4) + Dp(1:no,9,4)
  Pd(1:no, 4,10) = -BAz*Pp(1:no,4,4) + Dp(1:no,10,4)
  
  return
end subroutine calcPd0

!***************************************************


subroutine calcSd0(Sd,Sp,Pp,no,nL,L,BAx,BAy,BAz)
  implicit none
  integer, intent(in)   :: no, nL, L
  real*8, intent(out)   :: Sd(1:no,1,5:10)
  real*8, intent(in)    :: Sp(1:no,1,2:4), Pp(1:no,2:4,2:4), BAx(1:no), BAy(1:no), BAz(1:no)
  
  Sd(1:no,1,5) = -BAx*Sp(1:no,1,2) + Pp(1:no,2,2)
  Sd(1:no,1,6) = -BAy*Sp(1:no,1,2) + Pp(1:no,3,2)
  Sd(1:no,1,7) = -BAy*Sp(1:no,1,3) + Pp(1:no,3,3)
  Sd(1:no,1,8) = -BAx*Sp(1:no,1,4) + Pp(1:no,2,4)
  Sd(1:no,1,9) = -BAy*Sp(1:no,1,4) + Pp(1:no,3,4)
  Sd(1:no,1,10) = -BAz*Sp(1:no,1,4) + Pp(1:no,4,4)
  
  return
end subroutine calcSd0

!***************************************************


subroutine calcSf0(Sf,Sd,Pd,no,nL,L,BAx,BAy,BAz)
  implicit none
  integer, intent(in)   :: no, nL, L
  real*8, intent(out)   :: Sf(1:no,1,11:20)
  real*8, intent(in)    :: Sd(1:no,1,5:10), Pd(1:no,2:4,5:10),  BAx(1:no), BAy(1:no), BAz(1:no)
  
  Sf(1:no,1,11) = -BAx*Sd(1:no,1,5) + Pd(1:no,2,5)
  Sf(1:no,1,12) = -BAx*Sd(1:no,1,6) + Pd(1:no,2,6)
  Sf(1:no,1,13) = -BAx*Sd(1:no,1,7) + Pd(1:no,2,7)
  Sf(1:no,1,14) = -BAy*Sd(1:no,1,7) + Pd(1:no,3,7)
  Sf(1:no,1,15) = -BAx*Sd(1:no,1,8) + Pd(1:no,2,8)
  Sf(1:no,1,16) = -BAx*Sd(1:no,1,9) + Pd(1:no,2,9)
  Sf(1:no,1,17) = -BAy*Sd(1:no,1,9) + Pd(1:no,3,9)
  Sf(1:no,1,18) = -BAx*Sd(1:no,1,10) + Pd(1:no,2,10)
  Sf(1:no,1,19) = -BAy*Sd(1:no,1,10) + Pd(1:no,3,10)
  Sf(1:no,1,20) = -BAz*Sd(1:no,1,10) + Pd(1:no,4,10)
  return
end subroutine calcSf0



!***********************************
subroutine calckeSp(keSp,Sp,Sf,no,beta)
  implicit none
  real*8, intent(out)              :: keSp(1:3*no)
  real*8, intent(in)               :: Sf(1:no,1,11:20),  Sp(1:no,1,2:4),  beta(1:no)
  integer, intent(in)              :: no 
  real*8                           :: pf, pf2
  integer                          :: i
  
  keSp = 0d0;
  
  do i=1,no
     pf = beta(i)*(3d0+2*1);
     keSp(3*i-2) = pf*Sp(i,1,2)
     keSp(3*i-1) = pf*Sp(i,1,3)
     keSp(3*i  ) = pf*Sp(i,1,4)
     
     pf2 = -2d0*beta(i)**2
     keSp(3*i-2) = keSp(3*i-2) + pf2*Sf(i,1,11) 
     keSp(3*i-2) = keSp(3*i-2) + pf2*Sf(i,1,13)
     keSp(3*i-2) = keSp(3*i-2) + pf2*Sf(i,1,18)
     
     keSp(3*i-1) = keSp(3*i-1) + pf2*Sf(i,1,12) 
     keSp(3*i-1) = keSp(3*i-1) + pf2*Sf(i,1,14)
     keSp(3*i-1) = keSp(3*i-1) + pf2*Sf(i,1,19)
     
     keSp(3*i  ) = keSp(3*i  ) + pf2*Sf(i,1,15) 
     keSp(3*i  ) = keSp(3*i  ) + pf2*Sf(i,1,17)
     keSp(3*i  ) = keSp(3*i  ) + pf2*Sf(i,1,20)
     
  end do
  return
end subroutine calckeSp
