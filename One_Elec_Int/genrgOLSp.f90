!************************************************************************************

subroutine genRgSpol(OLRgSp,no,RgSp,kmax)
  implicit none
  integer,parameter                ::  L = 1, nL = (L+1)*(L+2)*(L+3)/6
  integer, intent(in)              :: no, kmax
  real*8, intent(in)               :: RgSp(1:30,1:no)
  real*8, intent(out)              :: OLRgSp(1:3*no,1:4)
  real*8                           :: beta(1:no),  BAx(1:no), BAy(1:no), BAz(1:no), T(1:no)
  real*8                           :: Rs(1:no,nL,0:L), Fund(1:no,0:L), h(1:no,0:kmax,0:L)
  integer                          :: n(1:no), i
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates (S|s) using the general COLD-like algorithm (truncation of expression)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  OLRgSP = 0d0;
  do i=1,no
     OLRgSP(3*i-2:3*i,2) = RgSp(1,i)
     OLRgSP(3*i-2,3) = RgSp(2,i)
     OLRgSP(3*i-1,3) = RgSp(2,i)+1
     OLRgSP(3*i  ,3) = RgSp(2,i)+2
     OLRgSP(3*i-2:3*i,4) = RgSp(9,i)
  end do
  
  beta = RgSp(8,1:no) 
  BAx = RgSp(4,1:no) 
  BAy = RgSp(5,1:no) 
  BAz = RgSp(6,1:no) 
  T = beta * (BAx*BAx + BAy*BAy + BAz*BAz)
  n = int(RgSp(7,1:no))
  
  call calcSphk(no,h,kmax,L,T)                              !  Compute h_k(T) and their derivatives
  call calcgenRgolSs(Fund,h,kmax,L,n,beta,no)                         !  Compute [S|s]^(m) integrals
  call calcgenRgolPs(Rs,Fund,nL,L,beta,BAx,BAy,BAz,no)         !  Compute [P|s] integrals
  call calcgenRgolSp(OLRgSp,Rs,L,nL,BAx,BAy,BAz,no)            !  Compute [S|p] integrals
  
  return
end subroutine genRgSpol

!************************************************************************************************

subroutine calcSphk(no,h,kmax,L,T)
  implicit none
  integer,intent(in)  :: kmax, L, no
  real*8 ,intent(in)  :: T(1:no)
  real*8 ,intent(out) :: h(1:no,0:kmax,0:L)
  integer             :: k,m
  
  h(1:no,0,0) = exp(-T) * 2
  h(1:no,1,0) = (4*T-6) * h(1:no,0,0)
  do k = 2,kmax                                      !  Compute h_k(T) using forward recursion
     h(1:no,k,0) = (4*T-8*k+2) * h(1:no,k-1,0) - 8*(k-1)*(2*k-1) * h(1:no,k-2,0)
  end do
  do m = 1,L                                         !  Compute the mth derivatives of h_k(T)
     h(1:no,0,m) = - h(1:no,0,m-1)
     do k = 1,kmax
        h(1:no,k,m) = - h(1:no,k,m-1) - 4*k * h(1:no,k-1,m)
     end do
  end do
  return
end subroutine calcSphk

!************************************************************************************************

subroutine calcgenRgolSs(Fund,h,kmax,L,n,beta,no)
  implicit none
  integer,intent(in)  :: kmax,L,n(1:no), no
  real*8 ,intent(in)  :: h(1:no,0:kmax,0:L), beta(1:no)
  real*8 ,intent(out) :: Fund(1:no,0:L)
  integer             :: k,m
  real*8              :: coeff(1:no)
  
  Fund = 0d0
  coeff = 4 * dacos(-1d0) / (n+1) / beta
  do k = 0,kmax                                      !  Compute [S|s]^(m) from h_k^(m)(T)
     coeff = coeff * beta / ((n+2*k+2)*(n+2*k+3))
     do m = 0,L
        Fund(1:no,m) = Fund(1:no,m) + (k+1)*coeff * h(1:no,k,m)
     end do
  end do
  return
end subroutine calcgenRgolSs

!************************************************************************************************

subroutine calcgenRgolPs(Rs,Fund,nL,L,beta,BAx,BAy,BAz,no)
  implicit none
  integer,intent(in)  :: nL, L , no
  real*8 ,intent(in)  :: beta(1:no), BAx(1:no), BAy(1:no), BAz(1:no), Fund(1:no,0:L)
  real*8 ,intent(out) :: Rs(1:no,nL,0:L)
  integer             :: m, i
  real*8              :: r2beta(1:no),t000(1:no),t100(1:no),t010(1:no),t001(1:no)
  
  r2beta = 1d0 / (2*beta)
  
  Rs(1:no,1,0:L) = Fund(1:no,0:L)                                   !  [S|s]^(m)
  do m = L-1,0,-1                                         !  [P|s]^(m)  4L flops
     t000 =  Rs(1:no, 1,m) + Rs( 1:no,1,m+1)
     Rs(1:no, 2,m) = BAx * t000                      !  100
     Rs(1:no, 3,m) = BAy * t000                      !  010
     Rs(1:no, 4,m) = BAz * t000                      !  001
  end do
  
  return
end subroutine calcgenRgolPs

!**************************************************************************************************

subroutine calcgenRgolSp(OLRgSp,Rs,L,nL,BAx,BAy,BAz,no)
  implicit none
  integer, intent(in) :: nL, L, no
  real*8, intent(in)  :: BAx(1:no), BAy(1:no), BAz(1:no), Rs(1:no,nL,0:L)
  real*8, intent(out) :: OLRgSp(1:3*no,1:4)
  integer             :: i
  
  
  do i=1,no
     OLRgSP(3*i-2,1) = -BAx(i)*Rs(i,1,0)+ Rs(i,2,0)              ! [S|px]
     OLRgSP(3*i-1,1) = -BAy(i)*Rs(i,1,0)+ Rs(i,3,0)      ! [S|py]
     OLRgSP(3*i  ,1) = -BAz(i)*Rs(i,1,0)+ Rs(i,4,0)      ! [S|pz]
  end do
  
  return
end subroutine calcgenRgolSp

!**************************************************************************************************
