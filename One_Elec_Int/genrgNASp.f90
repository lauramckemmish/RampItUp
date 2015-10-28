!************************************************************************************

subroutine genrgSpNA(NArgSp,no,RGSpr, noatoms,kmax)
   implicit none
   integer,parameter             :: L = 1, nL = (L+1)*(L+2)*(L+3)/6
   integer, intent(in)           :: no, noatoms, kmax
   real*8, intent(in)            :: RGSpr(1:30,1:no)
   real*8, intent(out)           :: NArgSp(1:3*no,1:4,1:noatoms)
   real*8                        :: NAint(1:3*no), Rs(1:no,nL,0:L), Fund(1:no,0:L), h(1:no,0:kmax,0:L)
   real*8                        :: beta(1:no),  BAx(1:no), BAy(1:no), BAz(1:no), T(1:no)
   integer                       :: n(1:no), i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program calculates (S|s) using the general CNAD-like algorithm (truncation of expression)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   NArgSP = 0d0;

  do i=1,no
   NArgSP(3*i-2:3*i,2,1:noatoms) = RGSpr(1,i)
   NArgSP(3*i-2,3,1:noatoms)     = RGSpr(2,i)
   NArgSP(3*i-1,3,1:noatoms)     = RGSpr(2,i)+1
   NArgSP(3*i,3,1:noatoms)       = RGSpr(2,i)+2
   NArgSP(3*i-2:3*i,4,1:noatoms) = RGSpr(9,i)
  end do

   beta = RGSpr(8,1:no) 
   BAx  = RGSpr(4,1:no) 
   BAy  = RGSpr(5,1:no) 
   BAz  = RGSpr(6,1:no) 
   T    = beta * (BAx*BAx + BAy*BAy + BAz*BAz)  
   n    = int(RGSpr(7,1:no))

   call calcnaSphk(no,h,kmax,L,T)                              !  Compute h_k(T) and their derivatives
   call calcgenrgNASs(Fund,h,kmax,L,n,beta,no)                         !  Compute [S|s]^(m) integrals
   call calcgenrgNAPs(Rs,Fund,nL,L,beta,BAx,BAy,BAz,no)         !  Compute [P|s] integrals
   call calcgenrgNASp(NAint,Rs,L,nL,BAx,BAy,BAz,no)            !  Compute [S|p] integrals

 do i=1,no
  NargSP(3*i-2,1,int(RGSpr(3,i))) = NAint(3*i-2)
  NargSP(3*i-1,1,int(RGSpr(3,i))) = NAint(3*i-1)
  NargSP(3*i  ,1,int(RGSpr(3,i))) = NAint(3*i  )
 end do
 
return
end subroutine genrgSpNA

!************************************************************************************************

subroutine calcnaSphk(no,h,kmax,L,T)
   implicit none
   integer,intent(in)  :: kmax, L
   real*8 ,intent(in)  :: T(1:no)
   real*8 ,intent(out) :: h(1:no,0:kmax,0:L)
   integer             :: k,m, no

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
end subroutine calcnaSphk

!************************************************************************************************

subroutine calcgenrgNASs(Fund,h,kmax,L,n,beta,no)
   implicit none
   integer,intent(in)  :: kmax,L,n(1:no), no
   real*8 ,intent(in)  :: h(1:no,0:kmax,0:L), beta(1:no)
   real*8 ,intent(out) :: Fund(1:no,0:L)
   integer             :: k,m
   real*8              :: coeff(1:no)

   Fund = 0d0
   coeff = 6.283185307179586d0 / beta
   do k = 0,kmax                                     
      coeff = coeff * beta / ((n+2d0*k+1d0)*(n+2d0*k+2d0))
      do m=0,L
        Fund(1:no,m) = Fund(1:no,m) +  coeff * h(1:no,k,m)
      end do
   end do
   return
end subroutine calcgenrgNASs

!************************************************************************************************

subroutine calcgenrgNAPs(Rs,Fund,nL,L,beta,BAx,BAy,BAz,no)
   implicit none
   integer,intent(in)  :: nL, L, no
   real*8 ,intent(in)  :: beta(1:no), BAx(1:no), BAy(1:no), BAz(1:no), Fund(1:no,0:L)
   real*8 ,intent(out) :: Rs(1:no,nL,0:L)
   integer             :: m
   real*8              :: r2beta(1:no), t000(1:no)

   r2beta = 1d0 / (2*beta)

   Rs(1:no,1,0:L) = Fund(1:no,0:L)                                   !  [S|s]^(m)
   do m = L-1,0,-1                                         !  [P|s]^(m)  4L flops
      t000 =  Rs(1:no, 1,m) + Rs( 1:no,1,m+1)
      Rs(1:no, 2,m) = BAx * t000                      !  100
      Rs(1:no, 3,m) = BAy * t000                      !  010
      Rs(1:no, 4,m) = BAz * t000                      !  001
   end do

return
end subroutine calcgenrgNAPs

!**************************************************************************************************

subroutine calcgenrgNASp(NArgSp,Rs,L,nL,BAx,BAy,BAz,no)
 implicit none
 integer, intent(in) :: nL, L, no
 real*8, intent(in)  :: BAx(1:no), BAy(1:no), BAz(1:no),  Rs(1:no,nL,0:L)
 real*8, intent(out) :: NArgSp(1:3*no)
 integer             :: i


 do i=1,no
  NArgSP(3*i-2) = -BAx(i)*Rs(i,1,0)+ Rs(i,2,0)              ! [S|px]
  NArgSP(3*i-1) = -BAy(i)*Rs(i,1,0)+ Rs(i,3,0)              ! [S|py]
  NArgSP(3*i  ) = -BAz(i)*Rs(i,1,0)+ Rs(i,4,0)              ! [S|pz]
 end do


return
end subroutine calcgenrgNASp