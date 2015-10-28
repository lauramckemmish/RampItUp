

!************************************************************************************

subroutine calcgenmpRgSs(mpRgSs,no,RgSs, L, nLCart,nLpure,kmax) 
 implicit none
 integer, intent(in)           :: no,kmax, L, nLCart, nLpure
 real*8, intent(in)            :: RgSs(1:30,1:no)
 real*8, intent(out)           :: mpRgSs(1:no,1:nLpure+4)
 integer                       :: n(1:no), i, K, Lc, Lz, Ly, Lx, m, j
 real*8                        :: CartSs(1:no,1:nLCart,0:L), PureMPDefSs(1:no,1:nLpure), Fund(1:no,0:L)
 real*8                        :: CartMP(1:no,1:nLCart), beta(1:no), BAx(1:no), BAy(1:no), BAz(1:no)
 real*8                        :: T(1:no), h(1:no,0:kmax,0:L), nc
 real*8  ::  chargeSs(1:no,1:4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program calculates pure real multipole moments of the 
!! Ss shell pair up to L-poles.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 mpRgSs = 0d0;

   beta =RgSs(8,1:no) ;
   BAx = RgSs(4,1:no) ;   BAy = RgSs(5,1:no) ;   BAz = RgSs(6,1:no) 
   T = beta * (BAx*BAx + BAy*BAy + BAz*BAz);
   n = int(RgSs(7,1:no));

   call calcmphk(h,kmax,L,T,no)                                   !  Compute h_k(T) and their derivatives
   call calcintmpSs(Fund,h,kmax,L,n,beta,no)                      !  Compute [S|s]^(m) integrals
   call calcmpCartSs(CartSs,Fund,nLCart,L,beta,BAx,BAy,BAz,no)        !  Compute [R|s] integrals

   PureMPDefSs = 0d0;
   call Cart2Pure(PureMPDefSs, CartSs,no,nLCart,nLPure, L)

   mpRgSs(1:no,1:nLPure) = PureMPDefSs;
   mpRgSs(1:no,nLPure+1) = RgSs(1,1:no);
   mpRgSs(1:no,nLPure+2) = RgSs(2,1:no);
   mpRgSs(1:no,nLPure+3) = RgSs(9,1:no);
   mpRgSs(1:no,nLPure+4) = RgSs(3,1:no);

do i=1,no
 if (beta(i) .gt. 10 .and. (BAx(i)*BAx(i) + BAy(i)*BAy(i) + BAz(i)*BAz(i)) .lt. 1d-10) then
  call newgenRgSsol(chargeSs(i,1),1,RgSs(1:30,i),40)
  mpRgSs(i,1)        = chargeSs(i,1)
  mpRgSs(i,2:nLPure) = 0d0;
  mpRgSs(i,nLPure+1:nLpure+4) = (/RgSs(1,i), RgSs(2,i), RgSs(9,i) , RgSs(3,i)/)
 end if
end do

return
end subroutine calcgenmpRgSs

!****************************************************************

subroutine calcmphk(h,kmax,L,T,no)
   implicit none
   integer,intent(in)            :: kmax, L, no
   real*8 ,intent(in)            :: T(1:no)
   real*8 ,intent(out)           :: h(1:no,0:kmax,0:L)
   integer                       :: k,m

   h = 0d0;
   h(1:no,0,0) = exp(-T) * 2
   h(1:no,1,0) = (4*T-6) * h(1:no,0,0)
   do k = 2,kmax                                      !  Compute h_k(T) using forward recuCartSsion
      h(1:no,k,0) = (4*T-8*k+2) * h(1:no,k-1,0) - 8*(k-1)*(2*k-1) * h(1:no,k-2,0)
   end do
   do m = 1,L                                         !  Compute the mth derivatives of h_k(T)
      h(1:no,0,m) = - h(1:no,0,m-1)
      do k = 1,kmax
         h(1:no,k,m) = - h(1:no,k,m-1) - 4*k * h(1:no,k-1,m)
      end do
   end do
   return
end subroutine calcmphk

!************************************************************************************************

subroutine calcintmpSs(Fund,h,kmax,L,n,beta,no)
   implicit none
   integer,intent(in)            :: kmax,L,n(1:no), no
   real*8 ,intent(in)            :: h(1:no,0:kmax,0:L), beta(1:no)
   real*8 ,intent(out)           :: Fund(1:no,0:L)
   integer                       :: k,m
   real*8                        :: coeff(1:no)

   Fund = 0d0
   coeff = 12.56637061435917d0 / (n+1) / beta
   do k = 0,kmax                                      !  Compute [S|s]^(m) from h_k^(m)(T)
      coeff = coeff * beta / ((n+2*k+2)*(n+2*k+3))
      do m = 0,L
         Fund(1:no,m) = Fund(1:no,m) + (k+1)*coeff * h(1:no,k,m)
      end do
   end do
   return
end subroutine calcintmpSs

!************************************************************************************************

subroutine calcmpCartSs(CartSs,Fund,nLCart,L,beta,BAx,BAy,BAz,no)
   implicit none
   integer,intent(in)  :: nLCart, L, no
   real*8 ,intent(in)  :: beta(1:no), BAx(1:no), BAy(1:no), BAz(1:no), Fund(1:no,0:L)
   real*8 ,intent(out) :: CartSs(1:no,1:nLCart,0:L)
   integer             :: m, i
   real*8              :: r2beta(1:no),t000(1:no),t100(1:no),t010(1:no)
   real*8              :: t001(1:no),t200(1:no),t110(1:no),t020(1:no),t002(1:no)
   real*8              :: t300(1:no),t120(1:no),t030(1:no),t201(1:no),t012(1:no),t003(1:no)
   real*8              :: t210(1:no),t102(1:no),t021(1:no),t400(1:no),t040(1:no),t004(1:no)
   real*8              :: t220(1:no),t022(1:no),t202(1:no),t031(1:no),t013(1:no),t310(1:no)
   real*8              :: t500(1:no),t050(1:no),t005(1:no),t320(1:no),t230(1:no),t410(1:no)
   real*8              :: t302(1:no),t130(1:no),t203(1:no),t122(1:no),t140(1:no),t103(1:no)
   real*8              :: t014(1:no),t023(1:no),t032(1:no),t041(1:no), t131(1:no)
   real*8              :: t600(1:no), t420(1:no), t330(1:no), t150(1:no), t060(1:no), t510(1:no), t240(1:no)
   real*8              :: t402(1:no), t132(1:no), t105(1:no), t042(1:no), t033(1:no), t024(1:no)
   real*8              :: t006(1:no), t015(1:no), t123(1:no), t204(1:no), t222(1:no), t303(1:no)
   real*8              :: t700(1:no), t520(1:no), t430(1:no), t250(1:no), t160(1:no), t070(1:no)
   real*8              :: t610(1:no), t502(1:no), t322(1:no), t232(1:no), t142(1:no), t052(1:no)
   real*8              :: t403(1:no), t223(1:no), t043(1:no), t304(1:no)
   real*8              :: t124(1:no), t034(1:no), t205(1:no), t025(1:no), t106(1:no)
   real*8              :: t016(1:no), t007(1:no), t051(1:no), t133(1:no),  t340(1:no)

   CartSs = 0d0;

   r2beta = 1d0 / (2*beta)

   CartSs(1:no,1,0:L) = Fund(1:no,0:L)                                   !  [S|s]^(m)

   do m = L-1,0,-1                                         !  [P|s]^(m)  4L flops
      t000 =  CartSs(1:no, 1,m) + CartSs(1:no, 1,m+1)
      CartSs(1:no, 2,m) = BAx * t000                      !  100
      CartSs(1:no, 3,m) = BAy * t000                      !  010
      CartSs(1:no, 4,m) = BAz * t000                      !  001
   end do

   do m = L-2,0,-1                                         !  [D|s]^(m)  14(L-1) flops
      t000 = (CartSs(1:no, 1,m) + CartSs(1:no, 1,m+1)) * r2beta
      t100 =  CartSs(1:no, 2,m) + CartSs(1:no, 2,m+1)
      t010 =  CartSs(1:no, 3,m) + CartSs(1:no, 3,m+1)
      t001 =  CartSs(1:no, 4,m) + CartSs(1:no, 4,m+1)
      CartSs(1:no, 5,m) = BAx * t100 + t000               !  200
      CartSs(1:no, 6,m) = BAy * t100                      !  110
      CartSs(1:no, 7,m) = BAy * t010 + t000               !  020
      CartSs(1:no, 8,m) = BAx * t001                      !  101
      CartSs(1:no, 9,m) = BAy * t001                      !  011
      CartSs(1:no,10,m) = BAz * t001 + t000               !  002
   end do

   do m = L-3,0,-1                                         !  [F|s]^(m)  26(L-2) flops
      t100 = (CartSs(1:no, 2,m) + CartSs(1:no, 2,m+1)) * r2beta
      t010 = (CartSs(1:no, 3,m) + CartSs(1:no, 3,m+1)) * r2beta
      t001 = (CartSs(1:no, 4,m) + CartSs(1:no, 4,m+1)) * r2beta

      t200 =  CartSs(1:no, 5,m) + CartSs(1:no, 5,m+1)
      t110 =  CartSs(1:no, 6,m) + CartSs(1:no, 6,m+1)
      t020 =  CartSs(1:no, 7,m) + CartSs(1:no, 7,m+1)
      t002 =  CartSs(1:no,10,m) + CartSs(1:no,10,m+1)

      CartSs(1:no,11,m) = BAx * t200 + t100 * 2           !  300
      CartSs(1:no,12,m) = BAy * t200                      !  210
      CartSs(1:no,13,m) = BAx * t020                      !  120
      CartSs(1:no,14,m) = BAy * t020 + t010 * 2           !  030
      CartSs(1:no,15,m) = BAz * t200                      !  201
      CartSs(1:no,16,m) = BAz * t110                      !  111
      CartSs(1:no,17,m) = BAz * t020                      !  021
      CartSs(1:no,18,m) = BAx * t002                      !  102
      CartSs(1:no,19,m) = BAy * t002                      !  012
      CartSs(1:no,20,m) = BAz * t002 + t001 * 2           !  003
   end do

   do m = L-4,0,-1                                         !  [G|s]^(m)  36(L-3) flops
      t200 = (CartSs(1:no, 5,m) + CartSs(1:no, 5,m+1)) * r2beta
      t020 = (CartSs(1:no, 7,m) + CartSs(1:no, 7,m+1)) * r2beta
      t002 = (CartSs(1:no,10,m) + CartSs(1:no,10,m+1)) * r2beta

      t300 =  CartSs(1:no,11,m) + CartSs(1:no,11,m+1)
      t120 =  CartSs(1:no,13,m) + CartSs(1:no,13,m+1)
      t030 =  CartSs(1:no,14,m) + CartSs(1:no,14,m+1)
      t201 =  CartSs(1:no,15,m) + CartSs(1:no,15,m+1)
      t012 =  CartSs(1:no,19,m) + CartSs(1:no,19,m+1)
      t003 =  CartSs(1:no,20,m) + CartSs(1:no,20,m+1)

      CartSs(1:no,21,m) = BAx * t300 + t200 * 3           !  400
      CartSs(1:no,22,m) = BAy * t300                      !  310
      CartSs(1:no,23,m) = BAx * t120 + t020 * 1           !  220
      CartSs(1:no,24,m) = BAx * t030                      !  130
      CartSs(1:no,25,m) = BAy * t030 + t020 * 3           !  040
      CartSs(1:no,26,m) = BAz * t300                      !  301
      CartSs(1:no,27,m) = BAy * t201                      !  211
      CartSs(1:no,28,m) = BAz * t120                      !  121
      CartSs(1:no,29,m) = BAz * t030                      !  031
      CartSs(1:no,30,m) = BAz * t201 + t200 * 1           !  202
      CartSs(1:no,31,m) = BAx * t012                      !  112
      CartSs(1:no,32,m) = BAy * t012 + t002 * 1           !  022
      CartSs(1:no,33,m) = BAx * t003                      !  103
      CartSs(1:no,34,m) = BAy * t003                      !  013
      CartSs(1:no,35,m) = BAz * t003 + t002 * 3           !  004
   end do


   do m = L-5,0,-1                                         !  [H|s]^(m)  
      t300 = (CartSs(1:no,11,m) + CartSs(1:no,11,m+1)) * r2beta
      t210 = (CartSs(1:no,12,m) + CartSs(1:no,12,m+1)) * r2beta
      t120 = (CartSs(1:no,13,m) + CartSs(1:no,13,m+1)) * r2beta
      t030 = (CartSs(1:no,14,m) + CartSs(1:no,14,m+1)) * r2beta
      t201 = (CartSs(1:no,15,m) + CartSs(1:no,15,m+1)) * r2beta
      t021 = (CartSs(1:no,17,m) + CartSs(1:no,17,m+1)) * r2beta
      t102 = (CartSs(1:no,18,m) + CartSs(1:no,18,m+1)) * r2beta
      t012 = (CartSs(1:no,19,m) + CartSs(1:no,19,m+1)) * r2beta
      t003 = (CartSs(1:no,20,m) + CartSs(1:no,20,m+1)) * r2beta

      t400 =  CartSs(1:no,21,m) + CartSs(1:no,21,m+1)
      t310 =  CartSs(1:no,22,m) + CartSs(1:no,22,m+1)
      t220 =  CartSs(1:no,23,m) + CartSs(1:no,23,m+1)
      t040 =  CartSs(1:no,25,m) + CartSs(1:no,25,m+1)
      t031 =  CartSs(1:no,29,m) + CartSs(1:no,29,m+1)
      t202 =  CartSs(1:no,30,m) + CartSs(1:no,30,m+1)
      t022 =  CartSs(1:no,32,m) + CartSs(1:no,32,m+1)
      t013 =  CartSs(1:no,34,m) + CartSs(1:no,34,m+1)
      t004 =  CartSs(1:no,35,m) + CartSs(1:no,35,m+1)

      CartSs(1:no,36,m) = BAx * t400 + t300 * 4          !  500
      CartSs(1:no,37,m) = BAy * t400                     !  410
      CartSs(1:no,38,m) = BAx * t220 + t120 * 2          !  320 
      CartSs(1:no,39,m) = BAy * t220 + t210 * 2          !  230
      CartSs(1:no,40,m) = BAx * t040                     !  140
      CartSs(1:no,41,m) = BAy * t040 + t030 * 4          !  050
      CartSs(1:no,42,m) = BAz * t400                     !  401
      CartSs(1:no,43,m) = BAz * t310                     !  311
      CartSs(1:no,44,m) = BAz * t220                     !  221
      CartSs(1:no,45,m) = BAx * t031                     !  131
      CartSs(1:no,46,m) = BAz * t040                     !  041
      CartSs(1:no,47,m) = BAx * t202 + t102 * 2          !  302
      CartSs(1:no,48,m) = BAy * t202                     !  212
      CartSs(1:no,49,m) = BAx * t022                     !  122
      CartSs(1:no,50,m) = BAy * t022 + t012 * 2          !  032
      CartSs(1:no,51,m) = BAz * t202 + t201 * 2          !  203
      CartSs(1:no,52,m) = BAx * t013                     !  113
      CartSs(1:no,53,m) = BAz * t022 + t021 * 2          !  023
      CartSs(1:no,54,m) = BAx * t004                     !  104
      CartSs(1:no,55,m) = BAy * t004                     !  014
      CartSs(1:no,56,m) = BAz * t004 + t003 * 4          !  005



   end do 

   do m = L-6,0,-1

      t400 = (CartSs(1:no,21,m) + CartSs(1:no,21,m+1)) * r2beta
      t040 = (CartSs(1:no,25,m) + CartSs(1:no,25,m+1)) * r2beta
      t004 = (CartSs(1:no,35,m) + CartSs(1:no,35,m+1)) * r2beta
      t130 = (CartSs(1:no,24,m) + CartSs(1:no,24,m+1)) * r2beta
      t103 = (CartSs(1:no,33,m) + CartSs(1:no,33,m+1)) * r2beta
      t013 = (CartSs(1:no,34,m) + CartSs(1:no,34,m+1)) * r2beta
      t220 = (CartSs(1:no,23,m) + CartSs(1:no,23,m+1)) * r2beta
      t202 = (CartSs(1:no,30,m) + CartSs(1:no,30,m+1)) * r2beta
      t022 = (CartSs(1:no,32,m) + CartSs(1:no,32,m+1)) * r2beta
 
      t500 =  CartSs(1:no,36,m) + CartSs(1:no,36,m+1)
      t050 =  CartSs(1:no,41,m) + CartSs(1:no,41,m+1)
      t005 =  CartSs(1:no,56,m) + CartSs(1:no,56,m+1)
      t410 =  CartSs(1:no,37,m) + CartSs(1:no,37,m+1)
      t140 =  CartSs(1:no,40,m) + CartSs(1:no,40,m+1)
      t230 =  CartSs(1:no,39,m) + CartSs(1:no,39,m+1)
      t320 =  CartSs(1:no,38,m) + CartSs(1:no,38,m+1)
      t302 =  CartSs(1:no,47,m) + CartSs(1:no,47,m+1)
      t014 =  CartSs(1:no,55,m) + CartSs(1:no,55,m+1)
      t023 =  CartSs(1:no,53,m) + CartSs(1:no,53,m+1)
      t203 =  CartSs(1:no,51,m) + CartSs(1:no,51,m+1)
      t122 =  CartSs(1:no,49,m) + CartSs(1:no,49,m+1)
      t032 =  CartSs(1:no,50,m) + CartSs(1:no,50,m+1)


      CartSs(1:no,57,m) = BAx * t500 + t400 * 5       ! 600
      CartSs(1:no,58,m) = BAy * t500                  ! 510
      CartSs(1:no,59,m) = BAx * t320 + t220 * 3       ! 420
      CartSs(1:no,60,m) = BAx * t230 + t130 * 2       ! 330
      CartSs(1:no,61,m) = BAy * t230 + t220 * 3       ! 240
      CartSs(1:no,62,m) = BAx * t050                  ! 150
      CartSs(1:no,63,m) = BAy * t050 + t040 * 5       ! 060
      CartSs(1:no,64,m) = BAz * t500                  ! 501
      CartSs(1:no,65,m) = BAz * t410                  ! 411
      CartSs(1:no,66,m) = BAz * t320                  ! 321
      CartSs(1:no,67,m) = BAz * t230                  ! 231
      CartSs(1:no,68,m) = BAz * t140                  ! 141
      CartSs(1:no,69,m) = BAz * t050                  ! 051
      CartSs(1:no,70,m) = BAx * t302 + t202 * 3       ! 402
      CartSs(1:no,71,m) = BAy * t302                  ! 312
      CartSs(1:no,72,m) = BAx * t122 + t022           ! 222
      CartSs(1:no,73,m) = BAx * t032                  ! 132
      CartSs(1:no,74,m) = BAy * t032 + t022 * 3       ! 042
      CartSs(1:no,75,m) = BAx * t203 + t103 * 2       ! 303
      CartSs(1:no,76,m) = BAy * t203                  ! 213
      CartSs(1:no,77,m) = BAx * t023                  ! 123
      CartSs(1:no,78,m) = BAy * t023 + t013 * 2       ! 033
      CartSs(1:no,79,m) = BAz * t203 + t202 * 3       ! 204
      CartSs(1:no,80,m) = BAx * t014                  ! 114
      CartSs(1:no,81,m) = BAz * t023 + t022 * 3       ! 024
      CartSs(1:no,82,m) = BAx * t005                  ! 105
      CartSs(1:no,83,m) = BAy * t005                  ! 015
      CartSs(1:no,84,m) = BAz * t005 + t004 * 5       ! 006
   end do


 

   do m = L-7,0,-1


      t500 = (CartSs(1:no,36,m) + CartSs(1:no,36,m+1)) * r2beta
      t320 = (CartSs(1:no,38,m) + CartSs(1:no,38,m+1)) * r2beta
      t230 = (CartSs(1:no,39,m) + CartSs(1:no,39,m+1)) * r2beta
      t050 = (CartSs(1:no,41,m) + CartSs(1:no,41,m+1)) * r2beta
      t005 = (CartSs(1:no,56,m) + CartSs(1:no,56,m+1)) * r2beta
      t302 = (CartSs(1:no,47,m) + CartSs(1:no,47,m+1)) * r2beta
      t203 = (CartSs(1:no,51,m) + CartSs(1:no,51,m+1)) * r2beta
      t122 = (CartSs(1:no,49,m) + CartSs(1:no,49,m+1)) * r2beta
      t032 = (CartSs(1:no,50,m) + CartSs(1:no,50,m+1)) * r2beta
      t023 = (CartSs(1:no,53,m) + CartSs(1:no,53,m+1)) * r2beta
      t041 = (CartSs(1:no,46,m) + CartSs(1:no,46,m+1)) * r2beta
      t131 = (CartSs(1:no,45,m) + CartSs(1:no,45,m+1)) * r2beta

      t600 =  CartSs(1:no,57,m) + CartSs(1:no,57,m+1)
      t420 =  CartSs(1:no,59,m) + CartSs(1:no,59,m+1)
      t330 =  CartSs(1:no,60,m) + CartSs(1:no,60,m+1)
      t150 =  CartSs(1:no,62,m) + CartSs(1:no,62,m+1)
      t060 =  CartSs(1:no,63,m) + CartSs(1:no,63,m+1)
      t510 =  CartSs(1:no,58,m) + CartSs(1:no,58,m+1)
      t240 =  CartSs(1:no,61,m) + CartSs(1:no,61,m+1)

      t402 =  CartSs(1:no,70,m) + CartSs(1:no,70,m+1)
      t132 =  CartSs(1:no,73,m) + CartSs(1:no,73,m+1)
      t105 =  CartSs(1:no,82,m) + CartSs(1:no,82,m+1)
      t042 =  CartSs(1:no,74,m) + CartSs(1:no,74,m+1)
      t033 =  CartSs(1:no,78,m) + CartSs(1:no,78,m+1)
      t024 =  CartSs(1:no,81,m) + CartSs(1:no,81,m+1)

      t006 =  CartSs(1:no,84,m) + CartSs(1:no,84,m+1)
      t015 =  CartSs(1:no,83,m) + CartSs(1:no,83,m+1)
      t123 =  CartSs(1:no,77,m) + CartSs(1:no,77,m+1)
      t204 =  CartSs(1:no,79,m) + CartSs(1:no,79,m+1)
      t222 =  CartSs(1:no,72,m) + CartSs(1:no,72,m+1)
      t303 =  CartSs(1:no,75,m) + CartSs(1:no,75,m+1)


      CartSs(1:no,85,m) = BAx * t600 + t500 * 6       ! 700
      CartSs(1:no,86,m) = BAy * t600                  ! 610
      CartSs(1:no,87,m) = BAx * t420 + t320 * 4       ! 520
      CartSs(1:no,88,m) = BAx * t330 + t230 * 3       ! 430
      CartSs(1:no,89,m) = BAy * t330 + t320 * 3       ! 340
      CartSs(1:no,90,m) = BAx * t150 + t050 * 1       ! 250
      CartSs(1:no,91,m) = BAx * t060                  ! 160
      CartSs(1:no,92,m) = BAy * t060 + t050 * 6       ! 070
      CartSs(1:no,93,m) = BAz * t600                  ! 601
      CartSs(1:no,94,m) = BAz * t510                  ! 511
      CartSs(1:no,95,m) = BAz * t420                  ! 421
      CartSs(1:no,96,m) = BAz * t330                  ! 331
      CartSs(1:no,97,m) = BAz * t240                  ! 241
      CartSs(1:no,98,m) = BAz * t150                  ! 151
      CartSs(1:no,99,m) = BAz * t060                  ! 061
      CartSs(1:no,100,m) = BAx * t402 + t302 * 4       ! 502
      CartSs(1:no,101,m) = BAy * t402                  ! 412
      CartSs(1:no,102,m) = BAx * t222 + t122 * 2       ! 322
      CartSs(1:no,103,m) = BAx * t132 + t032 * 1       ! 232
      CartSs(1:no,104,m) = BAx * t042                  ! 142
      CartSs(1:no,105,m) = BAy * t042 + t032 * 4       ! 052
      CartSs(1:no,106,m) = BAx * t303 + t203 * 3       ! 403
      CartSs(1:no,107,m) = BAy * t303                  ! 313
      CartSs(1:no,108,m) = BAx * t123 + t023 * 1       ! 223
      CartSs(1:no,109,m) = BAx * t033                  ! 133
      CartSs(1:no,110,m) = BAz * t042 + t041 * 2       ! 043

      CartSs(1:no,111,m) = BAz * t303 + t302 * 3       ! 304
      CartSs(1:no,112,m) = BAy * t204                  ! 214
      CartSs(1:no,113,m) = BAx * t024                  ! 124
      CartSs(1:no,114,m) = BAz * t033 + t032 * 3       ! 034

      CartSs(1:no,115,m) = BAx * t105 + t005 * 1       ! 205
      CartSs(1:no,116,m) = BAy * t105                  ! 115
      CartSs(1:no,117,m) = BAy * t015 + t005 * 1       ! 025
      CartSs(1:no,118,m) = BAx * t006                  ! 106
      CartSs(1:no,119,m) = BAy * t006                  ! 016
      CartSs(1:no,120,m) = BAz * t006 + t005 * 6       ! 007

   end do


   do m = L-8,0,-1

      t700 =  CartSs(1:no,85,m) + CartSs(1:no,85,m+1)
      t610 =  CartSs(1:no,86,m) + CartSs(1:no,86,m+1)
      t520 =  CartSs(1:no,87,m) + CartSs(1:no,87,m+1)
      t430 =  CartSs(1:no,88,m) + CartSs(1:no,88,m+1)
      t340 =  CartSs(1:no,89,m) + CartSs(1:no,89,m+1)
      t250 =  CartSs(1:no,90,m) + CartSs(1:no,90,m+1)
      t160 =  CartSs(1:no,91,m) + CartSs(1:no,91,m+1)
      t070 =  CartSs(1:no,92,m) + CartSs(1:no,92,m+1)

      t502 =  CartSs(1:no,100,m) + CartSs(1:no,100,m+1)
      t322 =  CartSs(1:no,102,m) + CartSs(1:no,102,m+1)
      t142 =  CartSs(1:no,104,m) + CartSs(1:no,104,m+1)
      t052 =  CartSs(1:no,105,m) + CartSs(1:no,105,m+1)

      t043 =  CartSs(1:no,110,m) + CartSs(1:no,110,m+1)
      t403 =  CartSs(1:no,106,m) + CartSs(1:no,106,m+1)

      t304 =  CartSs(1:no,111,m) + CartSs(1:no,111,m+1)
      t124 =  CartSs(1:no,113,m) + CartSs(1:no,113,m+1)
      t034 =  CartSs(1:no,114,m) + CartSs(1:no,114,m+1)

      t205 =  CartSs(1:no,115,m) + CartSs(1:no,115,m+1)
      t025 =  CartSs(1:no,117,m) + CartSs(1:no,117,m+1)

      t106 =  CartSs(1:no,118,m) + CartSs(1:no,118,m+1)
      t016 =  CartSs(1:no,119,m) + CartSs(1:no,119,m+1)
      t007 =  CartSs(1:no,120,m) + CartSs(1:no,120,m+1)


      t600 = (CartSs(1:no,57,m) + CartSs(1:no,57,m+1)) * r2beta
      t420 = (CartSs(1:no,59,m) + CartSs(1:no,59,m+1)) * r2beta
      t330 = (CartSs(1:no,60,m) + CartSs(1:no,60,m+1)) * r2beta
      t150 = (CartSs(1:no,62,m) + CartSs(1:no,62,m+1)) * r2beta
      t060 = (CartSs(1:no,63,m) + CartSs(1:no,63,m+1)) * r2beta

      t402 = (CartSs(1:no,70,m) + CartSs(1:no,70,m+1)) * r2beta
      t222 = (CartSs(1:no,72,m) + CartSs(1:no,72,m+1)) * r2beta
      t132 = (CartSs(1:no,73,m) + CartSs(1:no,73,m+1)) * r2beta
      t042 = (CartSs(1:no,74,m) + CartSs(1:no,74,m+1)) * r2beta


      t303 = (CartSs(1:no,75,m) + CartSs(1:no,75,m+1)) * r2beta
      t123 = (CartSs(1:no,77,m) + CartSs(1:no,77,m+1)) * r2beta
      t033 = (CartSs(1:no,78,m) + CartSs(1:no,78,m+1)) * r2beta
      t051 = (CartSs(1:no,69,m) + CartSs(1:no,69,m+1)) * r2beta

      t024 = (CartSs(1:no,81,m) + CartSs(1:no,81,m+1)) * r2beta
      t105 = (CartSs(1:no,82,m) + CartSs(1:no,82,m+1)) * r2beta
      t015 = (CartSs(1:no,83,m) + CartSs(1:no,83,m+1)) * r2beta

      t006 = (CartSs(1:no,84,m) + CartSs(1:no,84,m+1)) * r2beta

      CartSs(1:no,121,m) = BAx * t700 + t600 * 7       ! 800
      CartSs(1:no,122,m) = BAy * t700                  ! 710
      CartSs(1:no,123,m) = BAx * t520 + t420 * 5       ! 620
      CartSs(1:no,124,m) = BAx * t430 + t330 * 4       ! 530
      CartSs(1:no,125,m) = BAy * t430 + t420 * 3       ! 440
      CartSs(1:no,126,m) = BAx * t250 + t150 * 2       ! 350
      CartSs(1:no,127,m) = BAx * t160 + t060 * 1       ! 260
      CartSs(1:no,128,m) = BAx * t070                  ! 170
      CartSs(1:no,129,m) = BAy * t070 + t060 * 7       ! 080

      CartSs(1:no,130,m) = BAz * t700                  ! 701
      CartSs(1:no,131,m) = BAz * t610                  ! 611
      CartSs(1:no,132,m) = BAz * t520                  ! 521
      CartSs(1:no,133,m) = BAz * t430                  ! 431
      CartSs(1:no,134,m) = BAz * t340                  ! 341
      CartSs(1:no,135,m) = BAz * t250                  ! 251
      CartSs(1:no,136,m) = BAz * t160                  ! 161
      CartSs(1:no,137,m) = BAz * t070                  ! 071

      CartSs(1:no,138,m) = BAx * t502 + t402 * 5       ! 602
      CartSs(1:no,139,m) = BAy * t502                  ! 512
      CartSs(1:no,140,m) = BAx * t322 + t222 * 3       ! 422
      CartSs(1:no,141,m) = BAx * t232 + t132 * 2       ! 332
      CartSs(1:no,142,m) = BAx * t142 + t042 * 1       ! 242
      CartSs(1:no,143,m) = BAx * t052                  ! 152
      CartSs(1:no,144,m) = BAy * t052 + t042 * 5       ! 062

      CartSs(1:no,145,m) = BAx * t403 + t303 * 4       ! 503
      CartSs(1:no,146,m) = BAy * t403                  ! 413
      CartSs(1:no,147,m) = BAx * t223 + t123 * 2       ! 323
      CartSs(1:no,148,m) = BAx * t133 + t033 * 1       ! 233
      CartSs(1:no,149,m) = BAx * t043                  ! 143
      CartSs(1:no,150,m) = BAz * t052 + t051 * 2       ! 053

      CartSs(1:no,151,m) = BAz * t403 + t402 * 3       ! 404
      CartSs(1:no,152,m) = BAy * t304                  ! 314
      CartSs(1:no,153,m) = BAx * t124 + t024 * 1       ! 224
      CartSs(1:no,154,m) = BAx * t034                  ! 134
      CartSs(1:no,155,m) = BAz * t043 + t042 * 3       ! 044


      CartSs(1:no,156,m) = BAx * t205 + t105 * 2       ! 305
      CartSs(1:no,157,m) = BAy * t205                  ! 215
      CartSs(1:no,158,m) = BAx * t025                  ! 125
      CartSs(1:no,159,m) = BAy * t025 + t015 * 2       ! 035

      CartSs(1:no,160,m) = BAx * t106 + t006 * 1       ! 206
      CartSs(1:no,161,m) = BAy * t106                  ! 116
      CartSs(1:no,162,m) = BAy * t016 + t006 * 1       ! 026

      CartSs(1:no,163,m) = BAx * t007                  ! 107
      CartSs(1:no,164,m) = BAy * t007                  ! 017

      CartSs(1:no,165,m) = BAz * t007 + t006 * 7       ! 008
 
    end do


      return
end subroutine calcmpCartSs

!**********************************************************************


!****************************************************
subroutine Cart2Pure(PureMPDefSs, CartSs,no,nLCart,nLPure, L)
   implicit none
   integer, intent(in)   :: L, no, nLCart, nLPure
   real*8, intent(in)    :: CartSs(1:no,1:nLCart,0:L)
   real*8, intent(out)   :: PureMPDefSs(1:no,1:nLPure)
   real*8                :: t(1:no), invfact(0:2*L), rsq2
   integer               :: fact(0:2*L), lindex, mindex, cnt, i, x
   real*8                :: CY500(1:no), CY410(1:no), CY320(1:no), CY230(1:no), CY140(1:no)
   real*8                :: CY050(1:no), CY401(1:no), CY311(1:no), CY221(1:no), CY131(1:no)
   real*8                :: CY041(1:no), CY302(1:no), CY212(1:no), CY122(1:no), CY032(1:no)
   real*8                :: CY203(1:no), CY113(1:no), CY023(1:no), CY104(1:no), CY014(1:no), CY005(1:no)
   real*8                :: CY220(1:no), CY202(1:no), CY022(1:no), CY400(1:no), CY040(1:no), CY004(1:no)
   real*8                :: CY600(1:no),CY510(1:no),CY420(1:no),CY330(1:no),CY240(1:no),CY150(1:no),CY060(1:no)
   real*8                :: CY501(1:no),CY411(1:no),CY321(1:no),CY231(1:no),CY141(1:no),CY051(1:no)
   real*8  :: CY402(1:no),CY312(1:no),CY222(1:no),CY132(1:no),CY042(1:no)
   real*8  :: CY303(1:no),CY213(1:no),CY123(1:no),CY033(1:no)
   real*8  :: CY204(1:no),CY114(1:no),CY024(1:no)
   real*8  :: CY105(1:no),CY015(1:no),CY006(1:no)
   real*8  :: CY700(1:no), CY610(1:no), CY520(1:no), CY430(1:no), CY340(1:no), CY250(1:no), CY160(1:no), CY070(1:no)  
   real*8  :: CY601(1:no), CY511(1:no), CY421(1:no), CY331(1:no), CY241(1:no), CY151(1:no), CY061(1:no)
   real*8  :: CY502(1:no), CY412(1:no), CY322(1:no), CY232(1:no), CY142(1:no), CY052(1:no)
   real*8  :: CY403(1:no), CY313(1:no), CY223(1:no), CY133(1:no), CY043(1:no)
   real*8  :: CY304(1:no), CY214(1:no), CY124(1:no), CY034(1:no)
   real*8  :: CY205(1:no), CY115(1:no), CY025(1:no)
   real*8  :: CY106(1:no), CY016(1:no), CY007(1:no)



   PureMPDefSs(1:no,1) =  CartSs(1:no,1,0)

if (nLPure == 1) then; return; end if

   PureMPDefSs(1:no,2) = -0.5d0*CartSs(1:no,3,0)
   PureMPDefSs(1:no,3) =        CartSs(1:no,4,0)
   PureMPDefSs(1:no,4) = -0.5d0*CartSs(1:no,2,0)

if (nLPure == 4) then; return; end if

   PureMPDefSs(1:no,5) =  0.25d0* CartSs(1:no,6,0)
   PureMPDefSs(1:no,6) = - 0.5d0* CartSs(1:no,9,0)
   PureMPDefSs(1:no,7) =  0.25d0*(CartSs(1:no,10,0)*2-CartSs(1:no,5,0)-CartSs(1:no,7,0))
   PureMPDefSs(1:no,8) = - 0.5d0* CartSs(1:no,8,0)
   PureMPDefSs(1:no,9) = 0.125d0*(CartSs(1:no,5,0) - CartSs(1:no,7,0))

if (nLPure == 9) then; return; end if

   PureMPDefSs(1:no,10) = -0.020833333333333333333*(3*CartSs(1:no,12,0) - CartSs(1:no,14,0)) !l=3, m=-3
   PureMPDefSs(1:no,11) =  0.25d0*CartSs(1:no,16,0) !l=3, m=-2
   PureMPDefSs(1:no,12) = 0.0625d0*(CartSs(1:no,12,0)+CartSs(1:no,14,0) - 4d0*CartSs(1:no,19,0)) !l=3, m=-1
   PureMPDefSs(1:no,13) = 0.083333333333333333333d0*(-3d0*(CartSs(1:no,15,0)+CartSs(1:no,17,0)) &
           + 2d0*CartSs(1:no,20,0)) !l=3,m=0
   PureMPDefSs(1:no,14) = 0.0625d0*(CartSs(1:no,11,0)+CartSs(1:no,13,0)-4d0*CartSs(1:no,18,0))  !l=3, m=1
   PureMPDefSs(1:no,15) = 0.125d0*(CartSs(1:no,15,0)-CartSs(1:no,17,0)) !l=3, m=2
   PureMPDefSs(1:no,16) =  -0.020833333333333333333d0*(CartSs(1:no,11,0) - 3d0*CartSs(1:no,13,0)) !l=3, m=3

if (nLPure == 16) then; return; end if

   PureMPDefSs(1:no,17) = 0.010416666666666666667d0*(CartSs(1:no,22,0) - CartSs(1:no,24,0)) !l=4, m=-4
   PureMPDefSs(1:no,18) = -0.020833333333333333333d0*(3*CartSs(1:no,27,0) - CartSs(1:no,29,0)) !l=4, m=-3

   PureMPDefSs(1:no,19) = -0.020833333333333333333d0*(CartSs(1:no,22,0)+CartSs(1:no,24,0)-6d0*CartSs(1:no,31,0))  !l=4, m=-2
   PureMPDefSs(1:no,20) = 0.020833333333333333333d0*(3*CartSs(1:no,27,0)+3d0*CartSs(1:no,29,0) - 4*CartSs(1:no,34,0)) !l=4,m=-1

 CY220 = CartSs(1:no,23,0)
 CY202 = CartSs(1:no,30,0)
 CY022 = CartSs(1:no,32,0)
 CY400 = CartSs(1:no,21,0)
 CY040 = CartSs(1:no,25,0)
 CY004 = CartSs(1:no,35,0)

   PureMPDefSs(1:no,21) = 0.0052083333333333333333333d0*(8*CY004 - 24*CY022 + 3*CY040 -24*CY202 + 6*CY220 + 3*CY400) !l=4, m=0

   PureMPDefSs(1:no,22) = 0.020833333333333333333d0*(3*CartSs(1:no,26,0)+3d0*CartSs(1:no,28,0) -4*CartSs(1:no,33,0)) !l=4, m=1
   PureMPDefSs(1:no,23) = -0.010416666666666666667d0*(CartSs(1:no,21,0) -CartSs(1:no,25,0) &
       - 6*CartSs(1:no,30,0) + 6*CartSs(1:no,32,0)) !l=4, m=2  
   PureMPDefSs(1:no,24) = -0.020833333333333333333d0*(CartSs(1:no,26,0)-3d0*CartSs(1:no,28,0))!l=4, m=3
   PureMPDefSs(1:no,25) = 0.0026041666666666666667d0*(CartSs(1:no,21,0)-6d0*CartSs(1:no,23,0)+CartSs(1:no,25,0)) !l=4, m=4 

if (nLPure == 25) then; return; end if

CY500 = CartSs(1:no,36,0)
CY410 = CartSs(1:no,37,0)
CY320 = CartSs(1:no,38,0)
CY230 = CartSs(1:no,39,0)
CY140 = CartSs(1:no,40,0)
CY050 = CartSs(1:no,41,0)
CY401 = CartSs(1:no,42,0)

CY311 = CartSs(1:no,43,0)

CY221 = CartSs(1:no,44,0)

CY131 = CartSs(1:no,45,0)

CY041 = CartSs(1:no,46,0)
CY302 = CartSs(1:no,47,0)
CY212 = CartSs(1:no,48,0)
CY122 = CartSs(1:no,49,0)
CY032 = CartSs(1:no,50,0)
CY203 = CartSs(1:no,51,0)

CY113 = CartSs(1:no,52,0)

CY023 = CartSs(1:no,53,0)
CY104 = CartSs(1:no,54,0)
CY014 = CartSs(1:no,55,0)
CY005 = CartSs(1:no,56,0)

   PureMPDefSs(1:no,26) = -0.00026041666666666666667d0*(5*CY410 - 10*CY230 +   CY050 )           
   PureMPDefSs(1:no,27) =    0.010416666666666666667d0*(  CY311 -    CY131           )             
   PureMPDefSs(1:no,28) =   0.0013020833333333333333d0*(8*CY032 -    CY050 - 24*CY212 + 2*CY230 + 3*CY410)            
   PureMPDefSs(1:no,29) =   0.0208333333333333333333d0*(2*CY113 -    CY131 -    CY311)            
   PureMPDefSs(1:no,30) =  -0.0026041666666666666667d0*(8*CY014 - 12*CY032 +    CY050 -12*CY212 + 2*CY230 +   CY410)            
   PureMPDefSs(1:no,31) =   0.0010416666666666666667d0*(8*CY005 - 40*CY023 + 15*CY041 -40*CY203 +30*CY221 +15*CY401)            
   PureMPDefSs(1:no,32) =  -0.0026041666666666666667d0*(8*CY104 - 12*CY122 +    CY140 -12*CY302 + 2*CY320 +   CY500)              
   PureMPDefSs(1:no,33) =  -0.0104166666666666666667d0*(2*CY023 -    CY041 -  2*CY203 +   CY401)           
   PureMPDefSs(1:no,34) =   0.0013020833333333333333d0*(  CY500 + 24*CY122 -  3*CY140 - 8*CY302 - 2*CY320)           
   PureMPDefSs(1:no,35) =   0.0026041666666666666667d0*(  CY401 -  6*CY221 +   CY041  )              
   PureMPDefSs(1:no,36) = -0.00026041666666666666667d0*(5*CY140 - 10*CY320 +   CY500  )  

if (nLPure == 36) then; return; end if

 CY600 = CartSs(1:no,57,0)
 CY510 = CartSs(1:no,58,0)
 CY420 = CartSs(1:no,59,0)
 CY330 = CartSs(1:no,60,0)
 CY240 = CartSs(1:no,61,0)
 CY150 = CartSs(1:no,62,0)
 CY060 = CartSs(1:no,63,0)
 CY501 = CartSs(1:no,64,0)
 CY411 = CartSs(1:no,65,0)
 CY321 = CartSs(1:no,66,0)
 CY231 = CartSs(1:no,67,0)
 CY141 = CartSs(1:no,68,0)
 CY051 = CartSs(1:no,69,0)
 CY402 = CartSs(1:no,70,0)
 CY312 = CartSs(1:no,71,0)
 CY222 = CartSs(1:no,72,0)
 CY132 = CartSs(1:no,73,0)
 CY042 = CartSs(1:no,74,0)
 CY303 = CartSs(1:no,75,0)
 CY213 = CartSs(1:no,76,0)
 CY123 = CartSs(1:no,77,0)
 CY033 = CartSs(1:no,78,0)
 CY204 = CartSs(1:no,79,0)
 CY114 = CartSs(1:no,80,0)
 CY024 = CartSs(1:no,81,0)
 CY105 = CartSs(1:no,82,0)
 CY015 = CartSs(1:no,83,0)
 CY006 = CartSs(1:no,84,0)           

   PureMPDefSs(1:no,37) = 0.0000434027777777777778d0*( 3*CY150 - 10*CY330 +  3*CY510    )     
   PureMPDefSs(1:no,38) = -0.000260416666666666667d0*(   CY051 - 10*CY231 +  5*CY411     )           
   PureMPDefSs(1:no,39) = -0.000520833333333333333d0*(10*CY132 -    CY150 - 10*CY312 +    CY510     )            
   PureMPDefSs(1:no,40) = -0.000434027777777777778d0*(24*CY213 -  8*CY033 +  3*CY051 -  6*CY231 -  9*CY411)           
   PureMPDefSs(1:no,41) = -0.000651041666666666667d0*(16*CY132 -    CY150 + 16*CY312 -  2*CY330 -    CY510 - 16*CY114)            
   PureMPDefSs(1:no,42) = -0.000520833333333333333d0*( 8*CY015 - 20*CY033 +  5*CY051 - 20*CY213 + 10*CY231 +  5*CY411)               
   PureMPDefSs(1:no,43) = 0.0000868055555555555556d0*(16*CY006 -120*CY024 + 90*CY042 -  5*CY060 -120*CY204 +180*CY222 &
                               -15* CY240 + 90*CY402 -15*CY420 - 5*CY600)            
   PureMPDefSs(1:no,44) = -0.000520833333333333333d0*( 8*CY105 - 20*CY123 +  5*CY141 - 20*CY303 + 10*CY321 +  5*CY501)            
   PureMPDefSs(1:no,45) = -0.000325520833333333333d0*(16*CY024 - 16*CY042 +    CY060 - 16*CY204 +    CY240 + 16*CY402 &
                               - CY420 - CY600)             
   PureMPDefSs(1:no,46) = 0.0004340277777777777778d0*(24*CY123 -  9*CY141 -  8*CY303 -  6*CY321 +  3*CY501)           
   PureMPDefSs(1:no,47) = 0.0001302083333333333333d0*(10*CY042 -    CY060 - 60*CY222 +  5*CY240 + 10*CY402 &
                               + 5*CY420 - CY600)           
   PureMPDefSs(1:no,48) = -0.000260416666666666667d0*(5*CY141 - 10*CY321 +    CY501        )             
   PureMPDefSs(1:no,49) =  -0.00002170138888888889d0*(  CY060 - 15*CY240 + 15*CY420 -  CY600)                 

 
if (nLPure == 49) then; return; end if

 CY700 = CartSs(1:no,85,0)
 CY610 = CartSs(1:no,86,0)
 CY520 = CartSs(1:no,87,0)
 CY430 = CartSs(1:no,88,0)
 CY340 = CartSs(1:no,89,0)
 CY250 = CartSs(1:no,90,0)
 CY160 = CartSs(1:no,91,0)
 CY070 = CartSs(1:no,92,0)
 CY601 = CartSs(1:no,93,0)
 CY511 = CartSs(1:no,94,0)
 CY421 = CartSs(1:no,95,0)
 CY331 = CartSs(1:no,96,0)
 CY241 = CartSs(1:no,97,0)
 CY151 = CartSs(1:no,98,0)
 CY061 = CartSs(1:no,99,0)
 CY502 = CartSs(1:no,100,0)
 CY412 = CartSs(1:no,101,0)
 CY322 = CartSs(1:no,102,0)
 CY232 = CartSs(1:no,103,0)
 CY142 = CartSs(1:no,104,0)
 CY052 = CartSs(1:no,105,0)
 CY403 = CartSs(1:no,106,0)
 CY313 = CartSs(1:no,107,0)
 CY223 = CartSs(1:no,108,0)
 CY133 = CartSs(1:no,109,0)
 CY043 = CartSs(1:no,110,0)
 CY304 = CartSs(1:no,111,0)
 CY214 = CartSs(1:no,112,0)
 CY124 = CartSs(1:no,113,0)
 CY034 = CartSs(1:no,114,0)
 CY205 = CartSs(1:no,115,0)
 CY115 = CartSs(1:no,116,0)
 CY025 = CartSs(1:no,117,0)
 CY106 = CartSs(1:no,118,0)           
 CY016 = CartSs(1:no,119,0)           
 CY007 = CartSs(1:no,120,0)           

   PureMPDefSs(1:no,50) = -1.5500992063492063492d-6*( 7*CY610 -35*CY430 + 21*CY250 - CY070)
   PureMPDefSs(1:no,51) = 0.000043402777777777777d0*( 3*CY511 -10*CY331 + 3*CY151 ) 
   PureMPDefSs(1:no,52) = 0.000010850694444444444d0*( 5*CY610+120*CY232 -60*CY412 - 5*CY430 - 12*CY052 &
                               - 9*CY250 + CY070)
   PureMPDefSs(1:no,53) = 0.000173611111111111111d0*(10*CY313 - 3*CY511 - 10*CY133 + 3*CY151)
   PureMPDefSs(1:no,54) = 0.000010850694444444444d0*(-240*CY214-180*CY412 -9*CY610 +80*CY032 & 
                            + 120*CY232 - 15*CY430 - 60*CY052 - 3*CY250 + 3*CY070)
   PureMPDefSs(1:no,55) = 0.000043402777777777777d0*(48*CY115 - 80*CY313 + 15*CY511 - 80*CY133 + 30*CY331 + 15*CY151)       
   PureMPDefSs(1:no,56) = 0.000010850694444444444444d0*(-64*CY016 + 240*CY214 - 120*CY412 + 5*CY610 + 240*CY034 - 240*CY232 &
                                    + 15*CY430 - 120*CY052 + 15*CY250 + 5*CY070)
   PureMPDefSs(1:no,57) = 0.000012400793650793650794d0*( 16*CY007 - 168*CY205 + 210*CY403 - 35*CY601 &
                                   -168*CY025 + 420*CY223 - 105*CY421 + 210*CY043 - 105*CY241 - 35*CY061)
   PureMPDefSs(1:no,58) = 0.000010850694444444444d0*(-64*CY106 + 240*CY304 - 120*CY502 + 5*CY700 + 240*CY124 - 240*CY322 & 
                              + 15*CY520 - 120*CY142 + 15*CY340 + 5*CY160)
   PureMPDefSs(1:no,59) = 0.000021701388888888888d0*(48*CY205 - 80*CY403 + 15*CY601 - 48*CY025 + 15*CY421 &
                              + 80*CY043 - 15*CY241 - 15*CY061)
   PureMPDefSs(1:no,60) = 0.000010850694444444444d0*(-80*CY304 + 60*CY502 - 3*CY700 + 240*CY124 &
                            -120*CY322 + 3*CY520 - 180*CY142 + 15*CY340 + 9*CY160)          
   PureMPDefSs(1:no,61) = 0.000043402777777777777d0*(10*CY403-3*CY601-60*CY223+15*CY421 + 10*CY043 &
                            +15*CY241 - 3*CY061)
   PureMPDefSs(1:no,62) = 0.000010850694444444444d0*(CY700 - 12*CY502 +120*CY322 - 9*CY520 &
                                    -60*CY142 - 5*CY340 + 5*CY160)
   PureMPDefSs(1:no,63) = 0.000021701388888888888d0*(CY601 - 15*CY421 + 15*CY241 -   CY061)
   PureMPDefSs(1:no,64) = -1.55009920634920634921d0*(CY700 - 21*CY520 + 35*CY340 - 7*CY160)


if (nLPure == 64) then; return; end if

  
return
end subroutine Cart2Pure

!*********************************************************************************
