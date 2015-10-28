subroutine calcasympRgSp(mpRgSp,no,RgSp, nLpure) 
  implicit none
  integer, intent(in)           :: no, nLpure
  real*8, intent(in)            :: RgSp(1:30,1:no)
  real*8, intent(out)           :: mpRgSp(1:3*no,1:nLpure+4)
  real*8                        :: dipoleSp(1:no)
  integer                       :: i
  
  mpRgSp = 0d0;
  
  call asyRgSpdip(dipoleSp,no,RgSp)
  
  do i=1,no
     mpRgSp(3*i-2      , 4)       = -0.5d0*dipoleSp(i)
     mpRgSp(3*i-1      , 2)       = -0.5d0*dipoleSp(i)
     mpRgSp(3*i-0      , 3)       =        dipoleSp(i)
     mpRgSp(3*i-2:3*i  ,nLPure+1) = RgSp(1,i)
     mpRgSp(3*i-2      ,nLPure+2) = RgSp(2,i)
     mpRgSp(3*i-1      ,nLPure+2) = RgSp(2,i)+1
     mpRgSp(3*i-0      ,nLPure+2) = RgSp(2,i)+2
     mpRgSp(3*i-2:3*i  ,nLPure+3) = RgSp(9,i)
     mpRgSp(3*i-2:3*i  ,nLPure+4) = RgSp(3,i)
  end do
  
  return
end subroutine calcasympRgSp

!*****************************************
subroutine calcgenmpRgSp(mpRgSp,no,RgSp, L, nLCart,nLpure,kmax)
  implicit none
  integer, intent(in)           :: no, kmax,L, nLCart, nLpure
  real*8, intent(in)            :: RgSp(1:30,1:no)
  real*8, intent(out)           :: mpRgSp(1:3*no,1:nLpure+4)
  real*8                        :: RgSpy(1:30,1:no), RgSpz(1:30,1:no), CartSp(1:3*no,1:nLCart), CartMP(1:no,1:nLCart)
  real*8                        :: Rs(1:no,1:(L+2)*(L+3)*(L+4)/6,0:(L+1)), PureSp(1:3*no,1:nLpure), Fund(1:no,0:L+1)
  real*8                        :: beta(1:no), BAx(1:no), BAy(1:no), BAz(1:no), nc, T(1:no), h(1:no,0:kmax,0:(L+1))
  integer                       :: n(1:no), k, i, j, Lc, Lx, Ly, Lz,m
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This program calculates pure real multipole moments of the 
  !! Sp shell pair up to L-poles.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  mpRgSp = 0d0;
  beta =RgSp(8,1:no) ;  BAx = RgSp(4,1:no) ;   BAy = RgSp(5,1:no) ;   BAz = RgSp(6,1:no) ;
  T = beta * (BAx*BAx + BAy*BAy + BAz*BAz);  n = int(RgSp(7,1:no))

  call calcmphk(h,kmax,L+1,T,no)                                   !  Compute h_k(T) and their derivatives
  call calcintmpSs(Fund,h,kmax,L+1,n,beta,no)                      !  Compute [S|s]^(m) integrals
  call calcmpCartSs(Rs,Fund,(L+2)*(L+3)*(L+4)/6,L+1,beta,BAx,BAy,BAz,no)        !  Compute [R|s] integrals
  call calcCartmpSp(CartSp,Rs,L,nLCart,BAx,BAy,BAz,no)
  call Cart2PureSp(PureSp, CartSp,no,nLCart,nLpure, L)
  
  mpRgSp(1:3*no,1:nLPure) = PureSp
  
  do i=1,no
     mpRgSp(3*i-2,nLPure+1) = RgSp(1,i) 
     mpRgSp(3*i-1,nLPure+1) = RgSp(1,i)  
     mpRgSp(3*i  ,nLPure+1) = RgSp(1,i)   
     mpRgSp(3*i-2,nLPure+2) = RgSp(2,i)  
     mpRgSp(3*i-1,nLPure+2) = RgSp(2,i)+1  
     mpRgSp(3*i  ,nLPure+2) = RgSp(2,i)+2
     mpRgSp(3*i-2,nLPure+3) = RgSp(9,i)
     mpRgSp(3*i-1,nLPure+3) = RgSp(9,i)
     mpRgSp(3*i  ,nLPure+3) = RgSp(9,i)
     mpRgSp(3*i-2,nLPure+4) = RgSp(3,i) 
     mpRgSp(3*i-1,nLPure+4) = RgSp(3,i)  
     mpRgSp(3*i  ,nLPure+4) = RgSp(3,i) 
  end do

 do i=1,no
 if ((BAx(i)**2+BAy(i)**2+BAz(i)**2) .lt. 1d-8 .and. beta(i) > 10d0) then
   call calcasympRgSp(mpRgSp(3*i-2:3*i,1:nLpure+4),1,RgSp(1:30,i), nLpure) 
 end if
 end do;

return

end subroutine calcgenmpRgSp

!**************************************************************************************************

subroutine calcCartmpSp(CartSp,Rs,L,nLCart,BAx,BAy,BAz,no)
  implicit none
  integer, intent(in)    :: nLCart, L, no
  real*8, intent(in)     :: BAx(1:no), BAy(1:no), BAz(1:no), Rs(1:no,1:(L+2)*(L+3)*(L+4)/6,0:(L+1))
  real*8, intent(out)    :: CartSp(1:3*no,1:nLCart)
  integer                :: i
  
  CartSp = 0d0;
  
  do i=1,no
     
     ! Charges of Sp basis function pairs
     CartSp(3*i-2,1) = -BAx(i)*Rs(i,1,0)+ Rs(i,2,0)              ! [S|px]
     CartSp(3*i-1,1) = -BAy(i)*Rs(i,1,0)+ Rs(i,3,0)              ! [S|py]
     CartSp(3*i,1)   = -BAz(i)*Rs(i,1,0)+ Rs(i,4,0)              ! [S|pz]
     
     
     ! dipoleSps of Sp basis function pairs
     if (L .ge. 1) then;
     CartSp(3*i-2,2) = -BAx(i)*Rs(i,2,0)+ Rs(i,5,0)             
     CartSp(3*i-1,2) = -BAy(i)*Rs(i,2,0)+ Rs(i,6,0)             
     CartSp(3*i,2)   = -BAz(i)*Rs(i,2,0)+ Rs(i,8,0)                
     
     CartSp(3*i-2,3) = -BAx(i)*Rs(i,3,0)+ Rs(i,6,0)             
     CartSp(3*i-1,3) = -BAy(i)*Rs(i,3,0)+ Rs(i,7,0)             
     CartSp(3*i,3)   = -BAz(i)*Rs(i,3,0)+ Rs(i,9,0)                
     
     CartSp(3*i-2,4) = -BAx(i)*Rs(i,4,0)+ Rs(i,8,0)             
     CartSp(3*i-1,4) = -BAy(i)*Rs(i,4,0)+ Rs(i,9,0)             
     CartSp(3*i,4)   = -BAz(i)*Rs(i,4,0)+ Rs(i,10,0)                
     end if
     
     ! Quadrupoles of Sp basis function pairs
     
     if (L .ge. 2) then;
     CartSp(3*i-2,5) = -BAx(i)*Rs(i,5,0)+ Rs(i,11,0)             
     CartSp(3*i-1,5) = -BAy(i)*Rs(i,5,0)+ Rs(i,12,0)             
     CartSp(3*i,5)   = -BAz(i)*Rs(i,5,0)+ Rs(i,15,0)                
     
     CartSp(3*i-2,6) = -BAx(i)*Rs(i,6,0)+ Rs(i,12,0)             
     CartSp(3*i-1,6) = -BAy(i)*Rs(i,6,0)+ Rs(i,13,0)             
     CartSp(3*i,6)   = -BAz(i)*Rs(i,6,0)+ Rs(i,16,0)                
     
     CartSp(3*i-2,7) = -BAx(i)*Rs(i,7,0)+ Rs(i,13,0)             
     CartSp(3*i-1,7) = -BAy(i)*Rs(i,7,0)+ Rs(i,14,0)             
     CartSp(3*i  ,7) = -BAz(i)*Rs(i,7,0)+ Rs(i,17,0)                
     
     CartSp(3*i-2,8) = -BAx(i)*Rs(i,8,0)+ Rs(i,15,0)             
     CartSp(3*i-1,8) = -BAy(i)*Rs(i,8,0)+ Rs(i,16,0)             
     CartSp(3*i,8)   = -BAz(i)*Rs(i,8,0)+ Rs(i,18,0)  
     
     CartSp(3*i-2,9) = -BAx(i)*Rs(i,9,0)+ Rs(i,16,0)             
     CartSp(3*i-1,9) = -BAy(i)*Rs(i,9,0)+ Rs(i,17,0)             
     CartSp(3*i,9)   = -BAz(i)*Rs(i,9,0)+ Rs(i,19,0)                
     
     CartSp(3*i-2,10) = -BAx(i)*Rs(i,10,0)+ Rs(i,18,0)             
     CartSp(3*i-1,10) = -BAy(i)*Rs(i,10,0)+ Rs(i,19,0)             
     CartSp(3*i,10)   = -BAz(i)*Rs(i,10,0)+ Rs(i,20,0)                
     end if
     
     
     ! Octupoles of Sp basis function pairs

  if (L .ge. 3) then;
  CartSp(3*i-2,11) = -BAx(i)*Rs(i,11,0)+ Rs(i,21,0)             
  CartSp(3*i-1,11) = -BAy(i)*Rs(i,11,0)+ Rs(i,22,0)             
  CartSp(3*i,11)   = -BAz(i)*Rs(i,11,0)+ Rs(i,26,0) 

  CartSp(3*i-2,12) = -BAx(i)*Rs(i,12,0)+ Rs(i,22,0)             
  CartSp(3*i-1,12) = -BAy(i)*Rs(i,12,0)+ Rs(i,23,0)             
  CartSp(3*i,12)   = -BAz(i)*Rs(i,12,0)+ Rs(i,27,0)                
              
  CartSp(3*i-2,13) = -BAx(i)*Rs(i,13,0)+ Rs(i,23,0)             
  CartSp(3*i-1,13) = -BAy(i)*Rs(i,13,0)+ Rs(i,24,0)             
  CartSp(3*i,13)   = -BAz(i)*Rs(i,13,0)+ Rs(i,28,0) 

  CartSp(3*i-2,14) = -BAx(i)*Rs(i,14,0)+ Rs(i,24,0)             
  CartSp(3*i-1,14) = -BAy(i)*Rs(i,14,0)+ Rs(i,25,0)             
  CartSp(3*i,14)   = -BAz(i)*Rs(i,14,0)+ Rs(i,29,0) 

  CartSp(3*i-2,15) = -BAx(i)*Rs(i,15,0)+ Rs(i,26,0)             
  CartSp(3*i-1,15) = -BAy(i)*Rs(i,15,0)+ Rs(i,27,0)             
  CartSp(3*i,15)   = -BAz(i)*Rs(i,15,0)+ Rs(i,30,0) 

  CartSp(3*i-2,16) = -BAx(i)*Rs(i,16,0)+ Rs(i,27,0)             
  CartSp(3*i-1,16) = -BAy(i)*Rs(i,16,0)+ Rs(i,28,0)             
  CartSp(3*i,16)   = -BAz(i)*Rs(i,16,0)+ Rs(i,31,0) 

  CartSp(3*i-2,17) = -BAx(i)*Rs(i,17,0)+ Rs(i,28,0)             
  CartSp(3*i-1,17) = -BAy(i)*Rs(i,17,0)+ Rs(i,29,0)             
  CartSp(3*i,17)   = -BAz(i)*Rs(i,17,0)+ Rs(i,32,0) 

  CartSp(3*i-2,18) = -BAx(i)*Rs(i,18,0)+ Rs(i,30,0)             
  CartSp(3*i-1,18) = -BAy(i)*Rs(i,18,0)+ Rs(i,31,0)             
  CartSp(3*i,18)   = -BAz(i)*Rs(i,18,0)+ Rs(i,33,0) 

  CartSp(3*i-2,19) = -BAx(i)*Rs(i,19,0)+ Rs(i,31,0)             
  CartSp(3*i-1,19) = -BAy(i)*Rs(i,19,0)+ Rs(i,32,0)             
  CartSp(3*i,19)   = -BAz(i)*Rs(i,19,0)+ Rs(i,34,0)

  CartSp(3*i-2,20) = -BAx(i)*Rs(i,20,0)+ Rs(i,33,0)             
  CartSp(3*i-1,20) = -BAy(i)*Rs(i,20,0)+ Rs(i,34,0)             
  CartSp(3*i,20)   = -BAz(i)*Rs(i,20,0)+ Rs(i,35,0) 
  end if

! ********** L = 4 ***********

  if (L .ge. 4) then;
  CartSp(3*i-2,21) = -BAx(i)*Rs(i,21,0)+ Rs(i,36,0)             
  CartSp(3*i-1,21) = -BAy(i)*Rs(i,21,0)+ Rs(i,37,0)             
  CartSp(3*i,21)   = -BAz(i)*Rs(i,21,0)+ Rs(i,42,0) 

  CartSp(3*i-2,22) = -BAx(i)*Rs(i,22,0)+ Rs(i,37,0)             
  CartSp(3*i-1,22) = -BAy(i)*Rs(i,22,0)+ Rs(i,38,0)             
  CartSp(3*i  ,22) = -BAz(i)*Rs(i,22,0)+ Rs(i,43,0) 

  CartSp(3*i-2,23) = -BAx(i)*Rs(i,23,0)+ Rs(i,38,0)             
  CartSp(3*i-1,23) = -BAy(i)*Rs(i,23,0)+ Rs(i,39,0)             
  CartSp(3*i  ,23) = -BAz(i)*Rs(i,23,0)+ Rs(i,44,0) 

  CartSp(3*i-2,24) = -BAx(i)*Rs(i,24,0)+ Rs(i,39,0)             
  CartSp(3*i-1,24) = -BAy(i)*Rs(i,24,0)+ Rs(i,40,0)             
  CartSp(3*i  ,24) = -BAz(i)*Rs(i,24,0)+ Rs(i,45,0) 

  CartSp(3*i-2,25) = -BAx(i)*Rs(i,25,0)+ Rs(i,40,0)             
  CartSp(3*i-1,25) = -BAy(i)*Rs(i,25,0)+ Rs(i,41,0)             
  CartSp(3*i  ,25) = -BAz(i)*Rs(i,25,0)+ Rs(i,46,0) 

  CartSp(3*i-2,26) = -BAx(i)*Rs(i,26,0)+ Rs(i,42,0)             
  CartSp(3*i-1,26) = -BAy(i)*Rs(i,26,0)+ Rs(i,43,0)             
  CartSp(3*i  ,26) = -BAz(i)*Rs(i,26,0)+ Rs(i,47,0) 

  CartSp(3*i-2,27) = -BAx(i)*Rs(i,27,0)+ Rs(i,43,0)             
  CartSp(3*i-1,27) = -BAy(i)*Rs(i,27,0)+ Rs(i,44,0)             
  CartSp(3*i  ,27) = -BAz(i)*Rs(i,27,0)+ Rs(i,48,0) 

  CartSp(3*i-2,28) = -BAx(i)*Rs(i,28,0)+ Rs(i,44,0)             
  CartSp(3*i-1,28) = -BAy(i)*Rs(i,28,0)+ Rs(i,45,0)             
  CartSp(3*i  ,28) = -BAz(i)*Rs(i,28,0)+ Rs(i,49,0) 

  CartSp(3*i-2,29) = -BAx(i)*Rs(i,29,0)+ Rs(i,45,0)             
  CartSp(3*i-1,29) = -BAy(i)*Rs(i,29,0)+ Rs(i,46,0)             
  CartSp(3*i  ,29) = -BAz(i)*Rs(i,29,0)+ Rs(i,50,0) 

  CartSp(3*i-2,30) = -BAx(i)*Rs(i,30,0)+ Rs(i,47,0)             
  CartSp(3*i-1,30) = -BAy(i)*Rs(i,30,0)+ Rs(i,48,0)             
  CartSp(3*i  ,30) = -BAz(i)*Rs(i,30,0)+ Rs(i,51,0) 

  CartSp(3*i-2,31) = -BAx(i)*Rs(i,31,0)+ Rs(i,48,0)             
  CartSp(3*i-1,31) = -BAy(i)*Rs(i,31,0)+ Rs(i,49,0)             
  CartSp(3*i  ,31) = -BAz(i)*Rs(i,31,0)+ Rs(i,52,0) 

  CartSp(3*i-2,32) = -BAx(i)*Rs(i,32,0)+ Rs(i,49,0)             
  CartSp(3*i-1,32) = -BAy(i)*Rs(i,32,0)+ Rs(i,50,0)             
  CartSp(3*i  ,32) = -BAz(i)*Rs(i,32,0)+ Rs(i,53,0) 

  CartSp(3*i-2,33) = -BAx(i)*Rs(i,33,0)+ Rs(i,51,0)             
  CartSp(3*i-1,33) = -BAy(i)*Rs(i,33,0)+ Rs(i,52,0)             
  CartSp(3*i  ,33) = -BAz(i)*Rs(i,33,0)+ Rs(i,54,0) 

  CartSp(3*i-2,34) = -BAx(i)*Rs(i,34,0)+ Rs(i,52,0)             
  CartSp(3*i-1,34) = -BAy(i)*Rs(i,34,0)+ Rs(i,53,0)             
  CartSp(3*i  ,34) = -BAz(i)*Rs(i,34,0)+ Rs(i,55,0) 

  CartSp(3*i-2,35) = -BAx(i)*Rs(i,35,0)+ Rs(i,54,0)             
  CartSp(3*i-1,35) = -BAy(i)*Rs(i,35,0)+ Rs(i,55,0)             
  CartSp(3*i  ,35) = -BAz(i)*Rs(i,35,0)+ Rs(i,56,0) 
  end if


! ********** L = 5 ***********

  if (L .ge. 5) then;
  CartSp(3*i-2,36) = -BAx(i)*Rs(i,36,0)+ Rs(i,57,0)             
  CartSp(3*i-1,36) = -BAy(i)*Rs(i,36,0)+ Rs(i,58,0)             
  CartSp(3*i  ,36) = -BAz(i)*Rs(i,36,0)+ Rs(i,64,0) 

  CartSp(3*i-2,37) = -BAx(i)*Rs(i,37,0)+ Rs(i,58,0)             
  CartSp(3*i-1,37) = -BAy(i)*Rs(i,37,0)+ Rs(i,59,0)             
  CartSp(3*i  ,37) = -BAz(i)*Rs(i,37,0)+ Rs(i,65,0)

  CartSp(3*i-2,38) = -BAx(i)*Rs(i,38,0)+ Rs(i,59,0)             
  CartSp(3*i-1,38) = -BAy(i)*Rs(i,38,0)+ Rs(i,60,0)             
  CartSp(3*i  ,38) = -BAz(i)*Rs(i,38,0)+ Rs(i,66,0)

  CartSp(3*i-2,39) = -BAx(i)*Rs(i,39,0)+ Rs(i,60,0)             
  CartSp(3*i-1,39) = -BAy(i)*Rs(i,39,0)+ Rs(i,61,0)             
  CartSp(3*i  ,39) = -BAz(i)*Rs(i,39,0)+ Rs(i,67,0)

  CartSp(3*i-2,40) = -BAx(i)*Rs(i,40,0)+ Rs(i,61,0)             
  CartSp(3*i-1,40) = -BAy(i)*Rs(i,40,0)+ Rs(i,62,0)             
  CartSp(3*i  ,40) = -BAz(i)*Rs(i,40,0)+ Rs(i,68,0)

  CartSp(3*i-2,41) = -BAx(i)*Rs(i,41,0)+ Rs(i,62,0)             
  CartSp(3*i-1,41) = -BAy(i)*Rs(i,41,0)+ Rs(i,63,0)             
  CartSp(3*i  ,41) = -BAz(i)*Rs(i,41,0)+ Rs(i,69,0)

  CartSp(3*i-2,42) = -BAx(i)*Rs(i,42,0)+ Rs(i,64,0)             
  CartSp(3*i-1,42) = -BAy(i)*Rs(i,42,0)+ Rs(i,65,0)             
  CartSp(3*i  ,42) = -BAz(i)*Rs(i,42,0)+ Rs(i,70,0)

  CartSp(3*i-2,43) = -BAx(i)*Rs(i,43,0)+ Rs(i,65,0)             
  CartSp(3*i-1,43) = -BAy(i)*Rs(i,43,0)+ Rs(i,66,0)             
  CartSp(3*i  ,43) = -BAz(i)*Rs(i,43,0)+ Rs(i,71,0)

  CartSp(3*i-2,44) = -BAx(i)*Rs(i,44,0)+ Rs(i,66,0)             
  CartSp(3*i-1,44) = -BAy(i)*Rs(i,44,0)+ Rs(i,67,0)             
  CartSp(3*i  ,44) = -BAz(i)*Rs(i,44,0)+ Rs(i,72,0)

  CartSp(3*i-2,45) = -BAx(i)*Rs(i,45,0)+ Rs(i,67,0)             
  CartSp(3*i-1,45) = -BAy(i)*Rs(i,45,0)+ Rs(i,68,0)             
  CartSp(3*i  ,45) = -BAz(i)*Rs(i,45,0)+ Rs(i,73,0)

  CartSp(3*i-2,46) = -BAx(i)*Rs(i,46,0)+ Rs(i,68,0)             
  CartSp(3*i-1,46) = -BAy(i)*Rs(i,46,0)+ Rs(i,69,0)             
  CartSp(3*i  ,46) = -BAz(i)*Rs(i,46,0)+ Rs(i,74,0)

  CartSp(3*i-2,47) = -BAx(i)*Rs(i,47,0)+ Rs(i,70,0)             
  CartSp(3*i-1,47) = -BAy(i)*Rs(i,47,0)+ Rs(i,71,0)             
  CartSp(3*i  ,47) = -BAz(i)*Rs(i,47,0)+ Rs(i,75,0)

  CartSp(3*i-2,48) = -BAx(i)*Rs(i,48,0)+ Rs(i,71,0)             
  CartSp(3*i-1,48) = -BAy(i)*Rs(i,48,0)+ Rs(i,72,0)             
  CartSp(3*i  ,48) = -BAz(i)*Rs(i,48,0)+ Rs(i,76,0)

  CartSp(3*i-2,49) = -BAx(i)*Rs(i,49,0)+ Rs(i,72,0)             
  CartSp(3*i-1,49) = -BAy(i)*Rs(i,49,0)+ Rs(i,73,0)             
  CartSp(3*i  ,49) = -BAz(i)*Rs(i,49,0)+ Rs(i,77,0)

  CartSp(3*i-2,50) = -BAx(i)*Rs(i,50,0)+ Rs(i,73,0)             
  CartSp(3*i-1,50) = -BAy(i)*Rs(i,50,0)+ Rs(i,74,0)             
  CartSp(3*i  ,50) = -BAz(i)*Rs(i,50,0)+ Rs(i,78,0)

  CartSp(3*i-2,51) = -BAx(i)*Rs(i,51,0)+ Rs(i,75,0)             
  CartSp(3*i-1,51) = -BAy(i)*Rs(i,51,0)+ Rs(i,76,0)             
  CartSp(3*i  ,51) = -BAz(i)*Rs(i,51,0)+ Rs(i,79,0)

  CartSp(3*i-2,52) = -BAx(i)*Rs(i,52,0)+ Rs(i,76,0)             
  CartSp(3*i-1,52) = -BAy(i)*Rs(i,52,0)+ Rs(i,77,0)             
  CartSp(3*i  ,52) = -BAz(i)*Rs(i,52,0)+ Rs(i,80,0)

  CartSp(3*i-2,53) = -BAx(i)*Rs(i,53,0)+ Rs(i,77,0)             
  CartSp(3*i-1,53) = -BAy(i)*Rs(i,53,0)+ Rs(i,78,0)             
  CartSp(3*i  ,53) = -BAz(i)*Rs(i,53,0)+ Rs(i,81,0)

  CartSp(3*i-2,54) = -BAx(i)*Rs(i,54,0)+ Rs(i,79,0)             
  CartSp(3*i-1,54) = -BAy(i)*Rs(i,54,0)+ Rs(i,80,0)             
  CartSp(3*i  ,54) = -BAz(i)*Rs(i,54,0)+ Rs(i,82,0)

  CartSp(3*i-2,55) = -BAx(i)*Rs(i,55,0)+ Rs(i,80,0)             
  CartSp(3*i-1,55) = -BAy(i)*Rs(i,55,0)+ Rs(i,81,0)             
  CartSp(3*i  ,55) = -BAz(i)*Rs(i,55,0)+ Rs(i,83,0)

  CartSp(3*i-2,56) = -BAx(i)*Rs(i,56,0)+ Rs(i,82,0)             
  CartSp(3*i-1,56) = -BAy(i)*Rs(i,56,0)+ Rs(i,83,0)             
  CartSp(3*i  ,56) = -BAz(i)*Rs(i,56,0)+ Rs(i,84,0)
  end if


! ********** L = 6 ***********



  if (L .ge. 6) then;
  CartSp(3*i-2,57) = -BAx(i)*Rs(i,57,0)+ Rs(i,85,0)             
  CartSp(3*i-1,57) = -BAy(i)*Rs(i,57,0)+ Rs(i,86,0)             
  CartSp(3*i  ,57) = -BAz(i)*Rs(i,57,0)+ Rs(i,93,0)

  CartSp(3*i-2,58) = -BAx(i)*Rs(i,58,0)+ Rs(i,86,0)             
  CartSp(3*i-1,58) = -BAy(i)*Rs(i,58,0)+ Rs(i,87,0)             
  CartSp(3*i  ,58) = -BAz(i)*Rs(i,58,0)+ Rs(i,94,0)

  CartSp(3*i-2,59) = -BAx(i)*Rs(i,59,0)+ Rs(i,87,0)             
  CartSp(3*i-1,59) = -BAy(i)*Rs(i,59,0)+ Rs(i,88,0)             
  CartSp(3*i  ,59) = -BAz(i)*Rs(i,59,0)+ Rs(i,95,0)

  CartSp(3*i-2,60) = -BAx(i)*Rs(i,60,0)+ Rs(i,88,0)             
  CartSp(3*i-1,60) = -BAy(i)*Rs(i,60,0)+ Rs(i,89,0)             
  CartSp(3*i  ,60) = -BAz(i)*Rs(i,60,0)+ Rs(i,96,0)

  CartSp(3*i-2,61) = -BAx(i)*Rs(i,61,0)+ Rs(i,89,0)             
  CartSp(3*i-1,61) = -BAy(i)*Rs(i,61,0)+ Rs(i,90,0)             
  CartSp(3*i  ,61) = -BAz(i)*Rs(i,61,0)+ Rs(i,97,0)

  CartSp(3*i-2,62) = -BAx(i)*Rs(i,62,0)+ Rs(i,90,0)             
  CartSp(3*i-1,62) = -BAy(i)*Rs(i,62,0)+ Rs(i,91,0)             
  CartSp(3*i  ,62) = -BAz(i)*Rs(i,62,0)+ Rs(i,98,0)

  CartSp(3*i-2,63) = -BAx(i)*Rs(i,63,0)+ Rs(i,91,0)             
  CartSp(3*i-1,63) = -BAy(i)*Rs(i,63,0)+ Rs(i,92,0)             
  CartSp(3*i  ,63) = -BAz(i)*Rs(i,63,0)+ Rs(i,99,0)

  CartSp(3*i-2,64) = -BAx(i)*Rs(i,64,0)+ Rs(i,93,0)             
  CartSp(3*i-1,64) = -BAy(i)*Rs(i,64,0)+ Rs(i,94,0)             
  CartSp(3*i  ,64) = -BAz(i)*Rs(i,64,0)+ Rs(i,100,0)

  CartSp(3*i-2,65) = -BAx(i)*Rs(i,65,0)+ Rs(i,94,0)             
  CartSp(3*i-1,65) = -BAy(i)*Rs(i,65,0)+ Rs(i,95,0)             
  CartSp(3*i  ,65) = -BAz(i)*Rs(i,65,0)+ Rs(i,101,0)

  CartSp(3*i-2,66) = -BAx(i)*Rs(i,66,0)+ Rs(i,95,0)             
  CartSp(3*i-1,66) = -BAy(i)*Rs(i,66,0)+ Rs(i,96,0)             
  CartSp(3*i  ,66) = -BAz(i)*Rs(i,66,0)+ Rs(i,102,0)

  CartSp(3*i-2,67) = -BAx(i)*Rs(i,67,0)+ Rs(i,96,0)             
  CartSp(3*i-1,67) = -BAy(i)*Rs(i,67,0)+ Rs(i,97,0)             
  CartSp(3*i  ,67) = -BAz(i)*Rs(i,67,0)+ Rs(i,103,0)

  CartSp(3*i-2,68) = -BAx(i)*Rs(i,68,0)+ Rs(i,97,0)             
  CartSp(3*i-1,68) = -BAy(i)*Rs(i,68,0)+ Rs(i,98,0)             
  CartSp(3*i  ,68) = -BAz(i)*Rs(i,68,0)+ Rs(i,104,0)

  CartSp(3*i-2,69) = -BAx(i)*Rs(i,69,0)+ Rs(i,98,0)             
  CartSp(3*i-1,69) = -BAy(i)*Rs(i,69,0)+ Rs(i,99,0)             
  CartSp(3*i  ,69) = -BAz(i)*Rs(i,69,0)+ Rs(i,105,0)

  CartSp(3*i-2,70) = -BAx(i)*Rs(i,70,0)+ Rs(i,100,0)             
  CartSp(3*i-1,70) = -BAy(i)*Rs(i,70,0)+ Rs(i,101,0)             
  CartSp(3*i  ,70) = -BAz(i)*Rs(i,70,0)+ Rs(i,106,0)

  CartSp(3*i-2,71) = -BAx(i)*Rs(i,71,0)+ Rs(i,101,0)             
  CartSp(3*i-1,71) = -BAy(i)*Rs(i,71,0)+ Rs(i,102,0)             
  CartSp(3*i  ,71) = -BAz(i)*Rs(i,71,0)+ Rs(i,107,0)

  CartSp(3*i-2,72) = -BAx(i)*Rs(i,72,0)+ Rs(i,102,0)             
  CartSp(3*i-1,72) = -BAy(i)*Rs(i,72,0)+ Rs(i,103,0)             
  CartSp(3*i  ,72) = -BAz(i)*Rs(i,72,0)+ Rs(i,108,0)

  CartSp(3*i-2,73) = -BAx(i)*Rs(i,73,0)+ Rs(i,103,0)             
  CartSp(3*i-1,73) = -BAy(i)*Rs(i,73,0)+ Rs(i,104,0)             
  CartSp(3*i  ,73) = -BAz(i)*Rs(i,73,0)+ Rs(i,109,0)

  CartSp(3*i-2,74) = -BAx(i)*Rs(i,74,0)+ Rs(i,104,0)             
  CartSp(3*i-1,74) = -BAy(i)*Rs(i,74,0)+ Rs(i,105,0)             
  CartSp(3*i  ,74) = -BAz(i)*Rs(i,74,0)+ Rs(i,110,0)

  CartSp(3*i-2,75) = -BAx(i)*Rs(i,75,0)+ Rs(i,106,0)             
  CartSp(3*i-1,75) = -BAy(i)*Rs(i,75,0)+ Rs(i,107,0)             
  CartSp(3*i  ,75) = -BAz(i)*Rs(i,75,0)+ Rs(i,111,0)

  CartSp(3*i-2,76) = -BAx(i)*Rs(i,76,0)+ Rs(i,107,0)             
  CartSp(3*i-1,76) = -BAy(i)*Rs(i,76,0)+ Rs(i,108,0)             
  CartSp(3*i  ,76) = -BAz(i)*Rs(i,76,0)+ Rs(i,112,0)

  CartSp(3*i-2,77) = -BAx(i)*Rs(i,77,0)+ Rs(i,108,0)             
  CartSp(3*i-1,77) = -BAy(i)*Rs(i,77,0)+ Rs(i,109,0)             
  CartSp(3*i  ,77) = -BAz(i)*Rs(i,77,0)+ Rs(i,113,0)

  CartSp(3*i-2,78) = -BAx(i)*Rs(i,78,0)+ Rs(i,109,0)             
  CartSp(3*i-1,78) = -BAy(i)*Rs(i,78,0)+ Rs(i,110,0)             
  CartSp(3*i  ,78) = -BAz(i)*Rs(i,78,0)+ Rs(i,114,0)

  CartSp(3*i-2,79) = -BAx(i)*Rs(i,79,0)+ Rs(i,111,0)             
  CartSp(3*i-1,79) = -BAy(i)*Rs(i,79,0)+ Rs(i,112,0)             
  CartSp(3*i  ,79) = -BAz(i)*Rs(i,79,0)+ Rs(i,115,0)

  CartSp(3*i-2,80) = -BAx(i)*Rs(i,80,0)+ Rs(i,112,0)             
  CartSp(3*i-1,80) = -BAy(i)*Rs(i,80,0)+ Rs(i,113,0)             
  CartSp(3*i  ,80) = -BAz(i)*Rs(i,80,0)+ Rs(i,116,0)

  CartSp(3*i-2,81) = -BAx(i)*Rs(i,81,0)+ Rs(i,113,0)             
  CartSp(3*i-1,81) = -BAy(i)*Rs(i,81,0)+ Rs(i,114,0)             
  CartSp(3*i  ,81) = -BAz(i)*Rs(i,81,0)+ Rs(i,117,0)

  CartSp(3*i-2,82) = -BAx(i)*Rs(i,82,0)+ Rs(i,115,0)             
  CartSp(3*i-1,82) = -BAy(i)*Rs(i,82,0)+ Rs(i,116,0)             
  CartSp(3*i  ,82) = -BAz(i)*Rs(i,82,0)+ Rs(i,118,0)

  CartSp(3*i-2,83) = -BAx(i)*Rs(i,83,0)+ Rs(i,116,0)             
  CartSp(3*i-1,83) = -BAy(i)*Rs(i,83,0)+ Rs(i,117,0)             
  CartSp(3*i  ,83) = -BAz(i)*Rs(i,83,0)+ Rs(i,119,0)

  CartSp(3*i-2,84) = -BAx(i)*Rs(i,84,0)+ Rs(i,118,0)             
  CartSp(3*i-1,84) = -BAy(i)*Rs(i,84,0)+ Rs(i,119,0)             
  CartSp(3*i  ,84) = -BAz(i)*Rs(i,84,0)+ Rs(i,120,0)
  end if
 end do

return
end subroutine calcCartmpSp

!****************************************************

subroutine Cart2PureSp(PureSp, CartSp,no,nLCart,nLpure, L)
   implicit none
   integer, intent(in)   :: L, no, nLCart, nLPure
   real*8, intent(in)    :: CartSp(1:3*no,1:nLCart)
   real*8, intent(out)   :: PureSp(1:3*no,1:nLPure)
   real*8                :: t, invfact(0:2*L)
   integer               :: fact(0:2*L), lindex, mindex, cnt, i, k
   real*8  :: CartSpx(1:no,1:nLCart,0:L), CartSpy(1:no,1:nLCart,0:L), CartSpz(1:no,1:nLCart,0:L)
   real*8  :: PureSpx(1:no,1:nLCart    ), PureSpy(1:no,1:nLPure    ), PureSpz(1:no,1:nLPure    )

   PureSp = 0d0;  CartSpx = 0d0;  CartSpy = 0d0;  CartSpz = 0d0;
do i=1,no
   CartSpx(i,1:nLCart,0) = CartSp(3*i-2,1:nLCart)
   CartSpy(i,1:nLCart,0) = CartSp(3*i-1,1:nLCart)
   CartSpz(i,1:nLCart,0) = CartSp(3*i-0,1:nLCart)
end do

   call Cart2Pure(PureSpx,CartSpx,no,nLCart,nLpure,L)
   call Cart2Pure(PureSpy,CartSpy,no,nLCart,nLpure,L)
   call Cart2Pure(PureSpz,CartSpz,no,nLCart,nLpure,L)

do i=1,no
   PureSp(3*i-2,1:nLPure) = PureSpx(i,1:nLPure)
   PureSp(3*i-1,1:nLPure) = PureSpy(i,1:nLPure)
   PureSp(3*i  ,1:nLPure) = PureSpz(i,1:nLPure)
end do

return
end subroutine Cart2PureSp

!*********************************************************************************


!****************************************

subroutine asyRgSpdip(dipoleSp,no,RgSp)
 implicit none
 integer, intent(in) :: no
 real*8, intent(out) :: dipoleSp(1:no)
 real*8, intent(in)  :: RgSp(1:30,1:no)
 real*8              :: beta, b
 integer             :: n, i

do i=1,no
 n    = int(RgSp(7,i));
 beta = RgSp(8,i)
 b    = 1d0/sqrt(beta);
 select case (n)
   case (3)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.2566370614359172d1+(2.0881229988118903d1-  &
1.2566370614359172d1*b)*b))
   case (4)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.6755160819145565d1+b*(4.176245997623781d1+  &
b*(-5.026548245743669d1+2.4361434986138724d1*b))))
   case (5)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-2.0943951023931957d1+b*(6.960409996039635d1+  &
b*(-1.2566370614359172d2+(1.2180717493069362d2-5.026548245743669d1*b)*b))))
   case (6)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-2.5132741228718345d1+b*(1.0440614994059452d2+  &
b*(-2.5132741228718345d2+b*(3.6542152479208085d2+b*(-3.015928947446201d2+  &
1.0962645743762425d2*b))))))
   case (7)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-2.9321531433504737d1+b*(1.4616860991683234d2+  &
b*(-4.39822971502571d2+b*(8.526502245148553d2+b*(-1.0555751316061706d3+  &
(7.673852020633697d2-2.5132741228718345d2*b)*b))))))
   case (8)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-3.351032163829113d1+b*(1.9489147988910975d2+  &
b*(-7.037167544041136d2+b*(1.7053004490297106d3+b*(-2.8148670176164545d3+  &
b*(3.069540808253479d3+b*(-2.0106192982974678d3+6.029455159069333d2*b))))))))
   case (9)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-3.769911184307752d1+b*(2.5057475985742683d2+  &
b*(-1.0555751316061706d3+b*(3.069540808253479d3+b*(-6.333450789637023d3+  &
b*(9.208622424760437d3+b*(-9.047786842338605d3+(5.426509643162401d3-  &
1.5079644737231006d3*b)*b))))))))
   case (10)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-4.188790204786391d1+b*(3.1321844982178355d2+  &
b*(-1.5079644737231006d3+b*(5.115901347089132d3+b*(-1.2666901579274046d4+  &
b*(2.302155606190109d4+b*(-3.015928947446201d4+b*(2.7132548215812005d4+  &
b*(-1.5079644737231006d4+3.919145853395067d3*b))))))))))
   case (11)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-4.60766922526503d1+b*(3.8282254978217987d2+  &
b*(-2.0734511513692637d3+b*(8.039273545425779d3+b*(-2.322265289533575d4+  &
b*(5.06474233361824d4+b*(-8.293804605477055d4+b*(9.948601012464401d4+  &
b*(-8.293804605477055d4+(4.3110604387345735d4-  &
1.0555751316061706d4*b)*b))))))))))
   case (12)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-5.026548245743669d1+b*(4.593870597386159d2+  &
b*(-2.7646015351590183d3+b*(1.2058910318138667d4+b*(-3.981026210628986d4+  &
b*(1.0129484667236481d5+b*(-1.990513105314493d5+b*(2.98458030373932d5+  &
b*(-3.3175218421908212d5+b*(2.5866362632407442d5+b*(-1.2666901579274046d5+  &
2.9393593900462998d4*b))))))))))))
   case (13)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-5.445427266222308d1+b*(5.429119796910915d2+  &
b*(-3.5939819957067236d3+b*(1.7418426015089186d4+b*(-6.469167592272102d4+  &
b*(1.8811900096296321d5+b*(-4.312778394848068d5+b*(7.759908789722232d5+  &
b*(-1.078194598712017d6+b*(1.120875714070989d6+b*(-8.23348602652813d5+  &
(3.8211672070601903d5-8.444601052849364d4*b)*b))))))))))))
   case (14)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-5.8643062867009474d1+b*(6.333973096396067d2+  &
b*(-4.5741589036267385d3+b*(2.438579642112486d4+b*(-1.0063149587978826d5+  &
b*(3.292082516851856d5+b*(-8.625556789696136d5+b*(1.8106453842685208d6+  &
b*(-3.0189448763936477d6+b*(3.923064999248462d6+b*(-3.842293479046461d6+  &
b*(2.674817044942133d6+b*(-1.182244147398911d6+  &
2.498455481539355d5*b))))))))))))))
   case (15)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-6.283185307179586d1+b*(7.308430495841617d2+  &
b*(-5.717698629533424d3+b*(3.3253358756079354d4+b*(-1.5094724381968239d5+  &
b*(5.486804194753093d5+b*(-1.6172918980680255d6+b*(3.879954394861116d6+  &
b*(-7.547362190984119d6+b*(1.1769194997745385d7+b*(-1.4408600546424228d7+  &
b*(1.3374085224710666d7+b*(-8.866831105491833d6+(3.747683222309033d6-  &
7.6001409475644275d5*b)*b))))))))))))))
   case (16)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-6.702064327658226d1+b*(8.352491995247561d2+  &
b*(-7.037167544041136d3+b*(4.433781167477248d4+b*(-2.1955962737408345d5+  &
b*(8.77888671160495d5+b*(-2.875185596565379d6+b*(7.759908789722232d6+  &
b*(-1.7251113579392272d7+b*(3.13845199939877d7+b*(-4.610752174855753d7+  &
b*(5.349634089884266d7+b*(-4.728976589595644d7+b*(2.9981465778472263d7+  &
b*(-1.2160225516103085d7+2.3735327074623873d6*b))))))))))))))))
   case (17)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-7.120943348136865d1+b*(9.466157594613904d2+  &
b*(-8.545132017764237d3+b*(5.798021526701016d4+b*(-3.110428054466183d5+  &
b*(1.3567370372480378d6+b*(-4.887815514161144d6+b*(1.465760549169755d7+  &
b*(-3.6658616356208578d7+b*(7.621954855682727d7+b*(-1.3063797828757968d8+  &
b*(1.8188755905606504d8+b*(-2.0098150505781485d8+b*(1.6989497274467615d8+  &
b*(-1.0336191688687622d8+(4.035005602686058d7-  &
7.6001409475644275d6*b)*b))))))))))))))))
   case (18)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-7.539822368615504d1+b*(1.064942729394064d3+  &
b*(-1.0254158421317086d4+b*(7.454599105758449d4+b*(-4.306746536953176d5+  &
b*(2.0351055558720565d6+b*(-7.998243568627327d6+b*(2.638368988505559d7+  &
b*(-7.3317232712417155d7+b*(1.7149398425286133d8+b*(-3.3592622988234773d8+  &
b*(5.456626771681952d8+b*(-7.235334182081336d8+b*(7.645273773510427d8+  &
b*(-6.201715013212573d8+b*(3.631505042417453d8+b*(-1.3680253705615968d8+  &
2.492209342835507d7*b))))))))))))))))))
   case (19)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-7.958701389094143d1+b*(1.1902301093227776d3+  &
b*(-1.2176813125314039d4+b*(9.442492200627369d4+b*(-5.844870300150738d5+  &
b*(2.974385043197621d6+b*(-1.26638856503266d7+b*(4.5571827983277835d7+  &
b*(-1.393027421535926d8+b*(3.620428556449295d8+b*(-7.978247959705758d8+  &
b*(1.4810844094565299d9+b*(-2.2911891576590895d9+b*(2.9052040339339618d9+  &
b*(-2.945814631275972d9+b*(2.2999531935310533d9+b*(-1.2996241020335173d9+  &
(4.735197751387463d8-8.360155042320871d7*b)*b))))))))))))))))))
   case (20)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-8.377580409572783d1+b*(1.3224778992475306d3+  &
b*(-1.4325662500369458d4+b*(1.180311525078421d5+b*(-7.793160400200985d5+  &
b*(4.249121490282316d6+b*(-1.9482901000502462d7+b*(7.59530466387964d7+  &
b*(-2.53277713006532d8+b*(7.24085711289859d8+b*(-1.7729439910457239d9+  &
b*(3.7027110236413243d9+b*(-6.5462547361688275d9+b*(9.68401344644654d9+  &
b*(-1.1783258525103888d10+b*(1.1499765967655267d10+b*(-8.664160680223446d9+  &
b*(4.735197751387463d9+b*(-1.672031008464174d9+  &
2.8660407442608324d8*b))))))))))))))))))))
   case (21)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-8.79645943005142d1+b*(1.4616860991683234d3+  &
b*(-1.67132729170977d4+b*(1.4580318839204027d5+b*(-1.0228523025263794d6+  &
b*(5.948770086395242d6+b*(-2.922435150075369d7+b*(1.2269338303190187d8+  &
b*(-4.43235997761431d8+b*(1.3823454488260944d9+b*(-3.72318238119602d9+  &
b*(8.63965905516309d9+b*(-1.7183918682443173d10+b*(2.9052040339339618d10+  &
b*(-4.124140483786361d10+b*(4.829901706415212d10+b*(-4.548684357117311d10+  &
b*(3.314638425971224d10+b*(-1.755632558887383d10+(6.0186855629477485d9-  &
1.0032186050785044d9*b)*b))))))))))))))))))))
   case (22)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-9.21533845053006d1+b*(1.6078547090851556d3+  &
b*(-1.9352210746113128d4+b*(1.7820389692360474d5+b*(-1.3236912150341378d6+  &
b*(8.179558868793459d6+b*(-4.286238220110541d7+b*(1.9280388762156009d8+  &
b*(-7.5009168851934485d8+b*(2.534299989514506d9+b*(-7.44636476239204d9+  &
b*(1.90072499213588d10+b*(-4.200513455708331d10+b*(7.989311093318396d10+  &
b*(-1.2961584377614277d11+b*(1.770963959018911d11+b*(-2.0014211171316165d11+  &
b*(1.823051134284173d11+b*(-1.287463876517414d11+b*(6.620554119242524d10+  &
b*(-2.20708093117271d10+3.582550930326041d9*b))))))))))))))))))))))
   case (23)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-9.6342174710087d1+b*(1.7609837289980277d3+  &
b*(-2.2255042358030095d4+b*(2.157205068022584d5+b*(-1.6913832192102873d6+  &
b*(1.1066461998955854d7+b*(-6.161467441408903d7+b*(2.9563262768639214d8+  &
b*(-1.2322934882817806d9+b*(4.4837615199102805d9+b*(-1.427219912791808d10+  &
b*(3.9742431653750216d10+b*(-9.66118094812916d10+b*(2.041712834959146d11+  &
b*(-3.726455508564105d11+b*(5.8188815796335644d11+b*(-7.672114282337863d11+  &
b*(8.386035217707196d11+b*(-7.4029172899751305d11+b*(5.075758158085934d11+  &
b*(-2.5381430708486166d11+(8.239867139749894d10-  &
1.3041841866020558d10*b)*b))))))))))))))))))))))
   case (24)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.0053096491487339d2+b*(1.921073158906939d3+  &
b*(-2.5434334123462965d4+b*(2.5886460816271004d5+b*(-2.136484066370889d6+  &
b*(1.4755282665274474d7+b*(-8.698542270224333d7+b*(4.434489415295882d8+  &
b*(-1.9716695812508491d9+b*(7.686448319846195d9+b*(-2.634867531307953d10+  &
b*(7.948486330750043d10+b*(-2.1078940250463623d11+b*(4.90011080390195d11+  &
b*(-9.93721468950428d11+b*(1.7456644738900695d12+b*(-2.6304391825158384d12+  &
b*(3.354414087082879d12+b*(-3.553400299188063d12+b*(3.0454548948515607d12+  &
b*(-2.030514456678893d12+b*(9.887840567699872d11+b*(-3.130042047844934d11+  &
4.836443755940155d10*b))))))))))))))))))))))))
   case (25)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.0471975511965979d2+b*(2.0881229988118903d3+  &
b*(-2.8902652413026098d4+b*(3.0817215257465485d5+b*(-2.670605082963611d6+  &
b*(1.9414845612203255d7+b*(-1.208130870864491d8+b*(6.521307963670414d8+  &
b*(-3.0807337207044516d9+b*(1.2810747199743657d10+b*(-4.705120591621345d10+  &
b*(1.5285550636057774d11+b*(-4.391445885513255d11+b*(1.1136615463413522d12+  &
b*(-2.48430367237607d12+b*(4.849067983027971d12+b*(-8.220122445361996d12+  &
b*(1.198005031101028d13+b*(-1.480583457995026d13+b*(1.5227274474257804d13+  &
b*(-1.2690715354243083d13+b*(8.239867139749894d12+b*(-3.9125525598061675d12+  &
(1.2091109389850387d12-1.825857861242878d11*b)*b))))))))))))))))))))))))
   case (26)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.0890854532444616d2+b*(2.2621332487128814d3+  &
b*(-3.267256359733385d4+b*(3.642034530427739d5+b*(-3.3064634360501852d6+  &
b*(2.5239299295864233d7+b*(-1.6532317180250926d8+b*(9.419667058635042d8+  &
b*(-4.711710396371514d9+b*(2.0817464199583444d10+b*(-8.15554235881033d10+  &
b*(2.838745118125015d11+b*(-8.78289177102651d11+b*(2.412933350406263d12+  &
b*(-5.871990498343438d12+b*(1.2607576755872723d13+b*(-2.3747020397712433d13+  &
b*(3.893516351078341d13+b*(-5.499309986838669d13+b*(6.598485605511715d13+  &
b*(-6.599171984206402d13+b*(5.35591364083743d13+b*(-3.390878885165345d13+  &
b*(1.5718442206805505d13+b*(-4.747230439231483d12+  &
7.012843446113225d11*b))))))))))))))))))))))))))
   case (27)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.1309733552923256d2+b*(2.443103908609912d3+  &
b*(-3.675663404700058d4+b*(4.275431840067346d5+b*(-4.057932398788864d6+  &
b*(3.2450527666111157d7+b*(-2.231862819333875d8+b*(1.338584266227085d9+  &
b*(-7.067565594557271d9+b*(3.3063031375809d10+b*(-1.3762477730492435d11+  &
b*(5.1097412126250275d11+b*(-1.6938434129836841d12+b*(5.011476958536085d12+  &
b*(-1.3211978621272735d13+b*(3.0945870218960323d13+b*(-6.4116955073823565d13+  &
b*(1.1680549053235025d14+b*(-1.8560171205580507d14+b*(2.5451301621259477d14+  &
b*(-2.9696273928928814d14+b*(2.892193366052213d14+b*(-2.288843247486608d14+  &
b*(1.4146597986124954d14+b*(-6.408761092962503d13+(1.8934677304505707d13-  &
2.7387867918643174d12*b)*b))))))))))))))))))))))))))
   case (28)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.1728612573401895d2+b*(2.631034978502982d3+  &
b*(-4.116743013264065d4+b*(4.988003813411903d5+b*(-4.940091615916878d6+  &
b*(4.130067157505056d7+b*(-2.9758170924451672d8+b*(1.874017972717919d9+  &
b*(-1.0415359823558084d10+b*(5.143138214014734d10+b*(-2.2667610379634597d11+  &
b*(8.942047122093799d11+b*(-3.1618410375695434d12+b*(1.002295391707217d13+  &
b*(-2.845656933812589d13+b*(7.220703051090743d13+b*(-1.6320679473336908d14+  &
b*(3.270553734905807d14+b*(-5.774275486180603d14+b*(8.907955567440816d14+  &
b*(-1.1878509571571525d15+b*(1.3496902374910327d15+b*(-1.2817522185925003d15+  &
b*(9.902618590287467d14+b*(-5.981510353431668d14+b*(2.6508548226307993d14+  &
b*(-7.668603017220088d13+1.08699073414755d13*b))))))))))))))))))))))))))))
   case (29)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.2147491593880535d2+b*(2.8259264583920918d3+  &
b*(-4.591751822486842d4+b*(5.786084423557808d5+b*(-5.969277369232895d6+  &
b*(5.207475981202027d7+b*(-3.922667985495902d8+b*(2.5879295813723644d9+  &
b*(-1.5102271744159224d10+b*(7.850053063496172d10+b*(-3.652003894496685d11+  &
b*(1.5254080384748245d12+b*(-5.730836880594797d12+b*(1.9377710906339527d13+  &
b*(-5.89457507718322d13+b*(1.6107722190894733d14+b*(-3.9441642060564197d14+  &
b*(8.622368937478944d14+b*(-1.6745398909923745d15+b*(2.870341238397596d15+  &
b*(-4.305959719694678d15+b*(5.591573841034279d15+b*(-6.195135723197086d15+  &
b*(5.743518782366731d15+b*(-4.33659500623796d15+b*(2.5624929952097726d15+  &
b*(-1.1119474374969127d15+(3.152273129027895d14-  &
4.3820588669829075d13*b)*b))))))))))))))))))))))))))))
   case (30)  
       dipoleSp(i) = b**5*(2.784163998415854d0+b*(-1.2566370614359172d2+b*(3.0277783482772413d3+  &
b*(-5.101946469429825d4+b*(6.676251257951317d5+b*(-7.163132843079473d6+  &
b*(6.5093449765025335d7+b*(-5.116523459342481d8+b*(3.528994883689588d9+  &
b*(-2.1574673920227463d10+b*(1.177507959524426d11+b*(-5.766321938678976d11+  &
b*(2.542346730791374d12+b*(-1.011324155399082d13+b*(3.6333207949386614d13+  &
b*(-1.1789150154366441d14+b*(3.4516547551917283d14+b*(-9.101917398591738d14+  &
b*(2.155592234369736d15+b*(-4.5669269754337485d15+b*(8.611023715192788d15+  &
b*(-1.4353199065648925d16+b*(2.096840190387854d16+b*(-2.6550581670844653d16+  &
b*(2.8717593911833657d16+b*(-2.601957003742776d16+b*(1.9218697464073293d16+  &
b*(-1.1119474374969127d16+b*(4.728409693541842d15+b*(-1.3146176600948722d15+  &
1.7935347113434574d14*b))))))))))))))))))))))))))))))

 case default
   print *, "Your value of n is not programmed"
 end select
end do
return
end subroutine asyRgSpdip

