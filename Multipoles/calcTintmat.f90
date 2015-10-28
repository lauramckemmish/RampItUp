subroutine calcTintmat(Tintmat,L,nLpure,noatoms,atoms,invintatomdist,Tintprint)
  implicit none
  integer, intent(in)              :: noatoms, L, nLpure, Tintprint
  real*8, intent(in)               :: atoms(1:5,1:noatoms), invintatomdist(1:noatoms,1:noatoms)
  real*8, intent(out)              :: Tintmat(1:nLpure,1:nLpure,0:noatoms,0:noatoms)
  real*8                           :: J(1:noatoms,1:noatoms,0:10,-10:10)
  integer                          :: m, k
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates the unit multipole-unit multiple interaction (real pure multipoles) 
  ! for each pair of nuclei. 
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  Tintmat = 0d0;
  
  ! Note that the indexing for J is very different from the indexing of Tintmat
print *, "HH"
  call calcJlm(J,10,atoms,noatoms, invintatomdist)                                             !  Compute Jlm array
print *, "AAA"
  call calcImat(Tintmat,J,L,nLpure,noatoms)                                          !  Scatter J into I
print *, "DDD"  
  if (Tintprint ==1) then
     !********** Print out Tintmat ******************
     print *, "Tintmat between atoms 1 & 2"
     do m=1,9; do k=1,9
        print "(a1,i3,a1,i3,a3,f20.15,f20.15)", "(", m, ",", k, ") = ", Tintmat(m,k,1,2), Tintmat(m,k,2,1)
     end do; end do;
  end if
  
  return
end subroutine calcTintmat


!***************************************************************************************************

subroutine calcJlm(J,L,atoms,noatoms, invintatomdist)
  implicit none
  integer               :: i, m, k
  integer , intent(in)  :: L, noatoms
  real*8  , intent(in)  :: atoms(1:5,1:noatoms), invintatomdist(1:noatoms,1:noatoms)
  real*8  , intent(out) :: J(1:noatoms,1:noatoms,0:L,-L:L)
  real*8                :: Rx(1:noatoms,1:noatoms), Ry(1:noatoms,1:noatoms)
  real*8                :: Rz(1:noatoms,1:noatoms), rR2(1:noatoms,1:noatoms), RxrR2(1:noatoms,1:noatoms)
  real*8                :: RyrR2(1:noatoms,1:noatoms), RzrR2(1:noatoms,1:noatoms)
  real*8                :: xtemp(1:noatoms,1:noatoms), ytemp(1:noatoms,1:noatoms), ztemp(1:noatoms,1:noatoms)
  real*8                :: rtemp(1:noatoms,1:noatoms)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculates RY_{l,m} with Peter's definition of the multipole moments
!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  J = 0d0;
  Rx = 0d0;
  Ry = 0d0;
  Rz = 0d0;
  do i=1,noatoms
     do k=i+1,noatoms
        Rx(i,k) = atoms(3,k) - atoms(3,i)
        Rx(k,i) = -Rx(i,k)
        Ry(i,k) = atoms(4,k) - atoms(4,i)
        Ry(k,i) = -Ry(i,k)
        Rz(i,k) = atoms(5,k) - atoms(5,i)
        Rz(k,i) = -Rz(i,k)
     end do
  end do
  
  rR2  = invintatomdist**2;
  J(1:noatoms,1:noatoms,0,0) = invintatomdist(1:noatoms,1:noatoms)                           !  Compute 1 / R        8+D+S
  RxrR2 = Rx * rR2
  RyrR2 = Ry * rR2
  RzrR2 = Rz * rR2                                                  !  Compute Ri / R*R       3
  J(1:noatoms,1:noatoms,1,-1) =  RyrR2 * J(1:noatoms,1:noatoms,0,0)
  J(1:noatoms,1:noatoms,1, 0) = -RzrR2 * J(1:noatoms,1:noatoms,0,0)
  J(1:noatoms,1:noatoms,1, 1) =  RxrR2 * J(1:noatoms,1:noatoms,0,0)                                         !  Compute J(1,m)         3
  do m = 2,L
     xtemp =   (2*m-1)   * RxrR2                                    !                        L-1
     ytemp =   (2*m-1)   * RyrR2                                    !                        L-1
     ztemp =   (2*m-1)   * RzrR2                                    !                        L-1
     rtemp = (m-1)*(m-1) *   rR2                                    !                        L-1
     J(1:noatoms,1:noatoms,m, 0  ) = - ztemp * J(1:noatoms,1:noatoms,m-1, 0  ) - rtemp * J(1:noatoms,1:noatoms,m-2, 0 )         !  Compute J(m, 0  )   3(L-1)
     J(1:noatoms,1:noatoms,m, m  ) =   xtemp * J(1:noatoms,1:noatoms,m-1,+m-1) - ytemp * J(1:noatoms,1:noatoms,m-1,-m+1)        !  Compute J(m,+m  )   3(L-1)
     J(1:noatoms,1:noatoms,m,-m  ) =   ytemp * J(1:noatoms,1:noatoms,m-1,+m-1) + xtemp * J(1:noatoms,1:noatoms,m-1,-m+1)        !  Compute J(m,-m  )   3(L-1)
     J(1:noatoms,1:noatoms,m,+m-1) = - ztemp * J(1:noatoms,1:noatoms,m-1,+m-1)                                                  !  Compute J(m,+m-1)     L-1
     J(1:noatoms,1:noatoms,m,-m+1) = - ztemp * J(1:noatoms,1:noatoms,m-1,-m+1)                                                  !  Compute J(m,-m+1)     L-1
  end do
  do m = 1,L-2
     do i = m+2,L
        ztemp =      (2*i-1)      * RzrR2    
        rtemp = ((i-1)*(i-1)-m*m) *   rR2                                                                                       !                    (L-1)(L-2)
        J(1:noatoms,1:noatoms,i,+m) = - ztemp * J(1:noatoms,1:noatoms,i-1,+m) - rtemp * J(1:noatoms,1:noatoms,i-2,+m)           !  Compute J(i,+m)
        J(1:noatoms,1:noatoms,i,-m) = - ztemp * J(1:noatoms,1:noatoms,i-1,-m) - rtemp * J(1:noatoms,1:noatoms,i-2,-m)           !  Compute J(i,-m)  3(L-1)(L-2)
     end do
  end do
  
  
end subroutine calcJlm

!***************************************************************************************************

subroutine calcImat(Tintmat,J,L,nLpure,noatoms)
   implicit none
   integer , intent(in)  :: L, nLpure, noatoms
   real*8  , intent(in)  :: J(1:noatoms,1:noatoms,0:10,-10:10)
   real*8  , intent(out) :: Tintmat(1:nLpure,1:nLpure,0:noatoms,0:noatoms)
   real*8                :: I((1+10)**2,(1+10)**2,0:noatoms,0:noatoms)
   integer               :: k,m

print *, "GGGG"
!I = 0d0;

!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!  L = 0  !!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!! 

! Charge-charge interaction (0 with 0)
   I(1,1,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,0, 0);  


print *, "ssss"
!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!  L = 1  !!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!! 

! Charge-dipole interaction (1 with 0)
   I(1,2,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,1,-1);          
   I(1,3,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,1, 0);  
   I(1,4,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,1, 1);

do k=2,4;
   I(k,1,1:noatoms,1:noatoms) = I(1,k,1:noatoms,1:noatoms)
end do

print *, "ssss"

!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!  L = 2  !!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!! 

! Charge-quadrupole interactions (0 with 2) 
   I(1,5,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,2,-2);
   I(1,6,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,2,-1);
   I(1,7,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,2, 0);
   I(1,8,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,2, 1);
   I(1,9,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,2, 2);
  
do k=5,9;
   I(k,1,1:noatoms,1:noatoms) = I(1,k,1:noatoms,1:noatoms)
end do



! Dipole-dipole interaction  ( 1 with 1) 
   I(2,2,1:noatoms,1:noatoms) = -2*J(1:noatoms,1:noatoms,2, 2)-2*J(1:noatoms,1:noatoms,2,0); 
   I(2,3,1:noatoms,1:noatoms) =  2*J(1:noatoms,1:noatoms,2,-1); 
   I(2,4,1:noatoms,1:noatoms) =  2*J(1:noatoms,1:noatoms,2,-2);
   I(3,2,1:noatoms,1:noatoms) =  2*J(1:noatoms,1:noatoms,2,-1);     
   I(3,3,1:noatoms,1:noatoms) =    J(1:noatoms,1:noatoms,2, 0); 
   I(3,4,1:noatoms,1:noatoms) =  2*J(1:noatoms,1:noatoms,2, 1);
   I(4,2,1:noatoms,1:noatoms) =  2*J(1:noatoms,1:noatoms,2,-2);        
   I(4,3,1:noatoms,1:noatoms) =  2*J(1:noatoms,1:noatoms,2, 1); 
   I(4,4,1:noatoms,1:noatoms) =  2*J(1:noatoms,1:noatoms,2, 2)-2*J(1:noatoms,1:noatoms,2,0);


!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!  L = 3  !!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!! 

! Charge-octopole interactions (0 with 3)
   I(1,10,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-3);
   I(1,11,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-2);
   I(1,12,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-1);
   I(1,13,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,3, 0);
   I(1,14,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3, 1);
   I(1,15,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3, 2);
   I(1,16,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3, 3);
  
do k=10,16;
   I(k,1,1:noatoms,1:noatoms) = I(1,k,1:noatoms,1:noatoms)
end do

! Dipole-quadrupole interactions (1 with 2)
   I(2,5,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,3, 1) - 2*J(1:noatoms,1:noatoms,3, 3);
   I(2,6,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,3, 0) - 2*J(1:noatoms,1:noatoms,3, 2);
   I(2,7,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-1);
   I(2,8,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-2);
   I(2,9,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-1) + 2*J(1:noatoms,1:noatoms,3,-3);

   I(3,5,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-2);
   I(3,6,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-1);
   I(3,7,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,3, 0);
   I(3,8,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3, 1);
   I(3,9,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3, 2);

   I(4,5,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,3,-1) + 2*J(1:noatoms,1:noatoms,3,-3);
   I(4,6,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3,-2);
   I(4,7,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,3, 1);
   I(4,8,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,3, 0) + 2*J(1:noatoms,1:noatoms,3, 2);
   I(4,9,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,3, 1) + 2*J(1:noatoms,1:noatoms,3, 3);

do k=2,4; do m=5,9;
   I(m,k,1:noatoms,1:noatoms) = I(k,m,1:noatoms,1:noatoms)
end do; end do



!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!  L = 4   !!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!
!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!

! Charge-4-pole interactions (0 with 4)
   I(1,17,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-4);
   I(1,18,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-3);
   I(1,19,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(1,20,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1);
   I(1,21,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,4, 0);
   I(1,22,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 1);
   I(1,23,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 2);
   I(1,24,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 3);
   I(1,25,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 4);
  
do k=17,25;
   I(k,1,1:noatoms,1:noatoms) = I(1,k,1:noatoms,1:noatoms)
end do



! 1-pole with 3-pole

   I(2,10,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 4) - 2*J(1:noatoms,1:noatoms,4, 2);
   I(2,11,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 3) - 2*J(1:noatoms,1:noatoms,4, 1);
   I(2,12,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 2) - 2*J(1:noatoms,1:noatoms,4, 0);
   I(2,13,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1);
   I(2,14,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(2,15,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1) + 2*J(1:noatoms,1:noatoms,4,-3);
   I(2,16,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2) + 2*J(1:noatoms,1:noatoms,4,-4);
   
   I(3,10,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-3);
   I(3,11,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(3,12,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1);
   I(3,13,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,4, 0);
   I(3,14,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 1);
   I(3,15,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 2);
   I(3,16,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 3);

   I(4,10,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-4) - 2*J(1:noatoms,1:noatoms,4,-2);
   I(4,11,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-3) - 2*J(1:noatoms,1:noatoms,4,-1);
   I(4,12,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(4,13,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 1);
   I(4,14,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 0) + 2*J(1:noatoms,1:noatoms,4, 2);
   I(4,15,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 1) + 2*J(1:noatoms,1:noatoms,4, 3);
   I(4,16,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 2) + 2*J(1:noatoms,1:noatoms,4, 4);

do k=2,4; do m=10,16;
   I(m,k,1:noatoms,1:noatoms) = I(k,m,1:noatoms,1:noatoms)
end do; end do



! Quadrupole-quadrupole interactions (2 with 2)
   I(5,5,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 0) - 2*J(1:noatoms,1:noatoms,4, 4);
   I(5,6,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 1) - 2*J(1:noatoms,1:noatoms,4, 3);
   I(5,7,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(5,8,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4,-1) + 2*J(1:noatoms,1:noatoms,4,-3);
   I(5,9,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-4);

   I(6,5,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 1) - 2*J(1:noatoms,1:noatoms,4, 3);
   I(6,6,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 0) - 2*J(1:noatoms,1:noatoms,4, 2);
   I(6,7,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1);
   I(6,8,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(6,9,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1) + 2*J(1:noatoms,1:noatoms,4,-3);

   I(7,5,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(7,6,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1);
   I(7,7,1:noatoms,1:noatoms) =   J(1:noatoms,1:noatoms,4, 0);
   I(7,8,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 1);
   I(7,9,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 2);

   I(8,5,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4,-1) + 2*J(1:noatoms,1:noatoms,4,-3);
   I(8,6,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-2);
   I(8,7,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 1);
   I(8,8,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 0) + 2*J(1:noatoms,1:noatoms,4, 2);
   I(8,9,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 1) + 2*J(1:noatoms,1:noatoms,4, 3);

   I(9,5,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-4);
   I(9,6,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4,-1) + 2*J(1:noatoms,1:noatoms,4,-3);
   I(9,7,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 2);
   I(9,8,1:noatoms,1:noatoms) =-2*J(1:noatoms,1:noatoms,4, 1) + 2*J(1:noatoms,1:noatoms,4, 3);
   I(9,9,1:noatoms,1:noatoms) = 2*J(1:noatoms,1:noatoms,4, 0) + 2*J(1:noatoms,1:noatoms,4, 4);


!**************!!!!!!!!!!!!!!!!!!!!!!!!!!************************!!!!!!!!!!!!!!!!!!!!!



I( 2: 4,1:nLpure,1:noatoms,1:noatoms) = -I(2:4,1:nLpure,1:noatoms,1:noatoms)
I(10:16,1:nLpure,1:noatoms,1:noatoms) = -I(10:16,1:nLpure,1:noatoms,1:noatoms)

!print *, "END OF CALCIMAT";
!print *, "noatoms = ", noatoms, "nLpure = ", nLpure
Tintmat = I(1:nLpure,1:nLpure,0:noatoms,0:noatoms);
!print *, "AAAA"
   return
end subroutine calcImat
















