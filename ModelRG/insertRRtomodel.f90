subroutine insertRRtomodel(shellmodelRGmat,lenshellmodel,shellmodelids,noRR,RR)
  implicit none
  integer, intent(in)            :: noRR
  real*8, intent(out)            :: shellmodelRGmat(1:noRR,1:3,1:1), shellmodelids(1:noRR,1:4)
  integer, intent(out)           :: lenshellmodel(1:noRR)
  real*8, intent(in)             :: RR(1:30,1:noRR)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This subroutine inserts the (SS| into the array that has all the different kinds of ramp models 
  !! so it can be dealt with in the same way
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  shellmodelids(1:noRR,1) = RR(1,1:noRR)
  shellmodelids(1:noRR,2) = RR(2,1:noRR)
  shellmodelids(1:noRR,3) = RR(3,1:noRR)
  shellmodelids(1:noRR,4) = RR(6,1:noRR)                     ! nc of ramp
  
  shellmodelRGmat(1:noRR,1,1) = 1;
  shellmodelRGmat(1:noRR,2,1) = RR(4,1:noRR)+RR(5,1:noRR)     ! ntot of ramp
  shellmodelRGmat(1:noRR,3,1) = 1;
  
  lenshellmodel(1:noRR) = 1;
  
  return
end subroutine insertRRtomodel

!**************************************************
subroutine insertRRSPtomodel(shellmodelRGmat,lenshellmodel,shellmodelids,noRRSP,RRSP)
  implicit none
  integer, intent(in)            :: noRRSP
  real*8, intent(out)            :: shellmodelRGmat(1:3*noRRSP,1:3,1:1), shellmodelids(1:3*noRRSP,1:4)
  integer, intent(out)           :: lenshellmodel(1:3*noRRSP)
  real*8, intent(in)             :: RRSP(1:30,1:noRRSP)
  integer                        :: i
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This subroutine inserts the (SS| into the array that has all the different kinds of ramp models 
  !! so it can be dealt with in the same way
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
do i=1,noRRSP
  shellmodelids(3*i-2:3*i,1) = RRSP(1,i)
  shellmodelids(3*i-2    ,2) = RRSP(2,i)
  shellmodelids(3*i-1    ,2) = RRSP(2,i)+1
  shellmodelids(3*i      ,2) = RRSP(2,i)+2
  shellmodelids(3*i-2:3*i,3) = RRSP(3,i)
  shellmodelids(3*i-2    ,4) = -RRSP(6,i)                       ! nc of ramp
  shellmodelids(3*i-1    ,4) = -RRSP(6,i)                       ! nc of ramp
  shellmodelids(3*i      ,4) =  RRSP(6,i)                       ! nc of ramp
  
  shellmodelRGmat(3*i-2:3*i,1,1) = 1;
  shellmodelRGmat(3*i-2:3*i,2,1) = RRSP(4,i)+RRSP(5,i)     ! ntot of ramp
  shellmodelRGmat(3*i-2    ,3,1) = 4;                                    ! k angular momentum of ramp                   
  shellmodelRGmat(3*i-1    ,3,1) = 2;                                    ! k angular momentum of ramp      
  shellmodelRGmat(3*i      ,3,1) = 3;                                    ! k angular momentum of ramp  
end do

  lenshellmodel(1:3*noRRSP) = 1;
  
  return
end subroutine insertRRSPtomodel



!**************************************************
subroutine insertRRPPtomodel(shellmodelRGmat,lenshellmodel,shellmodelids,noRRPP,RRPP)
  implicit none
  integer, intent(in)            :: noRRPP
  real*8, intent(out)            :: shellmodelRGmat(1:20*noRRPP,1:3,1:1), shellmodelids(1:20*noRRPP,1:4)
  integer, intent(out)           :: lenshellmodel(1:20*noRRPP)
  real*8, intent(in)             :: RRPP(1:30,1:noRRPP)
  integer                        :: i
  real*8                         :: sqrtpi
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This subroutine inserts the (SS| into the array that has all the different kinds of ramp models 
  !! so it can be dealt with in the same way
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sqrtpi = sqrt(dacos(-1d0))

do i=1,noRRPP

! PxPx model
  shellmodelids(20*i-19:20*i-15,1) = RRPP(1,i)
  shellmodelids(20*i-19:20*i-15,2) = RRPP(2,i)

! PxPy model
  shellmodelids(20*i-14    ,1) = RRPP(1,i)
  shellmodelids(20*i-14    ,2) = RRPP(2,i)+1

! PxPz model
  shellmodelids(20*i-13    ,1) = RRPP(1,i)
  shellmodelids(20*i-13    ,2) = RRPP(2,i)+2

! PyPx model
  shellmodelids(20*i-12    ,1) = 0 !RRPP(1,i)+1
  shellmodelids(20*i-12    ,2) = 0 !RRPP(2,i)

! PyPy model
  shellmodelids(20*i-11:20*i-7,1) = RRPP(1,i)+1
  shellmodelids(20*i-11:20*i-7,2) = RRPP(2,i)+1

! PyPz model
  shellmodelids(20*i- 6    ,1) = RRPP(1,i)+1
  shellmodelids(20*i- 6    ,2) = RRPP(2,i)+2

! PzPx model
  shellmodelids(20*i- 5    ,1) = 0!RRPP(1,i)+2
  shellmodelids(20*i- 5    ,2) = 0!RRPP(2,i)

! PzPy model
  shellmodelids(20*i- 4    ,1) = 0!RRPP(1,i)+2
  shellmodelids(20*i- 4    ,2) = 0!RRPP(2,i)+1

! PzPz model
  shellmodelids(20*i-3:20*i,1) = RRPP(1,i)+2
  shellmodelids(20*i-3:20*i,2) = RRPP(2,i)+2


  shellmodelids(20*i-19:20*i,3) = RRPP(3,i)
  shellmodelids(20*i-19:20*i,4) = RRPP(6,i)*2d0*sqrtpi                     ! nc of ramp
  
  shellmodelRGmat(20*i-19:20*i,2,1) = RRPP(4,i)+RRPP(5,i)     ! ntot of ramp
  shellmodelRgmat(20*i-19, 2, 1) = RRPP(4,i)+RRPP(5,i)+2
  shellmodelRgmat(20*i-18, 2, 1) = RRPP(4,i)+RRPP(5,i)+1
  shellmodelRgmat(20*i-11, 2, 1) = RRPP(4,i)+RRPP(5,i)+2
  shellmodelRgmat(20*i-10, 2, 1) = RRPP(4,i)+RRPP(5,i)+1
  shellmodelRgmat(20*i- 3, 2, 1) = RRPP(4,i)+RRPP(5,i)+2
  shellmodelRgmat(20*i- 2, 2, 1) = RRPP(4,i)+RRPP(5,i)+1

  shellmodelRGmat(20*i-16,1,1) = -0.12615662610100800241d0;            ! contraction coefficient
  shellmodelRGmat(20*i-15,1,1) = 0.21850968611841581411d0;             ! contraction coefficient
  shellmodelRGmat(20*i-17,1,1) = 0.28209479177387814347d0;             ! contraction coefficient  
  shellmodelRGmat(20*i-18,1,1) = -2d0*0.28209479177387814347d0;             ! contraction coefficient   
  shellmodelRGmat(20*i-19,1,1) = 0.28209479177387814347d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i-16,3,1) = 7;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i-15,3,1) = 9;                                    ! k angular momentum of ramp   
  shellmodelRGmat(20*i-17,3,1) = 1;                                    ! k angular momentum of ramp 
  shellmodelRGmat(20*i-18,3,1) = 1;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i-19,3,1) = 1;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i-14,1,1) = 0.21850968611841581411d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i-14,3,1) = 5;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i-13,1,1) = -0.21850968611841581411d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i-13,3,1) = 8;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i-12,1,1) = 0.21850968611841581411d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i-12,3,1) = 5;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i- 8,1,1) = -0.12615662610100800241d0;            ! contraction coefficient
  shellmodelRGmat(20*i- 7,1,1) = -0.21850968611841581411d0;            ! contraction coefficient
  shellmodelRGmat(20*i- 9,1,1) = 0.28209479177387814347d0;             ! contraction coefficient
  shellmodelRGmat(20*i-10,1,1) = -2d0*0.28209479177387814347d0;             ! contraction coefficient   
  shellmodelRGmat(20*i-11,1,1) = 0.28209479177387814347d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 8,3,1) = 7;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i- 7,3,1) = 9;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i- 9,3,1) = 1;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i- 10,3,1) = 1;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i- 11,3,1) = 1;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i- 6,1,1) = -0.21850968611841581411d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 6,3,1) = 6;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i- 5,1,1) = 0.21850968611841581411d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 5,3,1) = 8;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i- 4,1,1) = 0.21850968611841581411d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 4,3,1) = 6;                                    ! k angular momentum of ramp  

  shellmodelRGmat(20*i- 0,1,1) = 0.25231325220201600482d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 0,3,1) = 7;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i- 1,1,1) = 0.28209479177387814347d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 1,3,1) = 1;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i- 2,1,1) = -2d0*0.28209479177387814347d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 2,3,1) = 1;                                    ! k angular momentum of ramp  
  shellmodelRGmat(20*i- 3,1,1) = 0.28209479177387814347d0;             ! contraction coefficient                  
  shellmodelRGmat(20*i- 3,3,1) = 1;                                    ! k angular momentum of ramp  
end do

  lenshellmodel(1:20*noRRPP) = 1;
   
  return
end subroutine insertRRPPtomodel
