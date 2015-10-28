subroutine calcusefulSPinter(algchange,noRR, noRGSs, noRGSp, &
 noggss, noggsp, noggpp,RR, RGSs, RGSp,ggss,ggsp,ggpp)
 implicit none
 real*8, intent(inout)            :: RR(1:30,1:noRR), RGSs(1:30,1:noRGSs), RGSp(1:30,1:noRGSp)
 real*8, intent(inout)            :: ggss(1:30,1:noggss),  ggsp(1:30,1:noggsp), ggpp(1:30,1:noggpp)
 integer, intent(in)              :: noRR, noRGSs, noRGSp, noggss, noggsp,noggpp, algchange(1:30)
 integer                          :: i, ntot
 real*8                           :: pi, pf, zeta, cons1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!         This subroutine is used to calculate useful shell pair intermediate quantites
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pi  = dacos(-1d0);

 do i=1,noRR
   ntot = int(RR(4,i)+RR(5,i))
   RR(11,i) = 1d0/((1d0+ntot)*(2d0+ntot)*(3d0+ntot))
 end do
 
cons1 = 4d0/sqrt(dacos(-1d0));
 do i=1,noggss
  zeta = ggss(8,i)
  ggss(11,i) = 4d0*pi**(5d0/2d0)/(zeta**(3d0/2d0))

! Used in genSSss.f90
  ggss(12,i) = 1d0/sqrt(ggss(8,i))
  ggss(13,i) =  cons1*exp(-zeta)*sqrt(zeta);

  ggss(14,i) = 1d0/(zeta**(3d0/2d0))
 end do



return
end subroutine