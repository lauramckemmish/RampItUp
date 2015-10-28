subroutine overlap(normolmat,algchange,noRR, noRgSs, noRgSp, &
     RRr, RgSs, RgSpr,  nobasisfun,basis,olprint,moleculename, kmax, pts, &
      ggss, ggsp, ggpp, noggss, noggsp, noggpp)
  implicit none
  integer, intent(in)     :: nobasisfun, kmax, algchange(1:2), olprint, pts
  integer, intent(in)     :: noRR, noRgSs, noRgSp, noggss, noggsp, noggpp
  real*8, intent(in)      :: basis(1:nobasisfun,1:30)
  real*8, intent(in)      :: RRr(1:30,1:noRR),RgSs(1:30,1:noRgSs),RgSpr(1:30,1:noRgSp)
  real*8, intent(in)      :: ggss(1:30,1:noggss), ggsp(1:30,1:noggsp), ggpp(1:30,1:noggpp)
  real*8, intent(out)     :: normolmat(0:nobasisfun,0:nobasisfun,1:2)
  real*8                  :: olmat(1:8*nobasisfun+nobasisfun**2+9*noggpp+3*noggsp+noggss,1:4)
  real*8                  :: OLr2SS(1:noRR,1:4), OLRgSsg(1:noRgSs,1:4), OLRgSpg(1:3*(noRgSp-algchange(2)),1:4)
  real*8                  :: genSpr(1:30,1:(noRgSp-algchange(2)))
  integer                 :: totolread, i, j, noSP, st, fin
  real*8, allocatable, dimension(:,:) :: olggss, olggsp, olggpp
  character(len=30),intent(in)           :: moleculename

#include "DEBUG.h"
  olmat = 0d0; 
 normolmat = 0d0;
  !*******************
  !*********(S|S) !*********
  call r2calcSS(OLr2SS,RRr, noRR);
  OLr2SS(1:noRR,1) = OLr2SS(1:noRR,1)/(4d0*dacos(-1d0));  olmat(1:noRR,1:4) = OLr2SS;
  
  !*********************
  !********* (S|s) via quadrature!*********
  call newgenRgSsol(OLRgSsg, noRgSs, RgSs,pts)
  st = noRR+1;  fin = noRR+noRgSs;
  OLRgSsg(1:noRgSs,1) = OLRgSsg(1:noRgSs,1)/(2d0*sqrt(dacos(-1d0)));  olmat(st:fin,1:4) = OLRgSsg(1:noRgSs,1:4)
  
  !****************************
  
  !********* (S|p) via general algorithm (non-concentric only) !*********
  noSP = noRgSp-algchange(2);  genSpr = RgSpr(1:30,algchange(2)+1:noRgSp)
  call genRgSpol(OLRgSpg,noSP,genSpr,kmax)
  st = noRR+noRgSs+1;  fin = noRR+noRgSs+3*(noRgSp-algchange(2));
  OLRgSpg(1:3*noSP,1) = OLRgSpg(1:3*noSP,1)/(2d0*sqrt(dacos(-1d0)));
  olmat(st:fin,1:4) = OLRgSpg(1:3*noSP,1:4)

  !******* Reads integrals from Qchem
 allocate(olggss(1:noggss,1:4), olggsp(1:3*noggsp,1:4), olggpp(1:9*noggpp,1:4))
  call calcggssol(OLggss,noggss,ggss);  st = fin; fin = st + noggss;  olmat(st+1:fin,1:4) = OLggss;
  call calcggspol(OLggsp,noggsp,ggsp);  st = fin; fin = st + 3*noggsp;  olmat(st+1:fin,1:4) = OLggsp;
  call calcggppol(OLggpp,noggpp,ggpp);  st = fin; fin = st + 9*noggpp;  olmat(st+1:fin,1:4) = OLggpp;
 deallocate(olggss, olggsp,olggpp)

  !************ Puts all integral together with normalisation constants into the overlap matrix****
  call overlapdigestion(normolmat,olmat, nobasisfun,fin,noggss,noggsp,noggpp)

#ifdef WRITEONE
  call writeoverlap(normolmat,nobasisfun,basis,moleculename)
#endif
 
  !******* Printing options*****
#ifdef ONE_ELEC_PRT
    !print *, "These are primitive OL integrals calculated in Fortran"
    ! print "(a1,a3,a1,a3,a3,a20,a20)", "(", "bf1", ",", "bf2", ") = ", "unnorm ol", "norm factor (with contra coeff)"
    ! do i=1,fin
    !    print "(a1,i3,a1,i3,a3,f20.15,f20.15)", "(", int(olmat(i,2)), ",", int(olmat(i,3)), ") = ", olmat(i,1), olmat(i,4)
    ! end do
     
     !print *, "These are primitive OL integrals calculated in Qchem"
     !print "(a1,a3,a1,a3,a3,a20,a20)", "(", "bf1", ",", "bf2", ") = ", "norm ol", "contraction coef"
     !do i=noRR+noRgSs+3*(noRgSp-algchange(2))+1,noRR+noRgSs+3*(noRgSp-algchange(2))+totolread
        !   print "(a1,i3,a1,i3,a3,f20.15,f20.15)", "(", int(olmat(i,2)), ",", int(olmat(i,3)), ") = ", olmat(i,1), olmat(i,4)
     !end do
     
   !  print *, ""
   !  print *, "These are digested OL integrals"
    ! do i=1,nobasisfun
    !    do j=i,nobasisfun
        !   if (normolmat(i,j,2) == 1) then
           !   print "(a1,i3,a1,i3,a3,f20.15)", "(", i, ",", j, ") = ", normolmat(i,j,1)
        !   end if
    !    end do
    ! end do
     
     print *, ""
     
     print *, "This is the overlap matrix" 
     do j=1,nobasisfun,4
        do i=1,nobasisfun
             print "(4f20.15)", normolmat(i,j:min(nobasisfun,j+3),1)
        end do
        print *,""
     end do

#endif
  
  return
end subroutine overlap


!********************************************

subroutine overlapdigestion(normolmat,olmat, nobasisfun,nonzero,noggss,noggsp,noggpp)
  implicit none
  real*8, intent(in)           :: olmat(1:8*nobasisfun+nobasisfun**2+9*noggpp+3*noggsp+noggss,1:4)
  real*8, intent(out)          :: normolmat(0:nobasisfun,0:nobasisfun,1:2)
  integer, intent(in)          :: nobasisfun, nonzero,noggss,noggsp,noggpp
  integer                      :: i, bf1, bf2
  
  normolmat = 0d0;
  
  do i=1,nonzero
     bf1 = int(olmat(i,2));
     bf2 = int(olmat(i,3));
if (bf1 == 1 .and. bf2 == 1) then; print *, olmat(i,1),olmat(i,4); end if
     normolmat(bf1,bf2,1) = normolmat(bf1,bf2,1)+olmat(i,1)*olmat(i,4) 
     normolmat(bf2,bf1,1) = normolmat(bf1,bf2,1)   
  end do
  
  return
end subroutine overlapdigestion
