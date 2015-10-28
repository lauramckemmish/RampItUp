program main
   implicit none
  integer, parameter :: maxrampdegree=30, maxL=4;
  integer                     :: n1, n2, i,j,l
  real*8                      :: pt1, pt2, pt3, r1, r2,r3,  pi, ntot, f1, f2, pf
  real*8       :: rrrrpref(0:maxrampdegree,0:maxrampdegree,0:maxL)

call calcrrrrpref(rrrrpref,maxrampdegree,maxL)

do i=0,maxrampdegree
 do j=0,maxrampdegree
  do l=0,maxL
   if ((rrrrpref(i,j,l)) .gt. 1d-10) then
   print *, i,j,l,rrrrpref(i,j,l)
   end if
  end do
 end do
end do

end
