subroutine calcconstants(maxL, maxrampdegree, RRRRpreftable)
  implicit none
  integer, intent(in)     :: maxL, maxrampdegree
  real*8, intent(out)     :: rrrrpreftable(0:maxrampdegree,0:maxrampdegree,0:maxL)

  call calcrrrrpref(rrrrpreftable,maxrampdegree,maxL)
  
  return
end subroutine calcconstants


!***************************************************************

subroutine calcrrrrpref(rrrrpref,maxrampdegree,maxL)
  implicit none
  integer, intent(in)         :: maxrampdegree, maxL
  real*8, intent(out)         :: rrrrpref(0:maxrampdegree,0:maxrampdegree,0:maxL)
  integer                     :: n1, n2, L
  real*8                      :: pt1, pt2, pt3, r1, r2,r3,  pi, ntot, f1, f2, pf
  
  pi = dacos(-1d0)
  rrrrpref=100000d0;
  
  do n1=0,maxrampdegree
     do n2=0,maxrampdegree
        ntot = dble(n1+n2);
        do L=0,maxL
           pf = (4d0*pi/(2*l+1));
           f1 = (16d0+n1**2+2d0*l*(2d0+n1)*(2d0+n2)+3d0*n1*(3d0+n2)+n2*(9d0+n2))/((1d0+n1)*(2d0+n1)*(1d0+n2)*(2d0+n2));
           f2 = exp(dlgama(3d0+ntot)+ dlgama(3d0+2d0*L) - dlgama(6d0+2d0*L+ntot))
           rrrrpref(n1,n2,L) = pf*f1*f2
        end do
     end do
  end do
  
  return
end subroutine calcrrrrpref
