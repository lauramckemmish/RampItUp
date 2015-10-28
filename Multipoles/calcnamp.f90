subroutine calcnamp(namp,flatdigmp,nobasisfun,noatoms,Tintmat,nLpure,Lcalc)
  implicit none
  integer, intent(in)                     :: nobasisfun,noatoms, nLpure, Lcalc
  real*8, intent(in)                      :: flatdigmp(1:nobasisfun*(nobasisfun+1)/2,1:nLpure+2)
  real*8, intent(in)                      :: Tintmat(1:nLpure,1:nLpure,0:noatoms,0:noatoms)
  real*8, intent(out)                     :: namp(1:nobasisfun,1:nobasisfun,1:noatoms)
  real*8                                  :: nampint(1:nobasisfun*(nobasisfun+1)/2,1:noatoms)
  integer                                 :: multipolecenter, i, j, atid, nuclei, k, bf1, bf2
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! This program calculates nuclear attraction integrals via multipole moment interactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 nampint = 0d0;
  
  do nuclei = 1,noatoms
     j = 0;
     do bf1 = 1,nobasisfun;   do bf2 = bf1,nobasisfun
        j = j + 1;
        multipolecenter = int(flatdigmp(j,nLpure+1))
        do k=1,(Lcalc+1)**2
           nampint(j,nuclei) = nampint(j,nuclei) + Tintmat(1,k, nuclei, multipolecenter)*flatdigmp(j,k)
        end do;   
     end do;  end do;  
  end do
  
  j=0;
  do bf1 = 1,nobasisfun; do bf2 = bf1,nobasisfun
     j=j+1;
     namp(bf1,bf2,1:noatoms) = nampint(j,1:noatoms)
     namp(bf2,bf1,1:noatoms) = nampint(j,1:noatoms)
  end do; end do
  
  return
end subroutine calcnamp
