
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file contains three main subroutines and three minor subroutines
! - get molecule
! - get normalisation constants
! - calculate interatomic distances
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getmolecule(nobasisfun, noatoms,atomsr,basis,moleculeprint,moleculename,kmax, betachangealg)
  implicit none
  integer, intent(in)               :: nobasisfun, noatoms, kmax, moleculeprint
  real*8, intent(in)                :: betachangealg
  real*8, intent(out)               :: atomsr(1:5,1:noatoms), basis(1:nobasisfun,1:30)
  real*8                            :: basisin(1:nobasisfun,1:30)
  integer                           :: i, st, atid, j, toadd
  character(len=3)                  :: atomtype
  character(len=100)                :: filename
  character(len=40)                 :: basisset
  character(len=30), intent(in)     :: moleculename
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !This program produces data for individual molecules 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! atoms has (atomid, Z, Ax, Ay, Az)  
  ! basis has elements (atom centre, no primitives, exp1, exp2, exp3, coef1, coef2, coef3, type of basis function, blank) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
print *, "get molecule start"  
  filename = 'Data/geom.csv'
  open(unit=18, file = filename)
  read(18,*)   ! doesn't take into account line about noatoms
  read(18,*)   ! doesn't take into account line about noelec
  read(18,*)   ! doesn't take into account line about multiplicity
  do i=1,noatoms
     read(18,*) atomsr(1,i)
     read(18,*) atomsr(2,i)
     read(18,*) atomsr(3,i)
     read(18,*) atomsr(4,i)
     read(18,*) atomsr(5,i)
  end do
  close(18)
  
  st = 0;
  toadd = 0;
  
print *, "here1"
  open(unit = 16, file = "Data/basisset.csv")
  read(16,*) basisset
  close(16)
  
print *, "here2"
  filename = 'Data/geom.csv' 
  open(unit=18, file = filename)
  read(18,*) 
  read(18,*) 
  read(18,*) 
  do atid=1,noatoms
     read(18,*)                ! don't need this line here
     read(18,*) atomtype
     read(18,*)                ! don't need this line here
     read(18,*)                ! don't need this line here
     read(18,*)                ! don't need this line here
     filename = 'Data/'//trim(basisset)
     filename =trim(filename)//'/Basis_fortran/'
     filename = trim(filename)//trim(atomtype)
     filename = trim(filename)//'.bas'
     open(unit = 25, file = filename)
     read(25,*) toadd
     read(25,*) 
     do j=1,toadd
        basis(st+j,1) = atid
        read(25,*) basis(st+j,2) 
        read(25,*) basis(st+j,3) 
        basis(st+j,4) = atomsr(3,atid)  ! Ax
        basis(st+j,5) = atomsr(4,atid)  ! Ay
        basis(st+j,6) = atomsr(5,atid)  ! Az
        read(25,*) 
        do i = 1,int(basis(st+j,2))   
           read(25,*) basis(st+j,10+i)  
           read(25,*) basis(st+j,20+i)
        end do
     end do
     st = st + toadd
     close(25)  
  end do
  
print *, "here3"
  basisin = basis;
  call getnormcons(basis,basisin,nobasisfun,kmax, betachangealg)
 
open(unit=25,file='Integrals/rampids.txt')
do i=1,nobasisfun
  if (int(basis(i,3)) == 0 .or. int(basis(i,3)) == 5) then
    print *, "rampids = ", i
    write(25,"(i5)") i
  end if
end do
close(25)
 
print *, "here4"
  if (moleculeprint ==1) then
     do i=1,nobasisfun
        print "(f5.0,f5.0,f5.0,f10.5,f10.5,f10.5,f10.5,f10.5,f10.5)", basis(i,1:3), basis(i,11:13), basis(i,21:23)
     end do
  end if
  
print *, "here5"
  return
end subroutine getmolecule
!*************************************************************************************



subroutine getnormcons(basis,basisin, nobasisfun,kmax, betachangealg)
  implicit none
  real*8, intent(in)          :: basisin(1:nobasisfun,1:30), betachangealg
  real*8, intent(out)         :: basis(1:nobasisfun,1:30)
  integer, intent(in)         :: nobasisfun,kmax
  real*8                      :: selfSP(1:9*nobasisfun,1:30), selfoverlap(1:nobasisfun)
  integer                     :: index, tot, i
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates the normalisation coefficients of each basis function and 
  ! makes sure each basis function is completely normalised.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call getprelnc(basis,basisin,nobasisfun)
  call getselfSP(selfSP,basis,nobasisfun,tot)
  call getselfoverlap(selfoverlap,selfSP, nobasisfun,tot,kmax, betachangealg)
  
  do i=1,nobasisfun
     if (abs(selfoverlap(i)-1) > 0.01) then 
        print *, ":: WARNING:: There is probably a problem with the basis set somewhere &
             because it wasn't initially normalised ", i, int(basis(i,3)), selfoverlap(i)
        ! call exit
     end if
  end do

  do index=21,30
     basis(1:nobasisfun,index) =basis(1:nobasisfun,index)/sqrt(selfoverlap(1:nobasisfun)) 
  end do
  
  call getselfSP(selfSP,basis,nobasisfun,tot)
  call getselfoverlap(selfoverlap,selfSP, nobasisfun,tot,kmax, betachangealg)
  
  return
end subroutine getnormcons

!*************************************************************
subroutine getprelnc(basis,basisin, nobasisfun) 
  real*8, intent(in)        :: basisin(1:nobasisfun,1:30)
  real*8, intent(out)       :: basis(1:nobasisfun,1:30)
  integer, intent(in)       :: nobasisfun
  real*8                    :: pi, beta, pref1, pref2
  integer                   :: n, i, v
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Multiplies the normalisation constant for each basis function 
  ! to the contraction coefficient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  basis = basisin;
  
  pi = dacos(-1d0);
  pref1 = (2d0/pi)**(3d0/4d0);
  do i=1,nobasisfun
     if (int(basis(i,3)) == 0) then
        n = int(basis(i,11))
        basis(i,21) = sqrt((1d0+n)*(1d0+2d0*n)*(3d0+2d0*n))*basis(i,21)
     elseif (int(basis(i,3)) == 5) then
        n = int(basis(i,11))
        basis(i,21) = sqrt((1d0+n)*(1d0+2d0*n)*(3d0+2d0*n))*basis(i,21)
        beta = basis(i,12)
        pref2 = (beta*2d0/pi)**(3d0/4d0)
        basis(i,22) = basis(i,22)*pref2
     elseif (int(basis(i,3)) == 1) then
      do v =1,10
        beta = basis(i,10+v)
        pref2 = (beta*2d0/pi)**(3d0/4d0)
        basis(i,20+v) = basis(i,20+v)*pref2
      end do
     elseif (int(basis(i,3)) == 2 .or. int(basis(i,3)) == 3 .or. int(basis(i,3)) == 4) then
      do v=1,10
        beta = basis(i,10+v)
        pref2 =  pref1*2d0*beta**(5d0/4d0)  
        basis(i,20+v) = basis(i,20+v)*pref2 
      end do 
     end if
  end do
  
  
  return
end subroutine getprelnc

!*******************************************************************

subroutine getselfsp(selfSP,basis,nobasisfun,tot)
  implicit none
  integer, intent(in)    :: nobasisfun
  integer, intent(out)   :: tot
  real*8, intent(in)     :: basis(1:nobasisfun,1:30)
  real*8, intent(out)    :: selfSP(1:9*nobasisfun,1:30)
  integer                :: id1, k, k2, i
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Constructs basis function pair components from each basis function * itself
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  tot = 0;
  selfSP = 0d0;
  
  do id1=1,nobasisfun
     do k=1,int(basis(id1,2))
        do k2=1,int(basis(id1,2))
           tot = tot+1;
           selfSP(tot,1) = id1
           selfSP(tot,2) = basis(id1,10+k)
           selfSP(tot,3) = basis(id1,10+k2)
           selfSP(tot,4) = basis(id1,20+k)
           selfSP(tot,5) = basis(id1,20+k2)
           selfSP(tot,6) = basis(id1,3)
           selfSP(tot,7) = basis(id1,3)
           if (basis(id1,3) == 5) then
              selfSP(tot,6) = k-1
              selfSP(tot,7) = k2-1
           elseif (basis(id1,3) >=9 .and. basis(id1,3) <=11) then;
              if (k == 1 .and. k2 == 1) then; 
                 selfSP(tot,6) = basis(id1,3) -3
                 selfSP(tot,7) = basis(id1,3) -3
              elseif (k == 2 .and. k2 == 2) then; 
                 selfSP(tot,6) = basis(id1,3) -7
                 selfSP(tot,7) = basis(id1,3) -7
              elseif (k == 1 .and. k2 == 2) then; 
                 selfSP(tot,6) = basis(id1,3) -3
                 selfSP(tot,7) = basis(id1,3) -7
              elseif (k == 2 .and. k2 == 1) then; 
                 selfSP(tot,6) = basis(id1,3) -7
                 selfSP(tot,7) = basis(id1,3) -3
              end if
           end if
        end do
     end do
  end do

  return
end subroutine getselfsp

!*****************************************************************
subroutine getselfoverlap(selfoverlap,selfSP,nobasisfun,tot,kmax, betachangealg)
  implicit none
  integer, intent(in)         :: nobasisfun, tot, kmax
  real*8, intent(in)          :: selfSP(1:9*nobasisfun,1:30), betachangealg
  real*8, intent(out)         :: selfoverlap(1:nobasisfun)
  integer                     :: n, id1, i, k
  real*8                      :: coeftable(0:30), beta, invsqrtbeta, prefactor, h(0:20), coeff, RGSs(1:30,1:1), OL(1:1,1:4), pi, pi2
  real*8                      :: RGPp(1:30,1:1), OL2(1:9,1:4)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Calculates the overlap of a basis function with itself 
  !!  for the purpose of normalising the basis function correctly.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  selfoverlap = 0d0;
  pi = dacos(-1d0);
  pi2 = pi**(3d0/2d0)/2d0;
  
  do i=1,tot
     id1 = int(selfSP(i,1))
     if (int(selfSP(i,6)) == 0 .and. int(selfSP(i,7)) == 0) then
        n = int(selfSP(i,2)+selfSP(i,3))
        selfoverlap(id1) = selfoverlap(id1) +  selfSP(i,4)*selfSP(i,5)*2d0/((3d0+n)*(2d0+n)*(1d0+n))
     elseif (int(selfSP(i,6)) == 1 .and. int(selfSP(i,7)) == 1) then
        beta =(selfSP(i,2)+selfSP(i,3))
        selfoverlap(id1) = selfoverlap(id1) +  selfSP(i,4)*selfSP(i,5)*(pi/beta)**(3d0/2d0)
     elseif (int(selfSP(i,6)) == 1 .and. int(selfSP(i,7)) == 0) then
        beta =selfSP(i,2);
        n = int(selfSP(i,3)); 
        RGSs = 1; 
        RGSs(7,1) = n;
        RGSs(8,1) = beta;
        RGSs(4:6,1) = 0;
        call  newgenrgSsol(OL,1,RGSs,40)
        selfoverlap(id1) = selfoverlap(id1) + OL(1,1)*selfSP(i,4)*selfSP(i,5)/(2d0*sqrt(pi))
     elseif (int(selfSP(i,6)) == 0 .and. int(selfSP(i,7)) == 1) then
        beta =selfSP(i,3)
        n = int(selfSP(i,2)) 
        RGSs = 1
        RGSs(7,1) = n;
        RGSs(8,1) = beta;
        RGSs(4:6,1) = 0;
        call  newgenrgSsol(OL,1,RGSs,40)
        selfoverlap(id1) = selfoverlap(id1) + OL(1,1)*selfSP(i,4)*selfSP(i,5)/(2d0*sqrt(pi))
     elseif (int(selfSP(i,6)) >= 2 .and. int(selfSP(i,6)) <= 4 .and. int(selfSP(i,7)) >= 2 .and. int(selfSP(i,7)) <= 4 ) then
        beta =(selfSP(i,2)+selfSP(i,3))
        selfoverlap(id1) = selfoverlap(id1) +  selfSP(i,4)*selfSP(i,5)*pi2/(beta**(5d0/2d0))
     else 
        print *, "The self overlap of basis function with itself is not able to be calculated"; call exit
     end if
  end do
  return
end subroutine getselfoverlap





!*****************************************************

subroutine getinvinteratomdist(invintatomdist,noatoms,atoms,distprint)
  implicit none
  integer, intent(in)  :: noatoms,distprint
  real*8, intent(in)   :: atoms(1:5,1:noatoms)
  real*8, intent(out)  :: invintatomdist(1:noatoms,1:noatoms)
  real*8               :: R, Ax, Ay, Az, Bx, By, Bz
  integer              :: i, j
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program calculates all interatomic distances
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  invintatomdist = 0d0;
  do i=1,noatoms
     Ax = atoms(3,i)
     Ay = atoms(4,i)
     Az = atoms(5,i)
     do j=i+1,noatoms
        Bx = atoms(3,j)
        By = atoms(4,j)
        Bz = atoms(5,j)
        
        R = sqrt((Bx-Ax)**2 + (By-Ay)**2+(Bz-Az)**2)
        if (R > 0 .and. R < 1) then
           print *, "The atoms are too close - check coordinate inputs are in Bohr";
           print *, R, i, j, Ax, Ay, Az, Bx, By, Bz;
           call exit;
        else if (R>0 .and. atoms(2,i) > 1 .and. atoms(2,j) > 1 .and. R<2) then
           print *, "Two atoms with ramps are closer than 2 Bohr - check coordinate inputs are in Bohr"
           print *, R, i, j, Ax, Ay, Az, Bx, By, Bz;
           call exit;
        end if
        invintatomdist(i,j) = 1/R;
        invintatomdist(j,i) = invintatomdist(i,j)
     end do
  end do
  
  if (distprint ==1) then
     print *, "inverse interatomic distance"
     do i=1,noatoms
        do j=i+1,noatoms
           print "(i3,i3,f20.15)", i, j, invintatomdist(i,j)
        end do
     end do
  end if
  return
end subroutine getinvinteratomdist












