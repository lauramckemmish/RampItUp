subroutine concmdlctrl(shellmodelRGmat,lenshellmodel,shellmodelids, no,&
     RGin, lendmdlconc, mdlchk)
  implicit none
  real*8, intent(out)            :: shellmodelRGmat(1:no,1:3,1: lendmdlconc), shellmodelids(1:no,1:4)
  integer, intent(out)           :: lenshellmodel(1:no)
  integer, intent(in)            :: no,  lendmdlconc, mdlchk
  real*8, intent(in)             :: RGin(1:30,1:no)
  real*8                         :: RGr(1:30,1:no),  Rs(1:no), n(1:no), beta
  integer                        :: K, i, j, rn, m, index, L, len, actn
  integer                        :: mdlbasis(1:3,1:lendmdlconc), fin, st
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function chooses the ramp basis functions that are to be used in the modelling.
  ! This should probably be modified in the future to suggest the optimal ramp basis functions for a specific
  ! class of Rg shell pairs.
  !
  ! Then, we call the function that finds the coefficients
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  lenshellmodel = 0d0; shellmodelids = 0d0; shellmodelRGmat = 0d0;   
   
   !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
   !     MODEL BASIS FUNCTION SELECTION
   !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
do i=1,no
  beta = Rgin(8,i)
  actn = Rgin(7,i)

  index=0;

select case (actn)

case(4)
! Lithium
 if (beta .lt. 0.01d0) then
  index=3; mdlbasis(1,1:index) = (/4,5,6/)
 elseif (beta .lt. 0.04d0) then
  index=4; mdlbasis(1,1:index) = (/4,5,6,7/)
 elseif (beta .lt. 0.08d0) then
  index=5; mdlbasis(1,1:index) = (/4,5,6,7,8/)
 elseif (beta .lt. 0.65d0 .and. beta .gt. 0.6d0) then
  index=6; mdlbasis(1,1:index) = (/4,5,6,7,8,9/)
 elseif (beta .gt. 2.3d0 .and. beta .lt. 2.4d0) then
  index=7; mdlbasis(1,1:index) = (/4 , 5 , 6 , 7 , 9, 11, 15/)
 elseif (beta .gt. 1.4d0 .and. beta .lt. 1.5d0) then
  index=6; mdlbasis(1,1:index) = (/ 4  ,5 , 6 , 8, 11, 14/);
 end if

case(5)
! Be
 if (beta .lt. 0.03d0) then
   index=4; mdlbasis(1,1:index) = (/5,6,7,8/)
 elseif (beta .lt. 0.09d0) then
  index=5; mdlbasis(1,1:index) = (/5,6,7,8,9/)
 elseif (beta .lt. 0.22d0 .and. beta .gt. 0.2d0) then
  index=5; mdlbasis(1,1:index) = (/5,6,7,8,10/);
 elseif (beta .lt. 0.8d0 .and. beta .gt. 0.7d0) then
  index=6; mdlbasis(1,1:index) = (/5,6,7,8,9,10/);
 elseif (beta .gt. 2.2d0 .and. beta .lt. 2.3d0) then
  index=7; mdlbasis(1,1:index) = (/5,6,7,8,10,12,15/);
 elseif (beta .lt. 3.2d0 .and. beta .gt. 3.0d0) then
  index=7; mdlbasis(1,1:index) = (/5,6,7,9,10,13,16/);
 end if

case(6)
! Boron
 if (beta .lt. 0.05d0) then
   index=4; mdlbasis(1,1:index) = (/6,7,8,9/)
 elseif (beta .lt. 0.15d0) then
  index=5; mdlbasis(1,1:index) = (/6,7,8,9,10/)
 elseif (beta .lt. 1.2d0 .and. beta .gt. 1.1d0) then
  index=6; mdlbasis(1,1:index) = (/6,7,8,10,12,14/);
 elseif (beta .lt. 0.4d0 .and. beta .gt. 0.3d0) then 
  index=5; mdlbasis(1,1:index) = (/6,7,8,9,11/);
 else if (beta .lt. 3.4d0 .and. beta .gt. 3.3d0) then
!  index=7; mdlbasis(1,1:index) = (/6,7,8,9,12,14,16/)
  index=7; mdlbasis(1,1:index) = (/6,7,8,10,12,13,16/)
 elseif (beta .lt. 4.8d0 .and. beta .gt. 4.7d0) then
  index=7; mdlbasis(1,1:index) = (/6,8,9,11,12,14,18/)
 end if


case(7) 
! Carbon
 if (beta .lt. 0.05d0) then
  index=4; mdlbasis(1,1:index) = (/7,8,9,10/)
 elseif (beta .lt. 0.17d0) then
  index=5; mdlbasis(1,1:index) = (/7,8,9,10,11/)
 elseif (beta .lt. 0.6d0 .and. beta .gt. 0.5d0) then
  index=6; mdlbasis(1,1:index) = (/7,8,9,10,12,13/);
 elseif (beta .lt. 1.9d0 .and. beta .gt. 1.8d0) then
   index=6;  mdlbasis(1,1:index) = (/7,8,9,12,14,15/)
 elseif (beta .lt. 4.6d0 .and. beta .gt. 4.4d0) then
   index=7; mdlbasis(1,1:index) = (/8,9,10,14,18,19,23/)
 elseif (beta .lt. 7.9d0 .and. beta .gt. 7.8d0) then
  index=7; mdlbasis(1,1:index) = (/9,11,13,15,20,22,25/)
 end if

case(8)
! Nitrogen
 if (beta .lt. 0.08d0) then
  index=5;  mdlbasis(1,1:index) = (/8,9,10,11,12/) 
 elseif (beta .lt. 0.22d0) then
  index=5;  mdlbasis(1,1:index) = (/8,9,10,11,13/) 
 elseif (beta .lt. 0.8d0 .and. beta .gt. 0.7d0) then
  index=6; mdlbasis(1,1:index) = (/8,9,10,11,16,18/)
 elseif (beta .lt. 2.8d0 .and. beta .gt. 2.7d0) then;
  index=7; mdlbasis(1,1:index) =(/8,9,11,12,14,16,20/)
 elseif (beta .lt. 6.0d0 .and. beta .gt. 5.8d0) then
  index=7; mdlbasis(1,1:index) = (/9,11,13,16,18,20,25/)
 elseif (beta .lt. 11.7d0 .and. beta .gt. 11.6d0) then
  index=9; mdlbasis(1,1:index) = (/11,12,16,18,22,23,25,26,27/)
 end if

case(9) 
!Oxygen 
 if (beta .lt. 0.1d0) then
  index=5; mdlbasis(1,1:5) = (/9,10,11,12,13/)
 elseif (beta .lt. 0.3d0) then
  index=5;  mdlbasis(1,1:index) = (/9,10,11,12,14/) 
 elseif (beta .lt. 1.1d0 .and. beta .gt. 1.0d0) then
  index=6; mdlbasis(1,1:index) = (/9,10,11,13,14,15/)
 elseif (beta .lt. 3.6d0 .and. beta .gt. 3.5d0) then
   index=7; mdlbasis(1,1:index) = (/10,11,12,14,19,21,25/)
 elseif (beta .lt. 7.6d0 .and. beta .gt. 7.5d0) then
  index=7; mdlbasis(1,1:index) = (/10,11,14,18,20,23,26/)
 elseif (beta .lt. 15.6d0 .and. beta .gt. 15.5d0) then
  index=8; mdlbasis(1,1:index) = (/15,19,22,24,25,29,36,43/)
 end if

case(10)
!Flourine
 if (beta .lt. 0.11d0) then
  index=5; mdlbasis(1,1:index) = (/10,11,12,13,14/)
 elseif (beta .lt. 0.37d0) then
  index=5; mdlbasis(1,1:index) = (/10,11,12,13,15/)
 elseif (beta .gt. 1.3d0 .and. beta .lt. 1.4d0) then
   index =6; mdlbasis(1,1:index) = (/10,11,12,14,17,18/)
 elseif (beta .lt. 4.9d0 .and. beta .gt. 4.8d0) then
  index=7; mdlbasis(1,1:index) = (/10,12,14,17,19,23,28/)
 elseif (beta .gt. 9.2d0 .and. beta .lt. 9.3d0) then
  index=8; mdlbasis(1,1:index) = (/10,13,17,18,22,26,27,29/)
 elseif (beta .gt. 20d0 .and. beta .lt. 21d0) then
  index=9; mdlbasis(1,1:index) = (/17,18,20,21,24,28,29,34,39/)
 end if
 
case(11)
!Neon
 if (beta .lt. 0.14d0) then
  index=5; mdlbasis(1,1:5) = (/11,12,13,14,16/)
 elseif (beta .lt. 0.45d0) then
  index=5; mdlbasis(1,1:index) = (/11,12,14,16, 18/)
 elseif (beta .lt. 1.7d0 .and. beta .gt. 1.6d0) then
  index=6; mdlbasis(1,1:index) = (/ 11,12,13,16,17,19/);
 elseif (beta .gt. 6.1d0 .and. beta .lt. 6.2d0) then
  index=7; mdlbasis(1,1:index) = (/12,13,16,18,19,20,23/)
 elseif (beta .gt. 11.1d0 .and. beta .lt. 11.2d0) then
  index=9; mdlbasis(1,1:index) = (/14,15,18,21,23,24,26,28,34/);
 elseif (beta .gt. 26d0 .and. beta .lt. 27d0) then
  index=8; mdlbasis(1,1:index) = (/19,20,23,24,26,32,39,47/)
 end if

end select

mdlbasis(2,1:index) = 1;
mdlbasis(3,1:index) = 0;

if (index == 0) then; print *, "actn = ", actn, "beta = ", beta;
 print *, "Model basis set is not specified"; call exit; end if;
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     FINDING COEFFICIENTS OF MODEL BASIS FUNCTIONS
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************

   call calcconcmodelsnew(shellmodelRGmat(i:i,1:3,1: lendmdlconc),lenshellmodel(i:i),&
        shellmodelids(i:i,1:4), RGin(1:30,i:i), lendmdlconc,&
        mdlbasis(1:3,1:index),index)

end do
end subroutine concmdlctrl

!!!****************************************************************************************





!!!****************************************************************************************

subroutine calcconcmodelsnew(shellmodelRGmat,lenshellmodel,shellmodelids, RGin, lendmdlconc,&
     mdlbasis, index)
  implicit none
  real*8, intent(out)      :: shellmodelRGmat(1:1,1:3,1: lendmdlconc), shellmodelids(1:1,1:4)
  integer, intent(out)     :: lenshellmodel(1:1)
  integer, intent(in)      ::  lendmdlconc,  mdlbasis(1:3,1:index), index
  real*8, intent(in)       :: RGin(1:30,1:1)
  real*8                   :: pi, intRampGaus4pirK, Nhat, gpf, IR2,nc, antiColRRRs
  real*8                  :: qpMatnew(1:index+1,1:index+1,1:1),qpVec(1:index+1,1:1), qpDummyArray(1:index+1), Rs(1:1), RGr(1:30,1:1)
  real*8                  :: qpMat(1:index+1,1:index+1,1:1),qpDummyMat(1:index+1,1:index+1),qpDummyArray2(1:index+1),FERR, BERR, conditionno
  real*8                  :: WORK(1:10*(index+1)),qpVecnew(1:index+1,1:1),qpVecout(1:index+1,1:1),qpgelssVec(1:index+1,1:1), beta
  integer                  :: DummyIntArray(1:index+1), IWORK(1:index+1), rank, cnt, swit, nmdl, n(1:1)
  integer                  :: K,ntot(1:1), i, j, info, n1(1:1), K1, L1, n2(1:1), K2, L2, rn, m, L, kcurrent, warnswit(1:10)
  character(len=1)         :: outputchar1
  integer    ::  printRRRs
 
  lenshellmodel = 0d0; shellmodelids = 0d0; shellmodelRGmat = 0d0;
  qpgelssVec = 0d0; warnswit = 0;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This calculates the coefficients of the model for Ss in terms of a sum of S ramps of different degrees
  !
  ! It also does Sp -> P ramps (combined with the shelltobf function) for concentric shell pairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  
  pi = dacos(-1d0);
  n(1:1) = RGin(7,1); 
  beta = Rgin(8,1);
  nc = Rgin(9,1);
  RGr = RGin;
     
     !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     MATRIX CONSTRUCTION
! Calculate overlap of model basis functions with each other 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  qpMat = 0d0;
  do i=1, index; 
   n1 = mdlbasis(1,i);   K1 = mdlbasis(2,i);   L1 = mdlbasis(3,i);
     do j=1,index; 
        n2 = mdlbasis(1,j);   K2 = mdlbasis(2,j);   L2 = mdlbasis(3,j);
        if (K1 == K2) then;
  qpMat(i,j,1:1) = (8*Pi)/((1d0+n1)*(2d0+n1)*(3d0+n1)*(4+n1)*(5+n1))*(&
     (6d0*(4+n1)*(5+n1))/((1d0+n2)*(2d0+n2)*(3d0+n2)*(4d0+n2))&
     -(30d0+6d0*n1+4*n2)/((5d0+n1+n2)*(6d0+n1+n2)*(7d0+n1+n2))&
     +4d0/(2d0+3d0*n2+n2**2))
        end if
     end do; 
  end do; 


do i=1,index;
    n1 = mdlbasis(1,i)
     qpMat(i,index+1,1:1) = 2d0/((n1+1d0)*(n1+2)*(n1+3))*(2d0*sqrt(pi))
     qpMat(index+1,i,1:1) = 2d0/((n1+1d0)*(n1+2)*(n1+3))*(2d0*sqrt(pi))
end do  
qpMat(index+1,index+1,1) = 0d0;


!print *, "MATRIX"
do i=1,index+1;
 do j=1,index+1;
!   print *, i,j, qpMat(i,j,1)
 end do
end do
!print *, ""

!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     VECTOR CONSTRUCTION
! Calculate overlap of model basis functions with gaussian 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************

  

  do i=1, index; 
     nmdl = mdlbasis(1,i)
     call calcintRampGaus4pirK(intRampGaus4pirK ,int(n+nmdl),beta,2, 40,2)
     call calcvecantiCol(antiColRRRs, int(nmdl), int(n(1)), beta, 40, 2)
     qpVec(i,1) = antiColRRRs/(4d0*pi)
  end do
  call calcintRampGaus4pirK(intRampGaus4pirK ,n,beta,2, 40,2)
  qpVec(index+1,1) = intRampGaus4pirK/(Nhat(int(n))*gpf(beta))/(2*sqrt(dacos(-1d0)))

!********!!!!!!!!!!!!************!!!!!!!!!!!!**************
!     LINEAR SOLVE
! Calculate coefficients by solving Mc = v 
!********!!!!!!!!!!!!************!!!!!!!!!!!!**************


           qpVecnew = qpVec; qpDummyMat = 0d0;  qpMatnew = qpMat;  qpDummyArray = 0d0;

   call dgelss(index+1,index+1,1,qpMatnew(1:index+1,1:index+1,1), index+1, qpVecnew(1:index+1,1),index+1,&
                qpDummyMat(1:index+1,1:index+1),1d-16, rank, work(1:10*(index+1)), 10*(index+1), info) 

qpgelssVec = qpVecNew;

call DGEMV ( 'N', index+1, index+1, 1d0, qpMat(1:index+1,1:index+1,1), index+1, &
              qpVecnew(1:index+1,1), 1,-1d0, qpVec(1:index+1,1),1)


           if (info .ne. 0) then; call exit; end if
        cnt = 0;   
        do k=1,index ;   
           if (Abs(qpgelssVec(k,1)) > 1d-10) then
              cnt = cnt + 1;
              shellmodelRgmat(1,1,cnt) = dble(qpgelssVec(k,1))*2*sqrt(pi);  
              shellmodelRGmat(1,2,cnt) = mdlbasis(1,k); ! Degree of model ramp
              shellmodelRGmat(1,3,cnt) = mdlbasis(2,k);      ! Angular momentum of model ramp in K form
           end if;
        end do
        lenshellmodel(1) = cnt;
  
  shellmodelids(1:1,1) = RGr(1,1:1);   shellmodelids(1:1,2) = RGr(2,1:1);
  shellmodelids(1:1,3) = RGr(3,1:1);   shellmodelids(1:1,4) = RGr(9,1:1);


  return
end subroutine calcconcmodelsnew







!******************

subroutine calcvecantiCol(antiColRRRs, nmdl, n, beta, pts, type)
 implicit none
 integer, intent(in)   :: n, nmdl, pts, type
 real*8, intent(in)    :: beta
 real*8  :: antiColRRRs, pi, absc(1:2*pts+1), weights(1:2*pts+1), integrand(1:2*pts+1), r
 integer :: v

 pi = dacos(-1d0);

  call getquad(weights,absc,pts,type); 

  do v=1,2*pts+1; 
       r = absc(v);
       integrand(v) = exp(-beta*r**2)*(32*Pi**2*(1-r)**n*r*(4+20*r**2+9*nmdl*r**2+nmdl**2*r**2-(1-r)**(4+nmdl)*(4+r+nmdl*r)))&
                          /((1d0+nmdl)*(2d0+nmdl)*(3d0+nmdl)*(4d0+nmdl)*(5d0+nmdl))*weights(v)
  end do

  antiColRRRs = sum(integrand); 

end subroutine calcvecantiCol




!***********************

subroutine checkRsss(shellmodelRGmat,cnt,n,beta,eta,pts,type,nc,printRRRs)
 implicit none
 integer, intent(in)   :: n,cnt
 real*8, intent(in)    :: beta, shellmodelRGmat(1:1,1:3,1:cnt),eta
 real*8  :: actRsss, modelRsss, pi,mdlcoeff,nc,nhat
 integer :: i, pts, type, nmdl, v,printRRRs
 real*8  :: absc(1:2*pts+1), weights(1:2*pts+1), pf 
 real*8  :: integrand(1:2*pts+1), otherpart(1:2*pts+1), gpf

  call getquad(weights,absc,pts,type); 

pi = dacos(-1d0)
 
 modelRsss = 0d0;
 do i=1,cnt
    nmdl = shellmodelRGmat(1,2,i)
    mdlcoeff = shellmodelRGmat(1,1,i)
  do v=1,2*pts+1; 
       integrand(v) = absc(v)*weights(v)*erf(sqrt(eta)*absc(v))*(1-absc(v))**nmdl;  
  end do
    modelRsss = modelRsss + mdlcoeff*sum(integrand)
 end do

modelRsss=modelRsss/(2*sqrt(pi))

pf = 4*pi**(5d0/2d0)/(eta**(3d0/2d0))

  do v=1,2*pts+1; 
       integrand(v) = absc(v)*weights(v)*erf(sqrt(eta)*absc(v))*(1-absc(v))**n*exp(-beta*absc(v)**2);
  end do
  actRsss = sum(integrand); 
nc = nhat(n)*gpf(beta)*gpf(eta)
if (abs(nc*modelRsss-nc*actRsss) .gt. 1d-9) then;
   printRRRs = 0;
end if
if (abs(nc*modelRsss-nc*actRsss) .gt. 1d-10) then;
   print "(a15,i3,2f15.8,e10.1,f8.2,30i3)", "MDLPROB in Rsss",n,beta,eta,nc*(modelRsss-actRsss),&
       nc*(modelRsss-actRsss)*10**8, int(shellmodelRGmat(1,2,1:cnt));
 end if
if (abs(modelRsss-actRsss) .gt. 1d-12) then;
  print *, "n=",n, "beta = ", beta, "eta = ", eta
  print *,"modelRsss = ", modelRsss, "actRsss = ", actRsss, "dif = ", modelRsss-actRsss
end if

end subroutine checkRsss








!***********************

subroutine checkRRRs(shellmodelRGmat,cnt,n,beta,pts,type,nc)
 implicit none
 integer, intent(in)   :: n,cnt
 real*8, intent(in)    :: beta, shellmodelRGmat(1:1,1:3,1:cnt)
 real*8  :: actRRRs, modelRRRs, pi,mdlcoeff,nc,intpart
 integer :: i, pts, type, nmdl, v
 real*8  :: absc(1:2*pts+1), weights(1:2*pts+1), pf, nhat, integrand(1:2*pts+1), otherpart(1:2*pts+1), gpf

  call getquad(weights,absc,pts,type); 

pi = dacos(-1d0)
 
 modelRRRs = 0d0;
 do i=1,cnt
    nmdl = shellmodelRGmat(1,2,i)
    mdlcoeff = shellmodelRGmat(1,1,i)
  do v=1,2*pts+1;
       integrand(v) = 8*pi**2*absc(v)*weights(v)*(2d0-(-1+absc(v))**(2+2*n)*(2+absc(v)+2*n*absc(v)))*&
            (1-absc(v))**nmdl/((1d0+n)*(1d0+2*n)*(3d0+2*n));
  end do

    modelRRRs = modelRRRs + mdlcoeff*sum(integrand)
 end do

modelRRRs=modelRRRs/(2*sqrt(pi))

  do v=1,2*pts+1; 
       integrand(v) = 8*pi**2*absc(v)*weights(v)*(2d0-(-1+absc(v))**(2+2*n)*(2+absc(v)+2*n*absc(v)))*&  
           (1-absc(v))**n*exp(-beta*absc(v)**2)/((1d0+n)*(1d0+2*n)*(3d0+2*n));
  end do
  actRRRs = sum(integrand); 
nc = nhat(n)*gpf(beta)

!if (abs(nhat(n)**2*nc*modelRRRs-nhat(n)**2*nc*actRRRs) .gt. 1d-9) then;
   print "(a15,i3,f15.8,e10.1,f10.5,30i3)", "MDLPROB in RRRs",n,beta,nhat(n)**2*nc*(modelRRRs-actRRRs),&
       nhat(n)**2*nc*(modelRRRs-actRRRs)*10**8, int(shellmodelRGmat(1,2,1:cnt)); 

!end if

if (abs(modelRRRs-actRRRs) .gt. 1d-12) then;
  print *, "n=",n, "beta = ", beta
  print *,"modelRRRs = ", modelRRRs, "actRRRs = ", actRRRs, "dif = ", modelRRRs-actRRRs
end if

end subroutine checkRRRs

!***********************************************
