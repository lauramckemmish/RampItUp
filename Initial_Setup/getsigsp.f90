subroutine getsigspRx(algchange,noRR, noRGSs, noRGSp, RR, RgSs, &
     RgSp, nobasisfun,noatoms,atoms,basis,sigspprint,thresh)
  implicit none
  integer, intent(in)    :: nobasisfun, noatoms,thresh 
  integer, intent(out)   :: noRR, noRGSs, noRGSp
  integer, intent(inout) :: algchange(1:30)
  real*8, intent(in)     ::  basis(1:nobasisfun,1:30), atoms(1:5,1:noatoms)
  real*8, intent(out)    :: RR(1:30,1:4*nobasisfun), RgSs(1:30,1:20*nobasisfun), RgSp(1:30,1:20*nobasisfun)
  integer                :: bf1, bf2, sigspprint,  i, j, k, k2, noS(0:11), idS(0:11,1:nobasisfun), x
  real*8                 :: nc1, nc2,  n, Ax, Ay, Az, Bx, By, Bz, BAx, BAy, BAz, beta, alpha, RAB2, zeta
  real*8 :: gpf, nhat
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program constructs a list of significant SS, Ss and Sp shell pairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      RR contains info about the SS shell pair and has elements (shellbf1, shellbf2, atoms centre,ramp deg1, rampdeg2, nc1*n2, Ax, Ay, Az, blank)
  !     
  !      RGSs contains info about the itive Ss shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
  !      RGSp contains info about the itive Sp shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  noRR   = 0;  noRGSs = 0;  noRGSp = 0;  algchange = 0;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! algchange(1) is no of (Ss| shells concentric with large beta
!!! algchange(2) is tot no of (Sp| shells concentric 
!!! algchange(3) is no of (ss| which are concentric with large zeta
!!! algchange(4) is tot no of concentric (Ss| shell pairs (large and small beta)
!!! algchange(5) is no of (sp| which are concentric with large zeta
!!! algchange(6) is no of (pp| which are concentric with large zeta
!!! algchange(7) is no of concentric (Sp| shell pairs with large beta
!!! algchange(8) is
!!! algchange(9) is
!!! algchange(10) is
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  noS = 0;
  idS = 0;
  do i=1,nobasisfun
     noS(int(basis(i,3))) = noS(int(basis(i,3))) +1 
     idS(int(basis(i,3)), noS(int(basis(i,3)))) = i   
  end do
  
  !******************************************************
  !*** Constructing |SS) shell pairs  *******************
  !******************************************************
  
  ! (S|S) where both are only S shells (not Ss shells)
  
  do i=1,noS(0)
     bf1 = idS(0,i)
     do j=1,noS(0)
        bf2 = idS(0,j)
        if (basis(bf1,1) == basis(bf2,1)) then
           noRR = noRR + 1;
           RR(1:2,noRR) = (/bf1,bf2/)
           RR(3,noRR) = basis(bf1,1)   ! atoms center
           RR(4,noRR) = basis(bf1,11)   ! ramp deg1
           RR(5,noRR) = basis(bf2,11)   ! ramp deg2
           RR(6,noRR) = basis(bf1,21)*basis(bf2,21)   ! nc1*nc2
           RR(7:9,noRR) = atoms(3,int(basis(bf1,1)))  
        end if
     end do
  end do
  
  ! (S|S) where one is a Ss shell and one is a S shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     do j=1,noS(0)
        bf2 = idS(0,j)
        if (basis(bf1,1) == basis(bf2,1)) then
           noRR = noRR + 1;
           RR(1:2,noRR) = (/bf1,bf2/)
           RR(3,noRR) = basis(bf1,1)   ! atoms center
           RR(4,noRR) = basis(bf1,11)   ! ramp deg1
           RR(5,noRR) = basis(bf2,11)   ! ramp deg2
           RR(6,noRR) = basis(bf1,21)*basis(bf2,21)   ! nc1*nc2
           RR(7:9,noRR) = atoms(3,int(basis(bf1,1)))   
        end if
     end do
  end do
  
  ! (S|S) where one is a Ss shell and one is a Ss shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     do j=1,noS(5)
        bf2 = idS(5,j)
        if (basis(bf1,1) == basis(bf2,1)) then
           noRR = noRR + 1;
           RR(1:2,noRR) = (/bf1,bf2/)
           RR(3,noRR) = basis(bf1,1)   ! atoms center
           RR(4,noRR) = basis(bf1,11)   ! ramp deg1
           RR(5,noRR) = basis(bf2,11)   ! ramp deg2
           RR(6,noRR) = basis(bf1,21)*basis(bf2,21)   ! nc1*nc2
           RR(7:9,noRR) = atoms(3,int(basis(bf1,1))) 
        end if
     end do
  end do
  
  
  !******************************************************
  !*** Constructing |Ss) shell pairs  *******************
  !******************************************************
  
  noRGSs = 0;
 
  algchange(1) = noRGSs
  
  
  ! (S|s) where one is S shell and one is a s shell
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(1)
        bf2 = idS(1,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then 
              if (BAx ==0 .and. BAy == 0 .and. BAz ==0) then
                 noRGSs = noRGSs + 1;
                 RgSs(1:2,noRGSs) = (/bf1,bf2/)
                 RgSs(3,noRGSs) = basis(bf1,1)
                 RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
                 RgSs(7:9,noRGSs) = (/n,beta,nc1*nc2/)
                 RgSs(10,noRGSs) = 1
                 RgSs(11:13,noRgSs) = atoms(3:5,int(basis(bf1,1))); !Ax,Ay,Az
              end if
           end if
        end do
     end do
  end do
  
  ! (S|s) where first is S shell and second is a Ss shell
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(5)
        bf2 = idS(5,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        beta = basis(bf2,12)
        nc2 = basis(bf2,22) 
        if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then 
           if (BAx ==0 .and. BAy == 0 .and. BAz ==0) then
              noRGSs = noRGSs + 1;
                 RgSs(1:2,noRGSs) = (/bf1,bf2/)
                 RgSs(3,noRGSs) = basis(bf1,1)
                 RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
                 RgSs(7:9,noRGSs) = (/n,beta,nc1*nc2/)
                 RgSs(10,noRGSs) = 1
                 RgSs(11:13,noRgSs) = atoms(3:5,int(basis(bf1,1))); !Ax,Ay,Az
           end if
        end if
     end do
  end do
  
  
  ! (S|s) where first is Ss shell and second is a s shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(1)
        bf2 = idS(1,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then 
              if (BAx ==0 .and. BAy == 0 .and. BAz ==0) then
                 noRGSs = noRGSs + 1;
                 RgSs(1:2,noRGSs) = (/bf1,bf2/)
                 RgSs(3,noRGSs) = basis(bf1,1)
                 RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
                 RgSs(7:9,noRGSs) = (/n,beta,nc1*nc2/)
                 RgSs(10,noRGSs) = 1
                 RgSs(11:13,noRgSs) = atoms(3:5,int(basis(bf1,1))); !Ax,Ay,Az
              end if
           end if
        end do
     end do
  end do
  
  ! (S|s) where first is Ss shell and second is a Ss shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(5)
        bf2 = idS(5,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        beta = basis(bf2,12)
        nc2 = basis(bf2,22) 
        if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then 
           if (BAx ==0 .and. BAy == 0 .and. BAz ==0) then
              noRGSs = noRGSs + 1;
              RgSs(1:2,noRGSs) = (/bf1,bf2/)
              RgSs(3  ,noRgSs) = basis(bf1,1)
              RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)));
              RgSs(7:8,noRGSs) = (/n,beta/)
if (i == j) then;
              RgSs(9,noRGSs) = 2*nc1*nc2
else
              RgSs(9,noRGSs) = nc1*nc2
end if
              RgSs(10,noRGSs) = 1
              RgSs(11:13,noRgSs) = (/Ax,Ay,Az/);
           end if
        end if
     end do
  end do
  
  algchange(4) = noRGSs
  
  ! Algchange(4) counts how many concentric Ss shell pairs.
  
  
  ! (S|s) where one is S shell and one is a s shell
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(1)
        bf2 = idS(1,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35 .and.((BAx .ne. 0) .or. (BAy .ne. 0) .or. (BAz .ne. 0))) then
              noRGSs = noRGSs + 1;
                 RgSs(1:2,noRGSs) = (/bf1,bf2/)
                 RgSs(3,noRGSs) = basis(bf1,1)
                 RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
                 RgSs(7:9,noRGSs) = (/n,beta,nc1*nc2/)
                 RgSs(10,noRGSs) = 1
                 RgSs(11:13,noRgSs) = atoms(3:5,int(basis(bf1,1))); !Ax,Ay,Az
           end if
        end do
     end do
  end do
  
  
  
  ! (S|s) where first is S shell and second is a Ss shell
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(5)
        bf2 = idS(5,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        beta = basis(bf2,12)
        nc2 = basis(bf2,22) 
        if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then 
           if ((BAx .ne. 0) .or. (BAy .ne. 0) .or. (BAz .ne. 0)) then

              noRGSs = noRGSs + 1;
                 RgSs(1:2,noRGSs) = (/bf1,bf2/)
                 RgSs(3,noRGSs) = basis(bf1,1)
                 RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
                 RgSs(7:9,noRGSs) = (/n,beta,nc1*nc2/)
                 RgSs(10,noRGSs) = 1
                 RgSs(11:13,noRgSs) = atoms(3:5,int(basis(bf1,1))); !Ax,Ay,Az
           end if
        end if
     end do
  end do
  
  
  ! (S|s) where first is Ss shell and second is a s shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(1)
        bf2 = idS(1,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35 .and. ((BAx .ne. 0) .or. (BAy .ne. 0) .or. (BAz .ne. 0))) then

              noRGSs = noRGSs + 1;
                 RgSs(1:2,noRGSs) = (/bf1,bf2/)
                 RgSs(3,noRGSs) = basis(bf1,1)
                 RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
                 RgSs(7:9,noRGSs) = (/n,beta,nc1*nc2/)
                 RgSs(10,noRGSs) = 1
                 RgSs(11:13,noRgSs) = atoms(3:5,int(basis(bf1,1))); !Ax,Ay,Az
           end if
        end do
     end do
  end do

x = noRgSs;
  
  ! (S|s) where first is Ss shell and second is a Ss shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(5)
        bf2 = idS(5,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        beta = basis(bf2,12)
        nc2 = basis(bf2,22) 
        if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then 
           if ((BAx .ne. 0) .or. (BAy .ne. 0) .or. (BAz .ne. 0)) then
  print *, "IMPORTANT - haven't figured out how to deal with Ss shell pairs formed from &
       two non-concentric Ss shells in NA and multipole integrals"
print *, "IMPORTANT - partially neglected shell pair has n = ", int(n), "beta = ", beta, &
   "unnorm nc1*nc2 = ", nc1*nc2, "beta*(BAx**2+BAy**2+BAz**2) = ", beta*(BAx**2+BAy**2+BAz**2);
              noRGSs = noRGSs + 1;
              RgSs(1:2,noRGSs) = (/bf1,bf2/)
              RgSs(3  ,noRgSs) = basis(bf1,1)
              RgSs(4:6,noRGSs) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)));
              RgSs(7:8,noRGSs) = (/n,beta/)
if (i == j) then;
              RgSs(9,noRGSs) = 2*nc1*nc2
else
              RgSs(9,noRGSs) = nc1*nc2
end if
              RgSs(10,noRGSs) = 1
              RgSs(11:13,noRgSs) = (/Ax,Ay,Az/);
              RgSs(14,noRgSs) = 111;
           end if
        end if
     end do
  end do

if ((noRgSs-x) == 1) then;
  print *, "IMPORTANT - shell pair not neglected after all -  n = ", int(RgSs(7,noRgSs)),&
        "beta = ", RgSs(8,noRgSs), "bf1,2 = ", int(RgSs(1:2,noRgSs))
  RgSs(14,noRgSs) = 0;
end if;
  
  !******************************************************
  !*** Constructing |Sp) shell pairs  *******************
  !******************************************************
  
  noRGSp = 0;
  ! (S|p) where first is a S shell and second is a p shell

  algchange(7)= noRGSp
  
  
  ! (S|p) where first is a S shell and second is a p shell
  
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(2)
        bf2 = idS(2,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then
              if  (BAx ==0 .and. BAy == 0 .and. BAz == 0) then
                 noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 1
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
              end if
           end if
        end do
     end do
  end do
  
  
  
  ! (S|p) where first is a Ss shell and second is a p shell
  
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(2)
        bf2 = idS(2,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then
              if  (BAx ==0 .and. BAy == 0 .and. BAz == 0) then
                 noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 1
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
              end if
           end if
        end do
     end do
  end do
  


  
  ! (S|p) where first is a S shell and second is a Pp shell
  
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(9)
        bf2 = idS(9,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=2,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then
              if  (BAx ==0 .and. BAy == 0 .and. BAz == 0) then
                 noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 1
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
              end if
           end if
        end do
     end do
  end do
  
  
  
  ! (S|p) where first is a Ss shell and second is a Pp shell
  
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(9)
        bf2 = idS(9,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=2,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35) then
              if  (BAx ==0 .and. BAy == 0 .and. BAz == 0) then
                 noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 1
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
              end if
           end if
        end do
     end do
  end do
  


  algchange(2)= noRGSp
  
  ! (S|p) where first is a S shell and second is a p shell
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(2)
        bf2 = idS(2,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35 .and.&
                ((BAx .ne. 0) .or. (BAy .ne. 0) .or. (BAz .ne. 0))) then
              noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 0
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
           end if
        end do
     end do
  end do
  
  
  ! (S|p) where first is a Ss shell and second is a p shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(2)
        bf2 = idS(2,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=1,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35 .and. &
                (BAx  .ne. 0 .or. BAy .ne. 0 .or. BAz .ne. 0)) then
              noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 1
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
           end if
        end do
     end do
  end do
  
  ! (S|p) where first is a S shell and second is a Pp shell
  do i=1,noS(0)
     bf1 = idS(0,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(9)
        bf2 = idS(9,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=2,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35 .and.&
                ((BAx .ne. 0) .or. (BAy .ne. 0) .or. (BAz .ne. 0))) then
              noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 0
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
           end if
        end do
     end do
  end do
  
  
  ! (S|p) where first is a Ss shell and second is a Pp shell
  do i=1,noS(5)
     bf1 = idS(5,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the ramp 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     n = basis(bf1,11)
     nc1 = basis(bf1,21)
     do j=1,noS(9)
        bf2 = idS(9,j)
        Bx = atoms(3,int(basis(bf2,1)))
        By = atoms(4,int(basis(bf2,1)))
        Bz = atoms(5,int(basis(bf2,1)))
        BAx = Bx - Ax
        BAy = By - Ay
        BAz = Bz - Az
        do k=2,int(basis(bf2,2))
           beta = basis(bf2,k+10)
           nc2 = basis(bf2,k+20) 
           if (beta*(BAx*BAx+BAy*BAy+BAz*BAz)< 35 .and. &
                (BAx  .ne. 0 .or. BAy .ne. 0 .or. BAz .ne. 0)) then
              noRGSp = noRGSp + 1;
              RgSp(1:2,noRGSp) = (/bf1,bf2/)
              RgSp(  3,noRGSp) = basis(bf1,1)
              RgSp(4:6,noRGSp) = atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1)))
              RgSp(7:9,noRGSp) = (/n,beta,nc1*nc2/)
              RgSp(10,noRGSp) = 1
              RgSp(11:13,noRgSp) = (/Ax,Ay,Az/);
           end if
        end do
     end do
  end do

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (sigspprint ==1) then
     
     print *, "NoRR = ", noRR
     print *, "SS info"
     do i=1,noRR
        print "(a2,i3,a2,i3,a4,4f23.15)", "(", int(RR(1,i)), ",", int(RR(2,i)), ") = ", &
             RR(3,i), RR(4,i), RR(5,i),RR(6,i)
     end do
     print *, ""
     
     
     print *, "NoSs = ", noRGSs
     print *, "Ss info"
     do i=1,noRGss
        print "(a2,i3,a2,i3,a4,3f23.15,a20,3f10.5)", "(", int(RgSs(1,i)), ",", int(RgSs(2,i)), ") = ", &
             RgSs(7,i), RgSs(8,i), RgSs(9,i), "BA dist", RgSs(4,i), RgSs(5,i), RgSs(6,i)
     end do
     print *, ""
     
     print *, "NoSp = ", noRGSp
     print *, "Sp info"
     do i=1,noRGsp
        print "(a2,i3,a2,i3,a4,3f23.15,a20,3f10.5)", "(", int(RgSp(1,i)), ",", int(RgSp(2,i)), ") = ", &
             RgSp(7,i), RgSp(8,i), RgSp(9,i), "BA dist", RgSp(4,i), RgSp(5,i), RgSp(6,i)
     end do
     print *, ""

     
  end if
  
  return
  
end subroutine getsigspRx










!***************************************
subroutine getsigspss(algchange, noggss,ggss,  nobasisfun,noatoms,atoms,basis,sigspprint,thresh)
  implicit none
  integer, intent(in)    :: nobasisfun, noatoms , thresh
  integer, intent(out)   :: noggss
  integer, intent(inout) :: algchange(1:30)
  real*8, intent(in)     :: basis(1:nobasisfun,1:30), atoms(1:5,1:noatoms)
  real*8, intent(out)    :: ggss(1:30,1:20*nobasisfun**2)
  integer                :: bf1, bf2, sigspprint,  i, j, k, k2, noS(0:11), idS(0:11,1:nobasisfun)
  real*8                 :: nc1, nc2,  n, Ax, Ay, Az, Bx, By, Bz, BAx, BAy, BAz, beta, alpha, RAB2, zeta, pi, ol
  real*8 :: pi3r2

pi3r2 = dacos(-1d0)**(3d0/2d0)  
print *, "THRESH = ", thresh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program constructs a list of significant SS, Ss and Sp shell pairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      RR contains info about the SS shell pair and has elements (shellbf1, shellbf2, atoms centre,ramp deg1, rampdeg2, nc1*n2, Ax, Ay, Az, blank)
  !     
  !      RGSs contains info about the itive Ss shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
  !      RGSp contains info about the itive Sp shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  pi = dacos(-1d0);
  
noggss = 0;  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! algchange(1) is no of (Ss| shells concentric with large beta
!!! algchange(2) is tot no of (Sp| shells concentric 
!!! algchange(3) is no of (ss| which are concentric with large zeta
!!! algchange(4) is tot no of concentric (Ss| shell pairs (large and small beta)
!!! algchange(5) is no of (sp| which are concentric with large zeta
!!! algchange(6) is no of (pp| which are concentric with large zeta
!!! algchange(7) is no of concentric (Sp| shell pairs with large beta
!!! algchange(8) is
!!! algchange(9) is
!!! algchange(10) is
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  noS = 0;
  idS = 0;
  do i=1,nobasisfun
     noS(int(basis(i,3))) = noS(int(basis(i,3))) +1 
     idS(int(basis(i,3)), noS(int(basis(i,3)))) = i   
  end do
 

  !******************************************************
  !*** Constructing |ss) shell pairs  *******************
  !******************************************************
  
  
  noggss = 0;
  
  
  ! Both s shells
  do i=1,noS(1)
     bf1 = idS(1,i)
     Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
     Ay = atoms(4,int(basis(bf1,1)))
     Az = atoms(5,int(basis(bf1,1)))
     do k=1,int(basis(bf1,2))
        alpha = basis(bf1,10+k)
        nc1 = basis(bf1,20+k)
        do j=i,noS(1)
             bf2 = idS(1,j)
             Bx = atoms(3,int(basis(bf2,1)))
             By = atoms(4,int(basis(bf2,1)))
             Bz = atoms(5,int(basis(bf2,1)))
             do k2=1,int(basis(bf2,2))  
                beta = basis(bf2,10+k2)
                nc2 = basis(bf2,20+k2)  
                zeta = alpha + beta
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
            noggss = noggss+1;
            ggss(1:2,noggss) = (/bf1,bf2/);
            ggss(4,noggss) = (alpha * Ax + beta*Bx)/(zeta)
            ggss(5,noggss) = (alpha * Ay + beta*By)/(zeta)
            ggss(6,noggss) = (alpha * Az + beta*Bz)/(zeta)
            ggss(8,noggss) = zeta
            ggss(9,noggss) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)
            ggss(11:12,noggss) = (/basis(bf1,1),basis(bf2,1)/)
            ggss(13:16,noggss) = (/1d0/(zeta**(3d0/2d0)),alpha,beta,RAB2/);
         end if
            end do
         end do
      end do
   end do
   
! First s shell, Second Ss shell
   do i=1,noS(1)
      bf1 = idS(1,i)
      Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
      Ay = atoms(4,int(basis(bf1,1)))
      Az = atoms(5,int(basis(bf1,1)))
      do k=1,int(basis(bf1,2))
         alpha = basis(bf1,10+k)
         nc1 = basis(bf1,20+k)
         do j=1,noS(5)
            bf2 = idS(5,j)
            Bx = atoms(3,int(basis(bf2,1)))
            By = atoms(4,int(basis(bf2,1)))
            Bz = atoms(5,int(basis(bf2,1)))
            beta = basis(bf2,12)
            nc2 = basis(bf2,22) 
            zeta = alpha + beta
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
        
     if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
            noggss = noggss+1;
            ggss(1:2,noggss) = (/bf1,bf2/);
            ggss(4,noggss) = (alpha * Ax + beta*Bx)/(zeta)
            ggss(5,noggss) = (alpha * Ay + beta*By)/(zeta)
            ggss(6,noggss) = (alpha * Az + beta*Bz)/(zeta)
            ggss(8,noggss) = zeta
            ggss(9,noggss) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)
            ggss(11:12,noggss) = (/basis(bf1,1),basis(bf2,1)/)
            ggss(13:16,noggss) = (/1d0/(zeta**(3d0/2d0)),alpha,beta,RAB2/);
          end if
         end do
      end do
   end do
   
   ! First Ss shell, Second Ss shell
   do i=1,noS(5)
      bf1 = idS(5,i)
      Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
      Ay = atoms(4,int(basis(bf1,1)))
      Az = atoms(5,int(basis(bf1,1)))
      alpha = basis(bf1,12)
      nc1 = basis(bf1,22)
      do j=i,noS(5)
         bf2 = idS(5,j)
         Bx = atoms(3,int(basis(bf2,1)))
         By = atoms(4,int(basis(bf2,1)))
         Bz = atoms(5,int(basis(bf2,1)))
         beta = basis(bf2,12)
         nc2 = basis(bf2,22) 
         zeta = alpha + beta
         RAB2 = (Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2;
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
           noggss = noggss+1;
            ggss(1:2,noggss) = (/bf1,bf2/);
            ggss(4,noggss) = (alpha * Ax + beta*Bx)/(zeta)
            ggss(5,noggss) = (alpha * Ay + beta*By)/(zeta)
            ggss(6,noggss) = (alpha * Az + beta*Bz)/(zeta)
            ggss(8,noggss) = zeta
            ggss(9,noggss) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)
            ggss(11:12,noggss) = (/basis(bf1,1),basis(bf2,1)/)
            ggss(13:16,noggss) = (/1d0/(zeta**(3d0/2d0)),alpha,beta,RAB2/);
         end if
      end do
   end do
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (sigspprint ==1) then
     

     print *, "Noss = ", noggss
     print *, "ss info"
     do i=1,noggss
        print "(a2,i3,a2,i3,a4,2f23.15,a12,4f10.5,f10.5)", "(", int(ggss(1,i)), ",", int(ggss(2,i)), ") = ", ggss(8,i), ggss(9,i), &
             "centered", ggss(4,i), ggss(5,i), ggss(6,i), sqrt(ggss(4,i)**2+ggss(5,i)**2+ggss(6,i)**2), ggss(16,i)
     end do
     print *, ""
     

     
  end if
  
  return
  
end subroutine getsigspss







!***************************************
subroutine getsigspsp(algchange,noggsp, ggsp, nobasisfun,noatoms,atoms,basis,sigspprint,thresh)
  implicit none
  integer, intent(in)    :: nobasisfun, noatoms ,thresh
  integer, intent(out)   :: noggsp
  integer, intent(inout) :: algchange(1:30)
  real*8, intent(in)     ::  basis(1:nobasisfun,1:30), atoms(1:5,1:noatoms)
  real*8, intent(out)    :: ggsp(1:30,1:20*nobasisfun**2)
  integer                :: bf1, bf2, sigspprint,  i, j, k, k2, noS(0:11), idS(0:11,1:nobasisfun)
  real*8                 :: nc1, nc2,  n, Ax, Ay, Az, Bx, By, Bz, BAx, BAy, BAz, beta, alpha, RAB2, zeta, pi
  real*8 :: pi3r2

pi3r2 = dacos(-1d0)**(3d0/2d0)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program constructs a list of significant SS, Ss and Sp shell pairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      RR contains info about the SS shell pair and has elements (shellbf1, shellbf2, atoms centre,ramp deg1, rampdeg2, nc1*n2, Ax, Ay, Az, blank)
  !     
  !      RGSs contains info about the itive Ss shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
  !      RGSp contains info about the itive Sp shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  pi = dacos(-1d0);
  
 noggsp = 0;
!   ggsp   = 0d0; 
print *, "THRESH = ", thresh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! algchange(1) is no of (Ss| shells concentric with large beta
!!! algchange(2) is tot no of (Sp| shells concentric 
!!! algchange(3) is no of (ss| which are concentric with large zeta
!!! algchange(4) is tot no of concentric (Ss| shell pairs (large and small beta)
!!! algchange(5) is no of (sp| which are concentric with large zeta
!!! algchange(6) is no of (pp| which are concentric with large zeta
!!! algchange(7) is no of concentric (Sp| shell pairs with large beta
!!! algchange(8) is
!!! algchange(9) is
!!! algchange(10) is
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  noS = 0;
  idS = 0;
  do i=1,nobasisfun
     noS(int(basis(i,3))) = noS(int(basis(i,3))) +1 
     idS(int(basis(i,3)), noS(int(basis(i,3)))) = i   
  end do
 

   
   !******************************************************
   !*** Constructing |sp) shell pairs  *******************
   !******************************************************
   
   
   noggsp = 0;
    
    ! First s shells, second p shell
    do i=1,noS(1)
       bf1 = idS(1,i)
       Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
       Ay = atoms(4,int(basis(bf1,1)))
       Az = atoms(5,int(basis(bf1,1)))
       do k=1,int(basis(bf1,2))
          alpha = basis(bf1,10+k)
          nc1 = basis(bf1,20+k)
          do j=1,noS(2)
             bf2 = idS(2,j)
             Bx = atoms(3,int(basis(bf2,1)))
             By = atoms(4,int(basis(bf2,1)))
             Bz = atoms(5,int(basis(bf2,1)))
             do k2=1,int(basis(bf2,2))  
               beta = basis(bf2,10+k2)
               nc2 = basis(bf2,20+k2)
                 zeta = alpha + beta
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
                  noggsp = noggsp+1;
                ggsp(1:2,noggsp) = (/bf1,bf2/)
                ggsp(3,noggsp) = (alpha * Ax + beta*Bx)/(zeta);
                ggsp(4,noggsp) = (alpha * Ay + beta*By)/(zeta);
                ggsp(5,noggsp) = (alpha * Az + beta*Bz)/(zeta);
                ggsp(6:8,noggsp) = (/alpha,beta,zeta/);
                ggsp(9,noggsp) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)
                ggsp(11:13,noggsp) = atoms(3:5,int(basis(bf1,1)))-atoms(3:5,int(basis(bf2,1)));
                ggsp(14:16,noggsp) = (/sqrt(zeta),  1d0/sqrt(zeta),1d0/zeta/);
                ggsp(17:18,noggsp) = (/basis(bf1,1),basis(bf2,1)/);
                ggsp(19,noggsp) = 1d0/(zeta)**(5d0/2d0)
              end if
             end do
          end do
       end do
    end do

    ! First Ss shell, second p shell
    do i=1,noS(5)
       bf1 = idS(5,i)
       Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
       Ay = atoms(4,int(basis(bf1,1)))
       Az = atoms(5,int(basis(bf1,1)))
       alpha = basis(bf1,12)
       nc1 = basis(bf1,22)
       do j=1,noS(2)
          bf2 = idS(2,j)
          Bx = atoms(3,int(basis(bf2,1)))
          By = atoms(4,int(basis(bf2,1)))
          Bz = atoms(5,int(basis(bf2,1)))
          do k2=1,int(basis(bf2,2))  
             beta = basis(bf2,10+k2)
             nc2 = basis(bf2,20+k2)
                zeta = alpha + beta
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
                noggsp = noggsp+1;
                ggsp(1:2,noggsp) = (/bf1,bf2/)
                ggsp(3,noggsp) = (alpha * Ax + beta*Bx)/(zeta);
                ggsp(4,noggsp) = (alpha * Ay + beta*By)/(zeta);
                ggsp(5,noggsp) = (alpha * Az + beta*Bz)/(zeta);
                ggsp(6:8,noggsp) = (/alpha,beta,zeta/);
                ggsp(9,noggsp) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)
                ggsp(11:13,noggsp) = atoms(3:5,int(basis(bf1,1)))-atoms(3:5,int(basis(bf2,1)));
                ggsp(14:16,noggsp) = (/sqrt(zeta),  1d0/sqrt(zeta),1d0/zeta/);
                ggsp(17:18,noggsp) = (/basis(bf1,1),basis(bf2,1)/);
                ggsp(19,noggsp) = 1d0/(zeta)**(5d0/2d0)
             end if
          end do
       end do
    end do
    


    ! First s shells, second Pp shell
    do i=1,noS(1)
       bf1 = idS(1,i)
       Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
       Ay = atoms(4,int(basis(bf1,1)))
       Az = atoms(5,int(basis(bf1,1)))
       do k=1,int(basis(bf1,2))
          alpha = basis(bf1,10+k)
          nc1 = basis(bf1,20+k)
          do j=1,noS(9)
             bf2 = idS(9,j)
             Bx = atoms(3,int(basis(bf2,1)))
             By = atoms(4,int(basis(bf2,1)))
             Bz = atoms(5,int(basis(bf2,1)))
             do k2=2,int(basis(bf2,2))  
               beta = basis(bf2,10+k2)
               nc2 = basis(bf2,20+k2)
                zeta = alpha + beta
               
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
                  noggsp = noggsp+1;
                ggsp(1:2,noggsp) = (/bf1,bf2/)
                ggsp(3,noggsp) = (alpha * Ax + beta*Bx)/(zeta);
                ggsp(4,noggsp) = (alpha * Ay + beta*By)/(zeta);
                ggsp(5,noggsp) = (alpha * Az + beta*Bz)/(zeta);
                ggsp(6:8,noggsp) = (/alpha,beta,zeta/);
                ggsp(9,noggsp) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)
                ggsp(11:13,noggsp) = atoms(3:5,int(basis(bf1,1)))-atoms(3:5,int(basis(bf2,1)));
                ggsp(14:16,noggsp) = (/sqrt(zeta),  1d0/sqrt(zeta),1d0/zeta/);
                ggsp(17:18,noggsp) = (/basis(bf1,1),basis(bf2,1)/);
                ggsp(19,noggsp) = 1d0/(zeta)**(5d0/2d0)
             end if
             end do
          end do
       end do
    end do

    ! First Ss shell, second Pp shell
    do i=1,noS(5)
       bf1 = idS(5,i)
       Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
       Ay = atoms(4,int(basis(bf1,1)))
       Az = atoms(5,int(basis(bf1,1)))
       alpha = basis(bf1,12)
       nc1 = basis(bf1,22)
       do j=1,noS(9)
          bf2 = idS(9,j)
          Bx = atoms(3,int(basis(bf2,1)))
          By = atoms(4,int(basis(bf2,1)))
          Bz = atoms(5,int(basis(bf2,1)))
          do k2=2,int(basis(bf2,2))  
             beta = basis(bf2,10+k2)
             nc2 = basis(bf2,20+k2)
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
                zeta = alpha + beta
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
                noggsp = noggsp+1;
                ggsp(1:2,noggsp) = (/bf1,bf2/)
                ggsp(3,noggsp) = (alpha * Ax + beta*Bx)/(zeta);
                ggsp(4,noggsp) = (alpha * Ay + beta*By)/(zeta);
                ggsp(5,noggsp) = (alpha * Az + beta*Bz)/(zeta);
                ggsp(6:8,noggsp) = (/alpha,beta,zeta/);
                ggsp(9,noggsp) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)
                ggsp(11:13,noggsp) = atoms(3:5,int(basis(bf1,1)))-atoms(3:5,int(basis(bf2,1)));
                ggsp(14:16,noggsp) = (/sqrt(zeta),  1d0/sqrt(zeta),1d0/zeta/);
                ggsp(17:18,noggsp) = (/basis(bf1,1),basis(bf2,1)/);
                ggsp(19,noggsp) = 1d0/(zeta)**(5d0/2d0)
            end if
          end do
       end do
    end do
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (sigspprint ==1) then
     
     print *, "Nosp = ", noggsp
     print *, "sp info"
     do i=1,noggsp
        print "(a2,i3,a2,i3,a4,3f23.15,a12,4f10.5)", "(", int(ggsp(1,i)), ",", int(ggsp(2,i)), ") = ",&
             ggsp(6,i), ggsp(7,i), ggsp(9,i), &
             "centered" , ggsp(3,i), ggsp(4,i), ggsp(5,i), sqrt(ggsp(4,i)**2+ggsp(5,i)**2+ggsp(3,i)**2)
     end do
     print *, ""
     

  end if
  
  return
  
end subroutine getsigspsp








!***************************************
subroutine getsigsppp(algchange,noggpp,  ggpp, nobasisfun,noatoms,atoms,basis,sigspprint,thresh)
  implicit none
  integer, intent(in)    :: nobasisfun, noatoms , thresh
  integer, intent(out)   ::  noggpp
  integer, intent(inout) :: algchange(1:30)
  real*8, intent(in)     ::  basis(1:nobasisfun,1:30), atoms(1:5,1:noatoms)
  real*8, intent(out)    :: ggpp(1:30,1:20*nobasisfun**2)
  integer                :: bf1, bf2, sigspprint,  i, j, k, k2, noS(0:11), idS(0:11,1:nobasisfun)
  real*8                 :: nc1, nc2,  n, Ax, Ay, Az, Bx, By, Bz, BAx, BAy, BAz, beta, alpha, RAB2, zeta, pi
  real*8 :: sqrtpi, sqrtzeta, rzeta, Px, Py, Pz
  
  real*8 :: pi3r2

pi3r2 = dacos(-1d0)**(3d0/2d0)

print *, "THRESH = ", thresh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This program constructs a list of significant SS, Ss and Sp shell pairs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      RR contains info about the SS shell pair and has elements (shellbf1, shellbf2, atoms centre,ramp deg1, rampdeg2, nc1*n2, Ax, Ay, Az, blank)
  !     
  !      RGSs contains info about the itive Ss shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
  !      RGSp contains info about the itive Sp shell pair and has elements (shellbf1, shellbf2, atomscenter, BAx, BAy, BAzramp degree, exp, nc1*nc2, algflag)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  pi = dacos(-1d0);
  
  noggpp = 0;   
 !ggpp   = 0d0; 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! algchange(1) is no of (Ss| shells concentric with large beta
!!! algchange(2) is tot no of (Sp| shells concentric 
!!! algchange(3) is no of (ss| which are concentric with large zeta
!!! algchange(4) is tot no of concentric (Ss| shell pairs (large and small beta)
!!! algchange(5) is no of (sp| which are concentric with large zeta
!!! algchange(6) is no of (pp| which are concentric with large zeta
!!! algchange(7) is no of concentric (Sp| shell pairs with large beta
!!! algchange(8) is
!!! algchange(9) is
!!! algchange(10) is
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  noS = 0;
  idS = 0;
  do i=1,nobasisfun
     noS(int(basis(i,3))) = noS(int(basis(i,3))) +1 
     idS(int(basis(i,3)), noS(int(basis(i,3)))) = i   
  end do
 
sqrtpi = sqrt(dacos(-1d0))

    
    
    !******************************************************
    !*** Constructing |pp) shell pairs  *******************
    !******************************************************
    
    noggpp = 0;
    
    ! First p shells, second p shell
    do i=1,noS(2)
       bf1 = idS(2,i)
       Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
       Ay = atoms(4,int(basis(bf1,1)))
       Az = atoms(5,int(basis(bf1,1)))
       do k=1,int(basis(bf1,2))
          alpha = basis(bf1,10+k)
          nc1 = basis(bf1,20+k)
          do j=i,noS(2)
             bf2 = idS(2,j)
             Bx = atoms(3,int(basis(bf2,1)))
             By = atoms(4,int(basis(bf2,1)))
             Bz = atoms(5,int(basis(bf2,1)))
             do k2=1,int(basis(bf2,2))  
                beta = basis(bf2,10+k2)
                nc2 = basis(bf2,20+k2)
                               zeta = alpha + beta
 
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
                   
                   noggpp = noggpp+1;
                   ggpp(1:2,noggpp) = (/bf1,bf2/);
                   ggpp(3,noggpp) = (alpha * Ax + beta*Bx)/(zeta)
                   ggpp(4,noggpp) = (alpha * Ay + beta*By)/(zeta)
                   ggpp(5,noggpp) = (alpha * Az + beta*Bz)/(zeta)
                   ggpp(6:8,noggpp) = (/alpha,beta,zeta/);
                   ggpp(9,noggpp) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)                 
                   ggpp(11:13,noggpp) = atoms(3:5,int(basis(bf1,1))) - atoms(3:5,int(basis(bf2,1)));
                       rzeta = 1d0/zeta; sqrtzeta = sqrt(zeta)
                   ggpp(14,noggpp) =  sqrtpi/(4d0*zeta*sqrt(zeta))*rzeta**2

 Px = (alpha * Ax + beta*Bx)/(zeta);
 Py = (alpha * Ay + beta*By)/(zeta);
 Pz = (alpha * Az + beta*Bz)/(zeta);

                   ggpp(15,noggpp) = Px-Ax;
                   ggpp(16,noggpp) = Py-Ay;
                   ggpp(17,noggpp) = Pz-Az;
                   ggpp(18,noggpp) = Px-Bx;
                   ggpp(19,noggpp) = Py-By;
                   ggpp(20,noggpp) = Pz-Bz;
                end if
             end do
          end do
       end do
  end do


  ! First Pp shells, second p shell
    do i=1,noS(9)
       bf1 = idS(9,i)
       Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
       Ay = atoms(4,int(basis(bf1,1)))
       Az = atoms(5,int(basis(bf1,1)))
       do k=2,int(basis(bf1,2))
          alpha = basis(bf1,10+k)
          nc1 = basis(bf1,20+k)
          do j=1,noS(2)
             bf2 = idS(2,j)
             Bx = atoms(3,int(basis(bf2,1)))
             By = atoms(4,int(basis(bf2,1)))
             Bz = atoms(5,int(basis(bf2,1)))
             do k2=1,int(basis(bf2,2))  
                beta = basis(bf2,10+k2)
                nc2 = basis(bf2,20+k2)
                               zeta = alpha + beta
 
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
                   
                   noggpp = noggpp+1;
                   ggpp(1:2,noggpp) = (/bf1,bf2/);
                   ggpp(3,noggpp) = (alpha * Ax + beta*Bx)/(zeta)
                   ggpp(4,noggpp) = (alpha * Ay + beta*By)/(zeta)
                   ggpp(5,noggpp) = (alpha * Az + beta*Bz)/(zeta)
                   ggpp(6:8,noggpp) = (/alpha,beta,zeta/);
                   ggpp(9,noggpp) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)                 
                   ggpp(11:13,noggpp) = atoms(3:5,int(basis(bf1,1))) - atoms(3:5,int(basis(bf2,1)));
                       rzeta = 1d0/zeta; sqrtzeta = sqrt(zeta)
                   ggpp(14,noggpp) =  sqrtpi/(4d0*zeta*sqrt(zeta))*rzeta**2

 Px = (alpha * Ax + beta*Bx)/(zeta);
 Py = (alpha * Ay + beta*By)/(zeta);
 Pz = (alpha * Az + beta*Bz)/(zeta);

                   ggpp(15,noggpp) = Px-Ax;
                   ggpp(16,noggpp) = Py-Ay;
                   ggpp(17,noggpp) = Pz-Az;
                   ggpp(18,noggpp) = Px-Bx;
                   ggpp(19,noggpp) = Py-By;
                   ggpp(20,noggpp) = Pz-Bz;
                end if
             end do
          end do
       end do
  end do
  ! First Pp shells, second Pp shell
    do i=1,noS(9)
       bf1 = idS(9,i)
       Ax = atoms(3,int(basis(bf1,1)))   ! Coordinates of the first Gaussian 
       Ay = atoms(4,int(basis(bf1,1)))
       Az = atoms(5,int(basis(bf1,1)))
       do k=2,int(basis(bf1,2))
          alpha = basis(bf1,10+k)
          nc1 = basis(bf1,20+k)
          do j=i,noS(9)
             bf2 = idS(9,j)
             Bx = atoms(3,int(basis(bf2,1)))
             By = atoms(4,int(basis(bf2,1)))
             Bz = atoms(5,int(basis(bf2,1)))
             do k2=2,int(basis(bf2,2))  
                beta = basis(bf2,10+k2)
                nc2 = basis(bf2,20+k2)
                
            RAB2 = sum((atoms(3:5,int(basis(bf2,1)))-atoms(3:5,int(basis(bf1,1))))**2);
             zeta = alpha + beta
        if (abs(nc1*nc2*pi3r2*exp(-alpha*beta/zeta*RAB2)/(zeta**1.5d0)) .gt. 10d0**(-thresh)) then
                   
                   noggpp = noggpp+1;
                   ggpp(1:2,noggpp) = (/bf1,bf2/);
                   ggpp(3,noggpp) = (alpha * Ax + beta*Bx)/(zeta)
                   ggpp(4,noggpp) = (alpha * Ay + beta*By)/(zeta)
                   ggpp(5,noggpp) = (alpha * Az + beta*Bz)/(zeta)
                   ggpp(6:8,noggpp) = (/alpha,beta,zeta/);
                   ggpp(9,noggpp) = nc1*nc2* exp(-(alpha*beta/(zeta))*RAB2)                 
                   ggpp(11:13,noggpp) = atoms(3:5,int(basis(bf1,1))) - atoms(3:5,int(basis(bf2,1)));
                       rzeta = 1d0/zeta; sqrtzeta = sqrt(zeta)
                   ggpp(14,noggpp) =  sqrtpi/(4d0*zeta*sqrt(zeta))*rzeta**2


 Px = (alpha * Ax + beta*Bx)/(zeta);
 Py = (alpha * Ay + beta*By)/(zeta);
 Pz = (alpha * Az + beta*Bz)/(zeta);

                   ggpp(15,noggpp) = Px-Ax;
                   ggpp(16,noggpp) = Py-Ay;
                   ggpp(17,noggpp) = Pz-Az;
                   ggpp(18,noggpp) = Px-Bx;
                   ggpp(19,noggpp) = Py-By;
                   ggpp(20,noggpp) = Pz-Bz;
               end if
             end do
          end do
       end do
  end do
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (sigspprint ==1) then
     

     print *, ""
     
     print *, "Nopp = ", noggpp
     print *, "pp info"
     do i=1,noggpp
        print "(a2,i3,a2,i3,a4,3f23.15,a12,4f10.5)", "(", int(ggpp(1,i)), ",", int(ggpp(2,i)), ") = ",&
             ggpp(6,i), ggpp(7,i), ggpp(9,i), &
             "centered" , ggpp(3,i), ggpp(4,i), ggpp(5,i), sqrt(ggpp(4,i)**2+ggpp(5,i)**2+ggpp(3,i)**2)
     end do
     
  end if
  
  return
  
end subroutine getsigsppp

