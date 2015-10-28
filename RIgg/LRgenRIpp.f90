
!***********************************************
!***********************************************
!***********************************************

subroutine calcLRgenRIppdigestoneRperatom(noggpp, ggpp,  &
     noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, lenmodelcontnewsort, contnewsortRg,  &
     Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, intcnt,Gtot, temp, & 
           newtwo,lengthtwoelec)
  implicit none
  real*8, intent(in)      :: ggpp(1:20,1:noggpp), Tcutoff,  atoms(1:5,1:noatoms)
  integer, intent(in)     :: noggpp, a, nobasisfun, noatoms, maxrampdegree
  integer                 :: i,j, ntot, k, l, m, n, warnswit(1:10),  bf1, bf2, bf3, bf4, bb3, bb4,intcnt
  real*8                  ::  nc1, Ax, Ay, Az, Q, rQA, x, y, rzeta
  real*8                  :: Qx, Qy, Qz, Cx, Cy, Cz, Dx, Dy, Dz, CDx, CDy, CDz, QAx, QAy, QAz, QA, QA2
  real*8                  :: gamma, delta, zeta, pi, rQ(0:4), pf(0:4)
  real*8                  :: ang(1:25), dangdx(1:25),dangdy(1:25), dangdz(1:25)
  real*8                  :: d2angdx2(1:25), d2angdxy(1:25), d2angdxz(1:25)
  real*8                  :: d2angdy2(1:25), d2angdyz(1:25), d2angdz2(1:25)
  real*8                  :: T, gf, contnewsortRg(1: 25*maxrgpa )
  integer                 :: lenmodelcontnewsort(1:3,1:maxrgpa), maxrgpa
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  real*8                  :: angcontpart(1:9,1:25), acum(1:20), stor(1:20)
  real*8                  :: dangdth(1:25), dangdph(1:25), temp(1:9,1:maxrgpa)
  real*8                  :: d2angdth2(1:25), d2angdph2(1:25),d2angdthph(1:25)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, dimension(:), allocatable :: stor2, acum2
  integer :: lengthtwoelec,braketno,cnt,cnt2
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calcLRgenRIpp calculates integrals of the form (RI|pp) in which the T parameter
!     T = |Q-A|^2 * zeta > 20. 
! 
! This is the long-range component of the interaction and is true when the two shell pairs don't
! overlap completely. 
!
! All integrals are correct. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  warnswit = 0; 
allocate( stor2(1:maxrgpa), acum2(1:maxrgpa))
  
!*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Calculates constants (can be precomputed if required)
!*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
  pi = dacos(-1d0); 
 
  do l=0,4;     pf(l) = 16*pi**2/(2d0*l+1);    end do;
   
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Main loop over shell quartets
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
  
   bf1 = lenmodelcontnewsort(1,1)

  do i=1,maxrgpa
       bf2 = lenmodelcontnewsort(2,i)
       stor2(i) = 2d0*Ptot(bf2,bf1)
  end do
  acum2 = 0d0;

        Ax      = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)
  !********** Goes through every pp shell pair
  do j=1,noggpp
     zeta     = ggpp(8,j);     
        QAx     = ggpp(3,j) - Ax;          QAy     = ggpp(4,j) - Ay;          QAz     = ggpp(5,j) - Az;
        QA2 = (QAx)**2+(QAy)**2+(QAz)**2;           
        !******** Decides whether to treat (RI|pp) on this atom with this shell pair by long-range method or 
        !******** short-range method (a different code entirely). 
        T      = zeta*QA2;
        
        if (T .ge. Tcutoff) then;   
     gamma    = ggpp(6,j);         delta = ggpp(7,j)
     CDx       = ggpp(11,j);        CDy = ggpp(12,j);   CDz = ggpp(13,j)
     
     bf3      = ggpp(1,j);        bf4      = ggpp(2,j)
 QA = sqrt(QA2); 
         gf = ggpp(14,j)*ggpp(9,j)

#ifdef INTCNTPRT
intcnt = intcnt + maxrgpa
#endif
           call calcder2thph(ang,dangdx,dangdy,dangdz,d2angdx2,d2angdxy,d2angdxz,d2angdy2,d2angdyz,d2angdz2,QAx,QAy,QAz,QA2,&
                  dangdth,dangdph,d2angdth2,d2angdph2,d2angdthph)
           
acum = 0d0;

stor(1:3) = 2d0*Ptot(bf3  ,bf4:bf4+2)
stor(4:6) = 2d0*Ptot(bf3+1,bf4:bf4+2)
stor(7:9) = 2d0*Ptot(bf3+2,bf4:bf4+2)
stor(10:12) = Palpha(bf1,bf4:bf4+2);
stor(13:15) = Pbeta(bf1,bf4:bf4+2);

           !*** Calculates inverse powers of QA
           rQ(0) = 1d0/QA;          do l=1,4;   rQ(l) = rQ(l-1)*rQ(0);    end do
           
k=0; do l=0,4; do m=-l,l; k =k +1;
angcontpart(1,k) = pf(l)*gf*rQ(l)*(ang(k)*(   0.5d0*(-2*CDx**2*gamma*delta+zeta) &
                      -0.25d0*(1+l)*(1*rQ(1)-(3+l)*QAx**2*rQ(3))         &
                      +0.5d0*(l+1)*CDx*QAx*rQ(1)*(delta-gamma)              )+ &
                      dangdx(k)*(-0.25*( 2*delta*CDx+(l+1)*QAx*rQ(1))) + &
                      dangdx(k)*(-0.25*(-2*gamma*CDx+(l+1)*QAx*rQ(1))) + &
                      d2angdx2(k)*0.25d0)

angcontpart(2,k) =   pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDx*CDy         &
                      +0.25d0*(l+1)*(l+3)*QAx*QAy*rQ(3)             &
                      +0.5d0*(l+1)*(CDx*QAy*delta-CDy*QAx*gamma)*rQ(1)  ) +  &
                      dangdy(k)*(-0.25*( 2*delta*CDx+(l+1)*QAx*rQ(1))) + &
                      dangdx(k)*(-0.25*(-2*gamma*CDy+(l+1)*QAy*rQ(1))) + &
                      d2angdxy(k)*0.25d0     )    ;


angcontpart(3,k) =   pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDx*CDz &
                      +0.25d0*(l+1)*(l+3)*QAx*QAz*rQ(3)                            &
                      +0.5d0*(l+1)*(CDx*QAz*delta-CDz*QAx*gamma)*rQ(1)  ) + &
                      dangdz(k)*(-0.25*( 2*delta*CDx+(l+1)*QAx*rQ(1))) + &
                      dangdx(k)*(-0.25*(-2*gamma*CDz+(l+1)*QAz*rQ(1))) + &
                      d2angdxz(k)*0.25d0) 

angcontpart(4,k) =   pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDy*CDx         &
                      +0.25d0*(l+1)*(l+3)*QAy*QAx*rQ(3)             &
                      +0.5d0*(l+1)*(CDy*QAx*delta-CDx*QAy*gamma)*rQ(1)  ) +  &
                      dangdx(k)*(-0.25*( 2*delta*CDy+(l+1)*QAy*rQ(1))) + &
                      dangdy(k)*(-0.25*(-2*gamma*CDx+(l+1)*QAx*rQ(1))) + &
                      d2angdxy(k)*0.25d0 )

angcontpart(5,k) =    pf(l)*gf*rQ(l)*(ang(k)*(   0.5d0*(-2*CDy**2*gamma*delta+zeta) &
                      -0.25d0*(1+l)*(1*rQ(1)-(3+l)*QAy**2*rQ(3))         &
                      +0.5d0*(l+1)*CDy*QAy*rQ(1)*(delta-gamma)              )+&
                      dangdy(k)*(-0.25*( 2*delta*CDy+(l+1)*QAy*rQ(1))) + &
                      dangdy(k)*(-0.25*(-2*gamma*CDy+(l+1)*QAy*rQ(1))) + &
                      d2angdy2(k)*0.25d0  )

angcontpart(6,k) = pf(l)*gf*rQ(l)*(ang(k)*( -gamma*delta*CDy*CDz          &
                      +0.25d0*(l+1)*(l+3)*QAy*QAz*rQ(3)                            &
                      +0.5d0*(l+1)*(CDy*QAz*delta-CDz*QAy*gamma)*rQ(1)  )  + &
                      dangdz(k)*(-0.25*( 2*delta*CDy+(l+1)*QAy*rQ(1))) + &
                      dangdy(k)*(-0.25*(-2*gamma*CDz+(l+1)*QAz*rQ(1))) + &
                      d2angdyz(k)*0.25d0  )

angcontpart(7,k) = pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDz*CDx         &
                      +0.25d0*(l+1)*(l+3)*QAx*QAz*rQ(3)                           &
                      +0.5d0*(l+1)*(CDz*QAx*delta-CDx*QAz*gamma)*rQ(1)  )  +&
                      dangdx(k)*(-0.25*( 2*delta*CDz+(l+1)*QAz*rQ(1))) + &
                      dangdz(k)*(-0.25*(-2*gamma*CDx+(l+1)*QAx*rQ(1))) + &
                      d2angdxz(k)*0.25d0  )   ;

angcontpart(8,k) =  pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDy*CDz     &
                      +0.25d0*(l+1)*(l+3)*QAy*QAz*rQ(3)                           &
                      +0.5d0*(l+1)*(CDz*QAy*delta-CDy*QAz*gamma)*rQ(1)  )  +      &
                      dangdy(k)*(-0.25*( 2*delta*CDz+(l+1)*QAz*rQ(1))) + &
                      dangdz(k)*(-0.25*(-2*gamma*CDy+(l+1)*QAy*rQ(1))) + &
                      d2angdyz(k)*0.25d0 )     

angcontpart(9,k) =  pf(l)*gf*rQ(l)*(ang(k)*(   0.5d0*(-2*CDz**2*gamma*delta+zeta)  &
                      -0.25d0*(1+l)*(1*rQ(1)-(3+l)*QAz**2*rQ(3))         &
                      +0.5d0*(l+1)*CDz*QAz*rQ(1)*(delta-gamma)              )+ &
                      dangdz(k)*(-0.25*( 2*delta*CDz+(l+1)*QAz*rQ(1))) + &
                      dangdz(k)*(-0.25*(-2*gamma*CDz+(l+1)*QAz*rQ(1))) + &
                      d2angdz2(k)*0.25d0  )   
     
end do; end do

call dgemm('N','N',9,maxrgpa,25,1d0, angcontpart, &
       9, contnewsortRg(lenmodelcontnewsort(3,1)+1: lenmodelcontnewsort(3,maxrgpa)+25),&
       25,0d0, temp,9)

         do i=1,maxrgpa
            bf2 = lenmodelcontnewsort(2,i)



if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if;
do cnt=1,3;
 do cnt2 = 1,3;
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3+cnt-1,bf4+cnt2-1,temp(3*cnt+cnt2-3,i),braketno)
 end do
end do
end if

     acum2(i) = acum2(i) + sum(stor(1:9)*temp(1:9,i))

     acum(1:3) = acum(1:3) + stor2(i)*temp(1:3,i)
     acum(4:6) = acum(4:6) + stor2(i)*temp(4:6,i)
     acum(7:9) = acum(7:9) + stor2(i)*temp(7:9,i)

     acum(10)   = acum(10) -sum(Palpha(bf2,bf4:bf4+2)*temp(1:3,i))
     acum(11)   = acum(11) -sum(Pbeta(bf2,bf4:bf4+2) *temp(1:3,i))
     acum(12)   = acum(12) -sum(Palpha(bf2,bf4:bf4+2)*temp(4:6,i))
     acum(13)   = acum(13)  -sum(Pbeta(bf2,bf4:bf4+2) *temp(4:6,i))
     acum(14)   = acum(14) -sum(Palpha(bf2,bf4:bf4+2)*temp(7:9,i))
     acum(15)   = acum(15) -sum(Pbeta(bf2,bf4:bf4+2) *temp(7:9,i))

     Galpha(bf4:bf4+2 ,bf1) = Galpha(bf4:bf4+2 ,bf1) &
         -Palpha(bf2,bf3)*temp(1:3,i) -Palpha(bf2,bf3+1)*temp(4:6,i) -Palpha(bf2,bf3+2)*temp(7:9,i) 
     Gbeta(bf4:bf4+2 ,bf1)  = Gbeta(bf4:bf4+2 ,bf1)  &
         - Pbeta(bf2,bf3)*temp(1:3,i) - Pbeta(bf2,bf3+1)*temp(4:6,i) - Pbeta(bf2,bf3+2)*temp(7:9,i)

     Galpha(bf2,bf4:bf4+2) = Galpha(bf2,bf4:bf4+2) &
         -Palpha(bf3 ,bf1)*temp(1:3,i) -Palpha(bf3+1 ,bf1)*temp(4:6,i) -Palpha(bf3+2 ,bf1)*temp(7:9,i) 
     Gbeta(bf2,bf4:bf4+2)  = Gbeta(bf2,bf4:bf4+2)  &
         - Pbeta(bf3 ,bf1)*temp(1:3,i) - Pbeta(bf3+1 ,bf1)*temp(4:6,i) - Pbeta(bf3+2 ,bf1)*temp(7:9,i)


     Galpha(bf2,bf3)   = Galpha(bf2,bf3)   -sum(stor(10:12)*temp(1:3,i))
     Gbeta(bf2,bf3)    = Gbeta(bf2,bf3)    -sum(stor(13:15)*temp(1:3,i))
     Galpha(bf2,bf3+1) = Galpha(bf2,bf3+1) -sum(stor(10:12)*temp(4:6,i))
     Gbeta(bf2,bf3+1)  = Gbeta(bf2,bf3+1)  -sum(stor(13:15)*temp(4:6,i))
     Galpha(bf2,bf3+2) = Galpha(bf2,bf3+2) -sum(stor(10:12)*temp(7:9,i))
     Gbeta(bf2,bf3+2)  = Gbeta(bf2,bf3+2)  -sum(stor(13:15)*temp(7:9,i))

          end do


     Gtot(bf3  ,bf4:bf4+2) = Gtot(bf3  ,bf4:bf4+2) + acum(1:3)
     Gtot(bf3+1,bf4:bf4+2) = Gtot(bf3+1,bf4:bf4+2) + acum(4:6)
     Gtot(bf3+2,bf4:bf4+2) = Gtot(bf3+2,bf4:bf4+2) + acum(7:9)
     Galpha(bf3  ,bf1)  = Galpha(bf3   ,bf1) + acum(10)
     Gbeta(bf3   ,bf1)   = Gbeta(bf3   ,bf1)  + acum(11)
     Galpha(bf3+1 ,bf1)  = Galpha(bf3+1 ,bf1) + acum(12)
     Gbeta(bf3+1 ,bf1)   = Gbeta(bf3+1 ,bf1)  + acum(13)
     Galpha(bf3+2 ,bf1)  = Galpha(bf3+2 ,bf1) + acum(14)
     Gbeta(bf3+2 ,bf1)   = Gbeta(bf3+2 ,bf1)  + acum(15)


        end if            !*** Ends decision of LR vs SR code
  end do           !*** Ends loop over pp shell pairs
  
 do i=1,maxrgpa
     bf2 = lenmodelcontnewsort(2,i)
     Gtot(bf2,bf1) = Gtot(bf2,bf1) + acum2(i)
  end do

deallocate(stor2, acum2)
  
  return
end subroutine calcLRgenRIppdigestoneRperatom

!***********************************************


!***********************************************
!***********************************************
!***********************************************

subroutine calcLRgenRIppdigest(noggpp, ggpp,  &
     noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, lenmodelcontnewsort, contnewsortRg, &
     Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, &
     ang,dangdx, dangdy, dangdz, d2angdx2, d2angdxy,d2angdxz,d2angdy2,d2angdyz,d2angdz2, & 
     dangdth,dangdph,d2angdth2,d2angdph2,d2angdthph,intcnt,Gtot)
  implicit none
  real*8, intent(in)      :: ggpp(1:20,1:noggpp), Tcutoff,  atoms(1:5,1:noatoms)
  integer, intent(in)     :: noggpp, a, nobasisfun, noatoms, maxrampdegree
  integer                 :: i,j, ntot, k, l, m, n, warnswit(1:10),  bf1, bf2, bf3, bf4, bb3, bb4,intcnt
  real*8                  ::  nc1, nc2, Ax, Ay, Az, Q, rQA, x, y
  real*8                  :: Qx, Qy, Qz, Cx, Cy, Cz, Dx, Dy, Dz, CDx, CDy, CDz, QAx, QAy, QAz, QA, QA2
  real*8                  :: gamma, delta, zeta, sqrtzeta, pi, sqrtpi, rzeta
  real*8                  :: rQ(0:4), pf(0:4)
  real*8                  :: ang(1:25), dangdx(1:25),dangdy(1:25), dangdz(1:25)
  real*8                  :: d2angdx2(1:25), d2angdxy(1:25), d2angdxz(1:25)
  real*8                  :: d2angdy2(1:25), d2angdyz(1:25), d2angdz2(1:25)
  real*8                  :: T, gf, contnewsortRg(1: 25*maxrgpa )
  integer                 :: lenmodelcontnewsort(1:3,1:maxrgpa), maxrgpa, v, st
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  real*8                  :: angcontpart(1:9,1:25)
  real*8                  :: dangdth(1:25), dangdph(1:25), temp(1:9,1:maxrgpa)
  real*8                  :: d2angdth2(1:25), d2angdph2(1:25),d2angdthph(1:25)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calcLRgenRIpp calculates integrals of the form (RI|pp) in which the T parameter
!     T = |Q-A|^2 * zeta > 20. 
! 
! This is the long-range component of the interaction and is true when the two shell pairs don't
! overlap completely. 
!
! All integrals are correct. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  warnswit = 0; 
  
!*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Calculates constants (can be precomputed if required)
!*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
  pi = dacos(-1d0); sqrtpi = sqrt(pi); 
 
  do l=0,4;     pf(l) = 16*pi**2/(2d0*l+1);    end do;
   
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Main loop over shell quartets
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
  
  
  !********** Goes through every pp shell pair
  do j=1,noggpp
     bf3      = ggpp(1,j);        bf4      = ggpp(2,j)
     nc2      = ggpp(9,j)
     zeta     = ggpp(8,j);      sqrtzeta = sqrt(zeta);     rzeta = 1d0/zeta;
     gamma    = ggpp(6,j);         delta = ggpp(7,j)
     Qx       = ggpp(3,j);         Qy = ggpp(4,j);    Qz = ggpp(5,j)
     CDx       = ggpp(11,j);        CDy = ggpp(12,j);   CDz = ggpp(13,j)
     
        Ax      = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)
        QAx     = Qx - Ax;          QAy     = Qy - Ay;          QAz     = Qz - Az;
        QA = sqrt((QAx)**2+(QAy)**2+(QAz)**2);   QA2 = QA**2;    rQA  = 1d0/QA
        
        !******** Decides whether to treat (RI|pp) on this atom with this shell pair by long-range method or 
        !******** short-range method (a different code entirely). 
        T      = zeta*QA**2;
        
        if (T .ge. Tcutoff) then;
           

intcnt = intcnt + maxrgpa
           call calcder2thph(ang,dangdx,dangdy,dangdz,d2angdx2,d2angdxy,d2angdxz,d2angdy2,d2angdyz,d2angdz2,QAx,QAy,QAz,QA2,&
                  dangdth,dangdph,d2angdth2,d2angdph2,d2angdthph)
           
           !*** Calculates inverse powers of QA
           rQ(0) = 1d0*rQA;           do l=1,4;   rQ(l) = rQ(l-1)*rQA;    end do
           
           !*** Gaussian prefactor
           gf = sqrtpi/(4d0*zeta*sqrtzeta)*nc2*rzeta**2
           
k=0; do l=0,4; do m=-l,l; k =k +1;
angcontpart(1,k) = pf(l)*gf*rQ(l)*(ang(k)*(   0.5d0*(-2*CDx**2*gamma*delta+zeta) &
                      -0.25d0*(1+l)*(1*rQ(1)-(3+l)*QAx**2*rQ(3))         &
                      +0.5d0*(l+1)*CDx*QAx*rQ(1)*(delta-gamma)              )+ &
                      dangdx(k)*(-0.25*( 2*delta*CDx+(l+1)*QAx*rQ(1))) + &
                      dangdx(k)*(-0.25*(-2*gamma*CDx+(l+1)*QAx*rQ(1))) + &
                      d2angdx2(k)*0.25d0)

angcontpart(2,k) =   pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDx*CDy         &
                      +0.25d0*(l+1)*(l+3)*QAx*QAy*rQ(3)             &
                      +0.5d0*(l+1)*(CDx*QAy*delta-CDy*QAx*gamma)*rQ(1)  ) +  &
                      dangdy(k)*(-0.25*( 2*delta*CDx+(l+1)*QAx*rQ(1))) + &
                      dangdx(k)*(-0.25*(-2*gamma*CDy+(l+1)*QAy*rQ(1))) + &
                      d2angdxy(k)*0.25d0     )    ;


angcontpart(3,k) =   pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDx*CDz &
                      +0.25d0*(l+1)*(l+3)*QAx*QAz*rQ(3)                            &
                      +0.5d0*(l+1)*(CDx*QAz*delta-CDz*QAx*gamma)*rQ(1)  ) + &
                      dangdz(k)*(-0.25*( 2*delta*CDx+(l+1)*QAx*rQ(1))) + &
                      dangdx(k)*(-0.25*(-2*gamma*CDz+(l+1)*QAz*rQ(1))) + &
                      d2angdxz(k)*0.25d0) 

angcontpart(4,k) =   pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDy*CDx         &
                      +0.25d0*(l+1)*(l+3)*QAy*QAx*rQ(3)             &
                      +0.5d0*(l+1)*(CDy*QAx*delta-CDx*QAy*gamma)*rQ(1)  ) +  &
                      dangdx(k)*(-0.25*( 2*delta*CDy+(l+1)*QAy*rQ(1))) + &
                      dangdy(k)*(-0.25*(-2*gamma*CDx+(l+1)*QAx*rQ(1))) + &
                      d2angdxy(k)*0.25d0 )

angcontpart(5,k) =    pf(l)*gf*rQ(l)*(ang(k)*(   0.5d0*(-2*CDy**2*gamma*delta+zeta) &
                      -0.25d0*(1+l)*(1*rQ(1)-(3+l)*QAy**2*rQ(3))         &
                      +0.5d0*(l+1)*CDy*QAy*rQ(1)*(delta-gamma)              )+&
                      dangdy(k)*(-0.25*( 2*delta*CDy+(l+1)*QAy*rQ(1))) + &
                      dangdy(k)*(-0.25*(-2*gamma*CDy+(l+1)*QAy*rQ(1))) + &
                      d2angdy2(k)*0.25d0  )

angcontpart(6,k) = pf(l)*gf*rQ(l)*(ang(k)*( -gamma*delta*CDy*CDz          &
                      +0.25d0*(l+1)*(l+3)*QAy*QAz*rQ(3)                            &
                      +0.5d0*(l+1)*(CDy*QAz*delta-CDz*QAy*gamma)*rQ(1)  )  + &
                      dangdz(k)*(-0.25*( 2*delta*CDy+(l+1)*QAy*rQ(1))) + &
                      dangdy(k)*(-0.25*(-2*gamma*CDz+(l+1)*QAz*rQ(1))) + &
                      d2angdyz(k)*0.25d0  )

angcontpart(7,k) = pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDz*CDx         &
                      +0.25d0*(l+1)*(l+3)*QAx*QAz*rQ(3)                           &
                      +0.5d0*(l+1)*(CDz*QAx*delta-CDx*QAz*gamma)*rQ(1)  )  +&
                      dangdx(k)*(-0.25*( 2*delta*CDz+(l+1)*QAz*rQ(1))) + &
                      dangdz(k)*(-0.25*(-2*gamma*CDx+(l+1)*QAx*rQ(1))) + &
                      d2angdxz(k)*0.25d0  )   ;

angcontpart(8,k) =  pf(l)*gf*rQ(l)*(ang(k)*(  -gamma*delta*CDy*CDz     &
                      +0.25d0*(l+1)*(l+3)*QAy*QAz*rQ(3)                           &
                      +0.5d0*(l+1)*(CDz*QAy*delta-CDy*QAz*gamma)*rQ(1)  )  +      &
                      dangdy(k)*(-0.25*( 2*delta*CDz+(l+1)*QAz*rQ(1))) + &
                      dangdz(k)*(-0.25*(-2*gamma*CDy+(l+1)*QAy*rQ(1))) + &
                      d2angdyz(k)*0.25d0 )     

angcontpart(9,k) =  pf(l)*gf*rQ(l)*(ang(k)*(   0.5d0*(-2*CDz**2*gamma*delta+zeta)  &
                      -0.25d0*(1+l)*(1*rQ(1)-(3+l)*QAz**2*rQ(3))         &
                      +0.5d0*(l+1)*CDz*QAz*rQ(1)*(delta-gamma)              )+ &
                      dangdz(k)*(-0.25*( 2*delta*CDz+(l+1)*QAz*rQ(1))) + &
                      dangdz(k)*(-0.25*(-2*gamma*CDz+(l+1)*QAz*rQ(1))) + &
                      d2angdz2(k)*0.25d0  )   
     
end do; end do

temp = 0d0;
call dgemm('N','N',9,maxrgpa,25,1d0, angcontpart, &
       9, contnewsortRg(lenmodelcontnewsort(3,1)+1: lenmodelcontnewsort(3,maxrgpa)+25),&
       25,1d0, temp,9)

         do i=1,maxrgpa
            bf1 = lenmodelcontnewsort(1,i)
            bf2 = lenmodelcontnewsort(2,i)

     Gtot(bf1,bf2) = Gtot(bf1,bf2) + 2d0*sum(Ptot(bf3,bf4:bf4+2)*temp(1:3,i))  & 
          + 2d0*sum(Ptot(bf3+1,bf4:bf4+2)*temp(4:6,i)) + 2d0*sum(Ptot(bf3+2,bf4:bf4+2)*temp(7:9,i))

     Gtot(bf3  ,bf4:bf4+2) = Gtot(bf3  ,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp(1:3,i)
     Gtot(bf3+1,bf4:bf4+2) = Gtot(bf3+1,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp(4:6,i)
     Gtot(bf3+2,bf4:bf4+2) = Gtot(bf3+2,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp(7:9,i)

     Galpha(bf1,bf4:bf4+2) = Galpha(bf1,bf4:bf4+2) &
         -Palpha(bf2,bf3)*temp(1:3,i) -Palpha(bf2,bf3+1)*temp(4:6,i) -Palpha(bf2,bf3+2)*temp(7:9,i) 
     Gbeta(bf1,bf4:bf4+2)  = Gbeta(bf1,bf4:bf4+2)  &
         - Pbeta(bf2,bf3)*temp(1:3,i) - Pbeta(bf2,bf3+1)*temp(4:6,i) - Pbeta(bf2,bf3+2)*temp(7:9,i)

     Galpha(bf2,bf4:bf4+2) = Galpha(bf2,bf4:bf4+2) &
         -Palpha(bf1,bf3)*temp(1:3,i) -Palpha(bf1,bf3+1)*temp(4:6,i) -Palpha(bf1,bf3+2)*temp(7:9,i) 
     Gbeta(bf2,bf4:bf4+2)  = Gbeta(bf2,bf4:bf4+2)  &
         - Pbeta(bf1,bf3)*temp(1:3,i) - Pbeta(bf1,bf3+1)*temp(4:6,i) - Pbeta(bf1,bf3+2)*temp(7:9,i)

     Galpha(bf1,bf3)   = Galpha(bf1,bf3)   -sum(Palpha(bf2,bf4:bf4+2)*temp(1:3,i))
     Gbeta(bf1,bf3)    = Gbeta(bf1,bf3)    -sum(Pbeta(bf2,bf4:bf4+2) *temp(1:3,i))
     Galpha(bf1,bf3+1) = Galpha(bf1,bf3+1) -sum(Palpha(bf2,bf4:bf4+2)*temp(4:6,i))
     Gbeta(bf1,bf3+1)  = Gbeta(bf1,bf3+1)  -sum(Pbeta(bf2,bf4:bf4+2) *temp(4:6,i))
     Galpha(bf1,bf3+2) = Galpha(bf1,bf3+2) -sum(Palpha(bf2,bf4:bf4+2)*temp(7:9,i))
     Gbeta(bf1,bf3+2)  = Gbeta(bf1,bf3+2)  -sum(Pbeta(bf2,bf4:bf4+2) *temp(7:9,i))

     Galpha(bf2,bf3)   = Galpha(bf2,bf3)   -sum(Palpha(bf1,bf4:bf4+2)*temp(1:3,i))
     Gbeta(bf2,bf3)    = Gbeta(bf2,bf3)    -sum(Pbeta(bf1,bf4:bf4+2) *temp(1:3,i))
     Galpha(bf2,bf3+1) = Galpha(bf2,bf3+1) -sum(Palpha(bf1,bf4:bf4+2)*temp(4:6,i))
     Gbeta(bf2,bf3+1)  = Gbeta(bf2,bf3+1)  -sum(Pbeta(bf1,bf4:bf4+2) *temp(4:6,i))
     Galpha(bf2,bf3+2) = Galpha(bf2,bf3+2) -sum(Palpha(bf1,bf4:bf4+2)*temp(7:9,i))
     Gbeta(bf2,bf3+2)  = Gbeta(bf2,bf3+2)  -sum(Pbeta(bf1,bf4:bf4+2) *temp(7:9,i))

          end do
        end if            !*** Ends decision of LR vs SR code
  end do           !*** Ends loop over pp shell pairs
  
  
  return
end subroutine calcLRgenRIppdigest

!***********************************************

