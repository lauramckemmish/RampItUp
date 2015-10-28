
subroutine calcLRgenRIspdigestoneRperatom(noggsp, ggsp, noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, &
    lenmodelcontnewsort, contnewsortRg, Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,intcnt,Gtot,temp, & 
           newtwo,lengthtwoelec)
  implicit none
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun),ggsp(1:20,1:noggsp), atoms(1:5,1:noatoms), Tcutoff
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  integer, intent(in)     :: noatoms, maxrampdegree, noggsp,a, nobasisfun
  integer                 :: i,j, ntot, k, l, m, n, swit, warnswit(1:10), bf1, bf2, bf3, bf4, bb4
  integer                 :: lenmodelcontnewsort(1:3,1:maxrgpa), maxrgpa, v, st,intcnt
  real*8                  :: nc2, Ax, Ay, Az, Q, rQA, intt, x, y,  CDx, CDy, CDz, QAx, QAy, QAz, QA, QA2
  real*8                  :: gamma, delta, zeta, sqrtzeta, cf, rzeta, sqrtpir4, pi
  real*8                  :: ang(1:25), dangdx(1:25), dangdy(1:25), dangdz(1:25)
  real*8                  :: dangdth(1:25), dangdph(1:25), rQ(0:4), pf(0:4), T, totrf, gf
  real*8                  :: contnewsortRg(1: 25*maxrgpa ), angpart(1:3,1:25), a1, a2, t1, t2
  real*8                  :: acum(1:11), stor(1:11), t11(1:4), temp(1:3,1:maxrgpa)
  real*8, dimension(:), allocatable :: acum2, stor2

  integer :: lengthtwoelec, braketno, cnt2
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calcLRgenRIsp calculates integrals of the form (RI|sp) in which the T parameter
  !     T = |Q-A|^2 * zeta > 20. 
  ! 
  ! This is the long-range component of the interaction and is true when the two shell pairs don't
  ! overlap completely. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(acum2(1:maxrgpa), stor2(1:maxrgpa))
  warnswit = 0; 
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Calculates constants (can be precomputed if required)
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
  pi = dacos(-1d0);
  sqrtpir4 = sqrt(pi)/4d0;
  
  do l=0,4;   pf(l) = 16*pi**2/(2d0*l+1);    end do;

  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Main loop over shell quartets
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
     
bf1 = lenmodelcontnewsort(1,1)
acum2 = 0d0;
do i=1,maxrgpa
    bf2 = lenmodelcontnewsort(2,i)
    stor2(i) = 2d0*Ptot(bf2,bf1)
end do
Ax      = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)     
  !********** Goes through every sp shell pair
  do j=1,noggsp
     zeta     = ggsp(8,j);     
        QAx = ggsp(3,j) - Ax;                      QAy = ggsp(4,j) - Ay;      QAz = ggsp(5,j) - Az;
        QA2 = (QAx)**2+(QAy)**2+(QAz)**2;  
        
        !******** Decides whether to treat (RI|sp) on this atom with this shell pair by long-range method or 
        !******** short-range method (a different code entirely). 
        T      = zeta*QA2;
        
 if (T .ge. Tcutoff) then;
  bf3   = ggsp(1,j);      bf4   = ggsp(2,j)
  gamma = ggsp(6,j);      
  CDx   = ggsp(11,j);      CDy   = ggsp(12,j);    CDz = ggsp(13,j);   

acum = 0d0
stor(1)    = Palpha(bf3,bf1)
stor(2)    = Pbeta(bf3,bf1)
stor(3:5)  = Palpha(bf4:bf4+2,bf1)
stor(6:8)  = Pbeta(bf4:bf4+2,bf1)
stor(9:11) = Ptot(bf3,bf4:bf4+2)

     !*** Gaussian prefactor
     gf = sqrtpir4/(ggsp(14,j))*ggsp(16,j)**2*ggsp(9,j)
     
#ifdef INTCNTPRT
intcnt = intcnt + maxrgpa
#endif
           QA  = sqrt(QA2);    rQA  = 1d0/QA
           
           call calcderthph(dangdx,dangdy,dangdz,ang,QAx,QAy,QAz,QA,QA2,dangdth,dangdph)
           
           !*** Calculates inverse powers of QA
           rQ(0) = rQA;           do l=1,4;    rQ(l) = rQ(l-1)*rQA;              end do
         
k=0; 
do l=0,4;
t11(1) = pf(l)*gf*rQ(l); 
t11(2) = (gamma*CDx - QAx*rQ(1)*(l+1)*0.5d0);
t11(3) = (gamma*CDy - QAy*rQ(1)*(l+1)*0.5d0);
t11(4) = (gamma*CDz - QAz*rQ(1)*(l+1)*0.5d0);
do m=-l,l; k=k+1;
  angpart(1,k) = t11(1)*(0.5d0*dangdx(k) + ang(k)*t11(2));
  angpart(2,k) = t11(1)*(0.5d0*dangdy(k) + ang(k)*t11(3));
  angpart(3,k) = t11(1)*(0.5d0*dangdz(k) + ang(k)*t11(4));
end do; 
end do;
        
call dgemm('N','N',3,maxrgpa,25,1d0, angpart, &
       3, contnewsortRg(lenmodelcontnewsort(3,1)+1: lenmodelcontnewsort(3,maxrgpa)+25),&
       25,0d0, temp,3)
             
 do i=1,maxrgpa
     bf2 = lenmodelcontnewsort(2,i)

     acum(1:3) = acum(1:3) + stor2(i)*temp(1:3,i)
     acum(4:6) = acum(4:6) -Palpha(bf2,bf3)*temp(1:3,i)
     acum(7:9) = acum(7:9) - Pbeta(bf2,bf3)*temp(1:3,i)
     acum(10)  = acum(10) - sum(Palpha(bf4:bf4+2,bf2)*temp(1:3,i))
     acum(11)  = acum(11)  - sum(Pbeta(bf4:bf4+2,bf2)*temp(1:3,i))

     Galpha(bf2,bf4:bf4+2) = Galpha(bf2,bf4:bf4+2) - stor(1)*temp(1:3,i)
     Gbeta(bf2,bf4:bf4+2)  = Gbeta(bf2,bf4:bf4+2)  - stor(2)*temp(1:3,i)

     Galpha(bf2,bf3) = Galpha(bf2,bf3) - sum(stor(3:5)*temp(1:3,i))
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - sum(stor(6:8)*temp(1:3,i))
     acum2(i) =acum2(i) + 2d0*sum(stor(9:11)*temp(1:3,i))


if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;        end if;
do cnt2=1,3; 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+cnt2-1,temp(cnt2,i),braketno)
end do; 
end if


 end do

 Gtot(bf4:bf4+2,bf3)   = Gtot(bf4:bf4+2,bf3) + acum(1:3)
 Galpha(bf4:bf4+2,bf1) = Galpha(bf4:bf4+2,bf1) +acum(4:6)
 Gbeta(bf4:bf4+2,bf1)  = Gbeta(bf4:bf4+2,bf1)  +acum(7:9)
 Galpha(bf1,bf3)       = Galpha(bf1,bf3) +acum(10)
 Gbeta(bf1,bf3)        = Gbeta(bf1,bf3)  +acum(11)

 end if            !*** Ends selection of LR vs SR code
end do           !*** Ends loop over sp shell pairs
    
do i=1,maxrgpa
 bf2           = lenmodelcontnewsort(2,i)
 Gtot(bf2,bf1) = Gtot(bf2,bf1) + acum2(i)
end do

deallocate(acum2,stor2)
  return
end subroutine calcLRgenRIspdigestoneRperatom

!***************************************************************





subroutine calcLRgenRIspdigest(noggsp, ggsp, noatoms, atoms, maxrampdegree, Tcutoff,a,maxrgpa, &
    lenmodelcontnewsort, contnewsortRg, Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, &
    ang,dangdx, dangdy, dangdz, dangdth,dangdph,intcnt,Gtot)
  implicit none
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun),ggsp(1:20,1:noggsp), atoms(1:5,1:noatoms), Tcutoff
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  integer, intent(in)     :: noatoms, maxrampdegree, noggsp,a, nobasisfun
  integer                 :: i,j, ntot, k, l, m, n, swit, warnswit(1:10), bf1, bf2, bf3, bf4, bb4
  integer                 :: lenmodelcontnewsort(1:3,1:maxrgpa), maxrgpa, v, st,intcnt
  real*8                  :: nc2, Ax, Ay, Az, Q, rQA, intt, x, y, temp(1:3,1:maxrgpa), CDx, CDy, CDz, QAx, QAy, QAz, QA, QA2
  real*8                  :: Qx, Qy, Qz, Cx, Cy, Cz, Dx, Dy, Dz, gamma, delta, zeta, sqrtzeta, cf, rzeta, sqrtpir4
  real*8                  :: ang(1:25), dangdx(1:25), dangdy(1:25), dangdz(1:25)
  real*8                  :: dangdth(1:25), dangdph(1:25), rQ(0:4), pf(0:4), T, totrf, gf
  real*8                  :: contnewsortRg(1: 25*maxrgpa ), angpart(1:3,1:25), a1, a2, t1, t2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calcLRgenRIsp calculates integrals of the form (RI|sp) in which the T parameter
  !     T = |Q-A|^2 * zeta > 20. 
  ! 
  ! This is the long-range component of the interaction and is true when the two shell pairs don't
  ! overlap completely. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  warnswit = 0; 
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Calculates constants (can be precomputed if required)
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
  sqrtpir4 = sqrt(dacos(-1d0))/4d0;
  
  do l=0,4;   pf(l) = 16*(dacos(-1d0)**2)/(2d0*l+1);    end do;

  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
!!!**** Main loop over shell quartets
  !*******!!!!!!!!!!!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!!
     
Ax      = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)     
  !********** Goes through every sp shell pair
  do j=1,noggsp
     bf3      = ggsp(1,j);        bf4      = ggsp(2,j)
     nc2      = ggsp(9,j)
     zeta     = ggsp(8,j);      sqrtzeta = sqrt(zeta);     rzeta = 1d0/zeta;
     gamma = ggsp(6,j);         delta = ggsp(7,j)
     Qx    = ggsp(3,j);         Qy = ggsp(4,j);    Qz = ggsp(5,j)
     CDx       = ggsp(11,j);        CDy = ggsp(12,j);   CDz = ggsp(13,j)
     !*** Gaussian prefactor
     gf = sqrtpir4/(sqrtzeta)*rzeta**2*nc2
     
        QAx = Qx - Ax;                      QAy = Qy - Ay;      QAz = Qz - Az;
        QA2 = (QAx)**2+(QAy)**2+(QAz)**2;  
        
        !******** Decides whether to treat (RI|sp) on this atom with this shell pair by long-range method or 
        !******** short-range method (a different code entirely). 
        T      = zeta*QA2;
        
        if (T .ge. Tcutoff) then;

temp = 0d0;

intcnt = intcnt + maxrgpa
           QA  = sqrt(QA2);    rQA  = 1d0/QA
           
           call calcderthph(dangdx,dangdy,dangdz,ang,QAx,QAy,QAz,QA,QA2,dangdth,dangdph)
           
           !*** Calculates inverse powers of QA
           rQ(0) = rQA;           do l=1,4;    rQ(l) = rQ(l-1)*rQA;              end do
         
k=0; do l=0,4; do m=-l,l; k=k+1;
angpart(1,k) = pf(l)*gf*rQ(l)*(0.5d0*dangdx(k) + ang(k)*(gamma*CDx - QAx*rQ(1)*(l+1)*0.5d0));
angpart(2,k) = pf(l)*gf*rQ(l)*(0.5d0*dangdy(k) + ang(k)*(gamma*CDy - QAy*rQ(1)*(l+1)*0.5d0));
angpart(3,k) = pf(l)*gf*rQ(l)*(0.5d0*dangdz(k) + ang(k)*(gamma*CDz - QAz*rQ(1)*(l+1)*0.5d0));
end do; end do;
        
temp = 0d0;
call dgemm('N','N',3,maxrgpa,25,1d0, angpart, &
       3, contnewsortRg(lenmodelcontnewsort(3,1)+1: lenmodelcontnewsort(3,maxrgpa)+25),&
       25,1d0, temp,3)

        do i=1,maxrgpa
                 bf1 = lenmodelcontnewsort(1,i)
                 bf2 = lenmodelcontnewsort(2,i)

     Gtot(bf1,bf2) = Gtot(bf1,bf2) + 2d0*sum(Ptot(bf3,bf4:bf4+2)*temp(1:3,i))
 
     Gtot(bf3,bf4:bf4+2) = Gtot(bf3,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp(1:3,i)
     
     Galpha(bf1,bf4:bf4+2) = Galpha(bf1,bf4:bf4+2) -Palpha(bf2,bf3)*temp(1:3,i) 
     Gbeta(bf1,bf4:bf4+2)  = Gbeta(bf1,bf4:bf4+2)  - Pbeta(bf2,bf3)*temp(1:3,i)
     Galpha(bf2,bf4:bf4+2) = Galpha(bf2,bf4:bf4+2) -Palpha(bf1,bf3)*temp(1:3,i)
     Gbeta(bf2,bf4:bf4+2)  = Gbeta(bf2,bf4:bf4+2)  - Pbeta(bf1,bf3)*temp(1:3,i)

     Galpha(bf2,bf3) = Galpha(bf2,bf3) - sum(Palpha(bf1,bf4:bf4+2)*temp(1:3,i))
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - sum(Pbeta(bf1,bf4:bf4+2)*temp(1:3,i))
     Galpha(bf1,bf3) = Galpha(bf1,bf3) - sum(Palpha(bf2,bf4:bf4+2)*temp(1:3,i))
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - sum(Pbeta(bf2,bf4:bf4+2)*temp(1:3,i))

        end do
        end if            !*** Ends selection of LR vs SR code
  end do           !*** Ends loop over sp shell pairs
    
  return
end subroutine calcLRgenRIspdigest

!***************************************************************

