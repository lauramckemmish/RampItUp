
!******************
!******************
!******************
!******************
!******************

subroutine calcSRgenRIspdigestoneRperatom(noggsp, ggsp, noatoms, atoms,  &
      maxrampdegree, Tcutoff,pts,at, newsortRg,lenmodelnewsort,maxrgpa, Omrnmatrix,& 
      Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,integrand,integrandsp,intcnt,Gtot,absc,weights, RIssmain, RIspmain, &
     r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow,totpart, newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     :: noatoms, noggsp, maxrampdegree, pts, at, nobasisfun
  real*8, intent(in)      :: atoms(1:5,1:noatoms), Tcutoff, ggsp(1:20,1:noggsp), Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: i,j, m, cnt, ntot, k, atid, angm, warnswit(1:10), v, printopt, n,ntotlist(0:4),angklist(0:4)
  integer                 :: lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, w, angk, angl, bf1, bf2, bf3, bf4, bb4,intcnt
  real*8                  :: newsortRg(1:4,lenmodelnewsort(3,1)+1:lenmodelnewsort(4,maxrgpa))
  real*8                  :: RIspmain(0:maxrampdegree,0:4), RIssmain(0:maxrampdegree,0:4)
  real*8                  :: integrand(1:2*pts+1,0:4), integrandsp(1:2*pts+1,0:4) , Omrnmatrix(1:2*pts+1,0:maxrampdegree)
  real*8                  :: absc(1:2*pts+1), weights(1:2*pts+1), z, inttlist(1:3,1:maxrgpa)
  real*8                  :: Ax, Ay, Az, QA2, T, x, y, eeint, eeintpart(1:3,0:4), Gtot(1:nobasisfun,1:nobasisfun)
  real*8                  :: CDx, CDy, CDz, QAx, QAy, QAz, QA, tempPIsp, cpfs, cpfp, nc1, nc2
  real*8                  :: gamma, delta, zeta, sqrtzeta, rzeta ,pi, sqrtpi, rQA, gammaCDx, gammaCDy, gammaCDz
  real*8                  :: ang(1:25), dangdx(1:25),dangdy(1:25), dangdz(1:25), cpf
  real*8                  :: dangdth(1:25), dangdph(1:25), Ptotbf1bf2, temp(1:3), cpfmat(0:4), acum(1:12), stor(1:11)
  real*8  :: r1(1:2*pts+1), r2(1:2*pts+1), r3(1:2*pts+1), r4(1:2*pts+1),r6(1:2*pts+1), r8(1:2*pts+1), r10(1:2*pts+1),&
      r12(1:2*pts+1),expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1),&
      evenpow(1:2*pts+1), oddpow(1:2*pts+1)


  real*8, dimension(:), allocatable :: conc1, conc0
  real*8:: totpart(1:3,0:maxrampdegree,1:25)   !** The real memory here is larger, but it should be fine as it is only in this function 

  integer :: lengthtwoelec,braketno,cnt2
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 
 allocate(conc1(0:maxrampdegree),conc0(0:maxrampdegree))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! genRIsp calculates integrals of the form (SS|sp) where the |sp) 
  ! is not concentric and large exponent. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "DEBUG.h"
printopt = 0
#ifdef DEBUG_PRINT
print *, "in SR (RI|sp)"
printopt = 1;
#endif

  bf1 = lenmodelnewsort(1,1)

  warnswit = 0;    pi = dacos(-1d0);   sqrtpi = sqrt(pi); 
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
!!                 MAIN LOOP OVER SHELL QUARTETS
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!
  Ax      = atoms(3,at);    Ay      = atoms(4,at);    Az      = atoms(5,at)
  do j=1,noggsp
     zeta     = ggsp(8,j);   
     QAx     = ggsp(3,j) - Ax;          QAy     = ggsp(4,j) - Ay;          QAz     = ggsp(5,j) - Az;
     QA2 = QAx**2+QAy**2+QAz**2;     QA = sqrt(QA2);         rQA     = 1/QA;
     T     = zeta*QA2;
           
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
           !!          DECISION 2: T < 20
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!
           if (T < Tcutoff) then;

     gamma = ggsp(6,j);         delta = ggsp(7,j)
     CDx   = ggsp(11,j); CDy = ggsp(12,j); CDz = ggsp(13,j);
   sqrtzeta = ggsp(14,j);    z = ggsp(15,j); rzeta = ggsp(16,j)
     gammaCDx = gamma*CDx; gammaCDy = gamma*CDy; gammaCDz = gamma*CDz
     bf3      = ggsp(1,j);      bf4      = ggsp(2,j);     nc2      = ggsp(9,j);

#ifdef INTCNTPRT 
intcnt = intcnt + maxrgpa
#endif

!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
              !!          DECISION 3: Concentric vs non-concentric
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!
              if (QA2 > 1d-15) then;
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
!!          NON-CONCENTRIC
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!
!***** Calculating the angular component!!
  call calcderthph(dangdx,dangdy,dangdz,ang,QAx,QAy,QAz,QA,QA2,dangdth,dangdph)


    expQmr = dexp(-zeta*(QA-absc)**2)
    expQpr = dexp(-zeta*(QA+absc)**2)
     erfQmr = derf(sqrtzeta*(QA-absc))
     erfQpr = derf(sqrtzeta*(QA+absc))

 call calcr2lrall(integrand(1:2*pts+1,0:4), z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, oddpow,evenpow)
 call calcspall(integrandsp(1:2*pts+1,0:4), z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, oddpow, evenpow)

 do angL=0,4
   integrand(1:2*pts+1,angL) = integrand(1:2*pts+1,angL)*weights; 
   integrandsp(1:2*pts+1,angL) = integrandsp(1:2*pts+1,angL)*weights;
 end do

call dscal(25,0.5d0,dangdx,1);  call dscal(25,0.5d0,dangdy,1);  call dscal(25,0.5d0,dangdz,1) 

call dgemm('T','N',maxrampdegree+1,5,2*pts+1,1d0,Omrnmatrix,2*pts+1,integrand,2*pts+1,0d0,RIssmain,maxrampdegree+1)
call dgemm('T','N',maxrampdegree+1,5,2*pts+1,1d0,Omrnmatrix,2*pts+1,integrandsp,2*pts+1,0d0,RIspmain,maxrampdegree+1)

do angL=0,4
  cpf = 16d0*pi**2/(2d0*angL+1d0)*rzeta*nc2
  RIssmain(0:maxrampdegree,angL) = RIssmain(0:maxrampdegree,angL)*cpf;
  RIspmain(0:maxrampdegree,angL) = RIspmain(0:maxrampdegree,angL)*0.5d0*cpf*rQA;
end do

 angk = 0; do angl=0,4; 
  do angm = -angl,angl; angk = angk+1;
  totpart(1,0:maxrampdegree,angk) = RIssmain(0:maxrampdegree,angL)*(gammaCDx*ang(angk)+ dangdx(angk)) &
                                          + RIspmain(0:maxrampdegree,angL)*QAx*ang(angk)
  totpart(2,0:maxrampdegree,angk) = RIssmain(0:maxrampdegree,angL)*(gammaCDy*ang(angk)+ dangdy(angk)) &
                                          + RIspmain(0:maxrampdegree,angL)*QAy*ang(angk)
  totpart(3,0:maxrampdegree,angk) = RIssmain(0:maxrampdegree,angL)*(gammaCDz*ang(angk)+ dangdz(angk)) &
                                          + RIspmain(0:maxrampdegree,angL)*QAz*ang(angk)
end do; end do; 
          !*** Main loop over relevant ramp model components
acum = 0d0;
stor(1:3) = Palpha(bf4:bf4+2,bf1)
stor(4:6) = Pbeta(bf4:bf4+2,bf1)
stor(7:9) = 2d0*Ptot(bf3,bf4:bf4+2)
stor(10) = Palpha(bf3,bf1)
stor(11) = Pbeta(bf3,bf1)

 do i=1, maxrgpa
     eeintpart = 0d0;
     do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)-6,5
  do k=0,4
       eeintpart(1:3,k) = eeintpart(1:3,k)+ totpart(1:3,int(newsortRg(2,w+k)),int(newsortRg(3,w+k)))*newsortRg(1,w+k)
  end do
     end do     
v = w;
do k=1,3
temp(k) = sum(eeintpart(k,0:4))
end do
     do w= v, lenmodelnewsort(4,i)
       ntot = newsortRg(2,w);
       angk = newsortRg(3,w);
       temp = temp+ totpart(1:3,ntot,angk)*newsortRg(1,w)
     end do      
inttlist(1:3,i) = temp;
end do
 do i=1,maxrgpa
temp = inttlist(1:3,i)
     bf2 = lenmodelnewsort(2,i)
if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;        end if;
if (i .ne. j .and. braketno == 1) then
do cnt2=1,3; 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+cnt2-1,2*temp(cnt2),braketno)
end do; 
else
do cnt2=1,3; 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+cnt2-1,temp(cnt2),braketno)
end do;
end if 
end if

     acum(1:3) = acum(1:3) + 2d0*Ptot(bf2,bf1)*temp
     acum(4:6) = acum(4:6) -Palpha(bf2,bf3)*temp 
     acum(7:9) = acum(7:9)- Pbeta(bf2,bf3)*temp
     acum(10) = acum(10) -sum(Palpha(bf2,bf4:bf4+2)*temp)
     acum(11) = acum(11)  - sum(Pbeta(bf2,bf4:bf4+2)*temp)

     Galpha(bf2,bf3) = Galpha(bf2,bf3) -sum(stor(1:3)*temp)
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - sum(stor(4:6)*temp)
     Gtot(bf2,bf1)   = Gtot(bf2,bf1) + sum(stor(7:9)* temp) 

     Galpha(bf4:bf4+2,bf2) = Galpha(bf4:bf4+2,bf2) -stor(10)*temp 
     Gbeta(bf4:bf4+2,bf2)  = Gbeta(bf4:bf4+2,bf2)  -stor(11)*temp

  end do;           !** Ends the loop over the models 


 Gtot(bf3,bf4:bf4+2) = Gtot(bf3,bf4:bf4+2) + acum(1:3)
 Galpha(bf1,bf4:bf4+2) = Galpha(bf1,bf4:bf4+2)+acum(4:6) 
 Gbeta(bf1,bf4:bf4+2)  = Gbeta(bf1,bf4:bf4+2) +acum(7:9)
 Galpha(bf1,bf3) = Galpha(bf1,bf3) +acum(10)
 Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  +acum(11)

                 
              else 
                 
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
                 !!         CONCENTRIC CASE
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!
                   do v=1,2*pts+1
                      integrand(v,0) = weights(v)*r1(v)*erf(r1(v)*sqrtzeta);
                      integrand(v,1) = weights(v)*r1(v)*(-2*r1(v)*sqrtzeta*exp(-r2(v)*zeta) &
                          + 1.7724538509055160273d0*erf(r1(v)*sqrtzeta)   );
                   end do
              !*** Main loop over relevant ramp model components
              cpfs = 19.739208802178717238d0*z**3*gamma*rzeta*nc2;
              cpfp = 3.2148756679069161256d0*z**5*nc2;   

do ntot=0,maxrampdegree
   conc0(ntot) = sum(integrand(1:2*pts+1,0)*Omrnmatrix(1:2*pts+1,ntot))*cpfs
   conc1(ntot) = sum(integrand(1:2*pts+1,1)*Omrnmatrix(1:2*pts+1,ntot))*cpfp
end do
              do i=1, maxrgpa
                 temp = 0d0;
                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    angk  = newsortRg(3,w)
  if (angk == 1) then; 
                    ntot = newsortRg(2,w)
                    temp(1) = temp(1) + CDx*conc0(ntot)*newsortRg(1,w)
                    temp(2) = temp(2) + CDy*conc0(ntot)*newsortRg(1,w)
                    temp(3) = temp(3) + CDz*conc0(ntot)*newsortRg(1,w)
  else if (angk == 2) then
                    ntot = newsortRg(2,w)
                    temp(2) = temp(2) - conc1(ntot)*newsortRg(1,w)
  else if (angk == 3) then; 
                    ntot = newsortRg(2,w)
                    temp(3) = temp(3) + conc1(ntot)*newsortRg(1,w)
  else if (angk == 4) then; 
                    ntot = newsortRg(2,w)
                    temp(1) = temp(1) - conc1(ntot)*newsortRg(1,w)
  else
    exit
  end if
                 end do        

                       cnt = cnt + 1
                       bf2 = lenmodelnewsort(2,i)


if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;        end if;
if (i .ne. j .and. braketno == 1) then;
do cnt2=1,3; 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+cnt2-1,2*temp(cnt2),braketno)
end do; 
else
do cnt2=1,3; 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+cnt2-1,temp(cnt2),braketno)
end do; 
end if
end if

     Gtot(bf3,bf4:bf4+2) = Gtot(bf3,bf4:bf4+2) + 2d0*Ptot(bf2,bf1)*temp
     Gtot(bf2,bf1)       = Gtot(bf2,bf1)       + 2d0*sum(Ptot(bf3,bf4:bf4+2)*temp) 

     Galpha(bf1,bf4:bf4+2) = Galpha(bf1,bf4:bf4+2) -Palpha(bf2,bf3)*temp 
     Gbeta(bf1,bf4:bf4+2)  = Gbeta(bf1,bf4:bf4+2)  - Pbeta(bf2,bf3)*temp
     Galpha(bf2,bf4:bf4+2) = Galpha(bf2,bf4:bf4+2) -Palpha(bf1,bf3)*temp 
     Gbeta(bf2,bf4:bf4+2)  = Gbeta(bf2,bf4:bf4+2)  - Pbeta(bf1,bf3)*temp

     Galpha(bf1,bf3) = Galpha(bf1,bf3) -sum(Palpha(bf2,bf4:bf4+2)*temp)
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - sum(Pbeta(bf2,bf4:bf4+2)*temp)
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -sum(Palpha(bf1,bf4:bf4+2)*temp)
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - sum(Pbeta(bf1,bf4:bf4+2)*temp)


              end do;           !** Ends the loop over the models with K angular momentum

              end if;      !*** Ends selection concentric vs non-concentric
           end if;   !*** End selection of short-range shell quartets only (i.e. overlapping shell pairs)
  end do      !*** Ends loop over all sp shell pairs
 
  
deallocate(conc1,conc0)
  return
end subroutine calcSRgenRIspdigestoneRperatom



