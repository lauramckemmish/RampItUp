
!******************
!******************
!******************
!******************
!******************

subroutine calcSRgenRIspdigest(noggsp, ggsp, noatoms, atoms,  &
      maxrampdegree, Tcutoff,pts,at, newsortRg,lenmodelnewsort,maxrgpa, Omrnmatrix,& 
      Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,integrand,integrandsp,&
      ang,dangdx, dangdy, dangdz,dangdth,dangdph,intcnt,Gtot, & 
           newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     :: noatoms, noggsp, maxrampdegree, pts, at, nobasisfun
  real*8, intent(in)      :: atoms(1:5,1:noatoms), Tcutoff, ggsp(1:20,1:noggsp)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  integer                 :: i,j, m, cnt, ntot, k, atid, angm, warnswit(1:10), v, printopt
  integer                 :: lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, w, angk, angl, bf1, bf2, bf3, bf4, bb4,intcnt
  real*8                  :: Ax, Ay, Az, Qx, Qy, Qz, Cx, Cy, Cz, Dx, Dy, Dz, QA2, T, x, y, eeint
  real*8                  :: CDx, CDy, CDz, QAx, QAy, QAz, QA, tempPIsp, cpfs, cpfp, nc1, nc2
  real*8                  :: gamma, delta, zeta, sqrtzeta, rzeta ,pi, sqrtpi, rQA, gammaCDx, gammaCDy, gammaCDz
  real*8                  :: ang(1:25), dangdx(1:25),dangdy(1:25), dangdz(1:25), cpf
  real*8                  :: absc(1:2*pts+1), weights(1:2*pts+1), z, Omrnmatrix(1:2*pts+1,0:maxrampdegree)
  real*8                  :: Omrn(1:2*pts+1), r1(1:2*pts+1), r2(1:2*pts+1), r3(1:2*pts+1), r4(1:2*pts+1)
  real*8                  :: r6(1:2*pts+1), r8(1:2*pts+1), r10(1:2*pts+1),r12(1:2*pts+1), RIspmain, RIssmain
  real*8                  :: expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1)
  real*8                  :: integrand(1:2*pts+1,0:4), evenpow(1:2*pts+1), oddpow(1:2*pts+1), integrandsp(1:2*pts+1,0:4) 
  real*8                  :: newsortRg(1:4,lenmodelnewsort(3,1)+1:lenmodelnewsort(4,maxrgpa)), cpfmat(0:4)
  real*8                  :: totpart(1:3,1:25,0:maxrampdegree), temp(1:3)
  real*8                  :: dangdth(1:25), dangdph(1:25), Ptotbf1bf2
  integer :: lengthtwoelec, braketno, cnt2
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 

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

  warnswit = 0;    pi = dacos(-1d0);   sqrtpi = sqrt(pi); 
  
  call getquad(weights,absc,pts);  call evalallrk(r1,r2,r3,r4, r6, r8,absc,pts);
  r10=r8*r2;   r12 = r6*r6; 

!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
!!                 MAIN LOOP OVER SHELL QUARTETS
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!
           Ax      = atoms(3,at);    Ay      = atoms(4,at);    Az      = atoms(5,at)
  do j=1,noggsp
     bf3      = ggsp(1,j);      bf4      = ggsp(2,j)
     nc2      = ggsp(9,j);
     zeta     = ggsp(8,j);    
     gamma = ggsp(6,j);         delta = ggsp(7,j)
     CDx   = ggsp(11,j); CDy = ggsp(12,j); CDz = ggsp(13,j);
   sqrtzeta = ggsp(14,j);    z = ggsp(15,j); rzeta = ggsp(16,j)
     gammaCDx = gamma*CDx; gammaCDy = gamma*CDy; gammaCDz = gamma*CDz

           QAx     = ggsp(3,j) - Ax;          QAy     = ggsp(4,j) - Ay;          QAz     = ggsp(5,j) - Az;
           QA2 = QAx**2+QAy**2+QAz**2;     QA = sqrt(QA2);         rQA     = 1/QA;
           
           T     = zeta*QA2;
           
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
           !!          DECISION 2: T < 20
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!
           if (T < Tcutoff) then;

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
   call allerfexp(erfQmr,erfQpr,expQmr,expQpr,absc,pts,QA, sqrtzeta,zeta); 
   call calcr2lrall(integrand, z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, &
                   oddpow,evenpow)
   call calcspall(integrandsp, z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, &
         oddpow, evenpow)

              do v=1,2*pts+1;  
                  integrand(v,0:4) = integrand(v,0:4)*weights(v); 
                  integrandsp(v,0:4) = integrandsp(v,0:4)*weights(v);
               end do
 
do k=1,25
   dangdx(k) = 0.5d0* dangdx(k)
   dangdy(k) = 0.5d0* dangdy(k)
   dangdz(k) = 0.5d0* dangdz(k)
end do

do ntot = 0,maxrampdegree
 angk = 0; do angl=0,4; 
cpf = 16d0*pi**2/(2d0*angL+1d0)*rzeta*nc2
                       RIssmain = sum(integrand(1:2*pts+1,angl)*Omrnmatrix(1:2*pts+1,ntot))*cpf;
                       RIspmain = sum(integrandsp(1:2*pts+1,angl)*Omrnmatrix(1:2*pts+1,ntot))*0.5d0*cpf*rQA;
  do angm = -angl,angl; angk = angk+1;
  totpart(1,angk,ntot) = RIssmain*(gammaCDx*ang(angk)+ dangdx(angk)) + RIspmain*QAx*ang(angk)
  totpart(2,angk,ntot) = RIssmain*(gammaCDy*ang(angk)+ dangdy(angk)) + RIspmain*QAy*ang(angk)
  totpart(3,angk,ntot) = RIssmain*(gammaCDz*ang(angk)+ dangdz(angk)) + RIspmain*QAz*ang(angk)

end do; end do; end do;
          !*** Main loop over relevant ramp model components
              do i=1, maxrgpa
temp = 0d0;

                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    nc1  = newsortRg(1,w)
                    ntot = newsortRg(2,w)
                    angk  = newsortRg(3,w)
                    angl  = newsortRg(4,w)
 
temp = temp+ totpart(1:3,angk,ntot)*nc1
                 end do   
                bf1 = lenmodelnewsort(1,i);
                bf2 = lenmodelnewsort(2,i)

if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;        end if;
do cnt2=1,3; 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+cnt2-1,temp(cnt2),braketno)
end do; 
end if

     Gtot(bf3,bf4:bf4+2) = Gtot(bf3,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp
     Gtot(bf1,bf2)       = Gtot(bf1,bf2)       + 2d0*sum(Ptot(bf3,bf4:bf4+2)* temp) 

     Galpha(bf1,bf4:bf4+2) = Galpha(bf1,bf4:bf4+2) -Palpha(bf2,bf3)*temp 
     Gbeta(bf1,bf4:bf4+2)  = Gbeta(bf1,bf4:bf4+2)  - Pbeta(bf2,bf3)*temp
     Galpha(bf2,bf4:bf4+2) = Galpha(bf2,bf4:bf4+2) -Palpha(bf1,bf3)*temp 
     Gbeta(bf2,bf4:bf4+2)  = Gbeta(bf2,bf4:bf4+2)  - Pbeta(bf1,bf3)*temp

     Galpha(bf1,bf3) = Galpha(bf1,bf3) -sum(Palpha(bf2,bf4:bf4+2)*temp)
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - sum(Pbeta(bf2,bf4:bf4+2)*temp)
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -sum(Palpha(bf1,bf4:bf4+2)*temp)
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - sum(Pbeta(bf1,bf4:bf4+2)*temp)

              end do;           !** Ends the loop over the models with K angular momentum


                 
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


              do i=1, maxrgpa
                 temp = 0d0;
                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    nc1  = newsortRg(1,w)
                    ntot = newsortRg(2,w)
                    angk  = newsortRg(3,w)
  if (angk == 1) then; 
                    temp(1) = temp(1) + CDx*sum(integrand(1:2*pts+1,0)*Omrnmatrix(1:2*pts+1,ntot))*cpfs*nc1
                    temp(2) = temp(2) + CDy*sum(integrand(1:2*pts+1,0)*Omrnmatrix(1:2*pts+1,ntot))*cpfs*nc1
                    temp(3) = temp(3) + CDz*sum(integrand(1:2*pts+1,0)*Omrnmatrix(1:2*pts+1,ntot))*cpfs*nc1
  else if (angk == 2) then
                    temp(2) = temp(2) - sum(integrand(1:2*pts+1,1)*Omrnmatrix(1:2*pts+1,ntot))*cpfp*nc1
  else if (angk == 3) then; 
                    temp(3) = temp(3) + sum(integrand(1:2*pts+1,1)*Omrnmatrix(1:2*pts+1,ntot))*cpfp*nc1
  else if (angk == 4) then; 
                    temp(1) = temp(1) - sum(integrand(1:2*pts+1,1)*Omrnmatrix(1:2*pts+1,ntot))*cpfp*nc1
  end if
                 end do        
                       cnt = cnt + 1
                        bf1 = lenmodelnewsort(1,i);
                       bf2 = lenmodelnewsort(2,i)

if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;        end if;
do cnt2=1,3; 
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4+cnt2-1,temp(cnt2),braketno)
end do; 
end if

     Gtot(bf3,bf4:bf4+2) = Gtot(bf3,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp
     Gtot(bf1,bf2)       = Gtot(bf1,bf2)       + 2d0*sum(Ptot(bf3,bf4:bf4+2)* temp) 

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
 
  
  return
end subroutine calcSRgenRIspdigest



