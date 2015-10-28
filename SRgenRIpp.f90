
subroutine calcSRgenRIppdigestoneRperatom(noggpp, ggpp,  noatoms, atoms, &
     maxrampdegree, Tcutoff, pts,atid, newsortRg,lenmodelnewsort,maxrgpa, Omrnmatrix, & 
      Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun, integrand, integrandsp, integrandpp,intcnt, &
      Gtot,absc,weights, RIssmain, RIspmain, RIppmain, &
     r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow,totpart, & 
           newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     :: noatoms, noggpp, pts, maxrampdegree,atid, nobasisfun
  real*8, intent(in)      :: atoms(1:5,1:noatoms),ggpp(1:20,1:noggpp),Tcutoff
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  integer                 :: i,j, m, ntot, k, angk, angl, angm, warnswit(1:10), v, bf1, bf2, bf3, bf4, bb3, bb4
  integer                 :: lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, w, ksp, lsp,intcnt, printopt, ttt
  real*8                  :: Ax, Ay, Az, P, Qx, Qy, Qz, QA, QA2, Cx, Cy, Cz, Dx, Dy, Dz, tempSIss, tempPIsp
  real*8                  :: gamma, delta, zeta, sqrtzeta, pi, sqrtpi, pftot,  nc1, nc2, temp2
  real*8                  :: T, a, tempRIpp,  rQA, QAx, QAy, QAz, CDx, CDy, CDz, pi3r2, sq3r2, rzeta, D0Ipzpz, z, cpf, cpf2
  real*8                  :: ang(1:25), dangdx(1:25), dangdy(1:25), dangdz(1:25)
  real*8                  :: d2angdx2(1:25), d2angdxy(1:25), d2angdxz(1:25), d2angdy2(1:25), d2angdyz(1:25), d2angdz2(1:25)
  real*8                  :: RIspxmain, RIspymain, RIspzmain, absc(1:2*pts+1), weights(1:2*pts+1)
  real*8                  :: Omrn(1:2*pts+1), r1(1:2*pts+1), r2(1:2*pts+1), r3(1:2*pts+1), r4(1:2*pts+1)
  real*8                  ::  RIssmain(0:maxrampdegree,0:4),RIspmain(0:maxrampdegree,0:4),RIppmain(0:maxrampdegree,0:4)
  real*8                  :: r6(1:2*pts+1), r8(1:2*pts+1),r10(1:2*pts+1), r12(1:2*pts+1), r5(1:2*pts+1), r7(1:2*pts+1)
  real*8                  :: expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1), r14(1:2*pts+1)
  real*8                  :: integrand(1:2*pts+1,0:4), evenpow(1:2*pts+1), oddpow(1:2*pts+1)
  real*8                  :: integrandsp(1:2*pts+1,0:4), integrandpp(1:2*pts+1,0:4)
  real*8                  :: Omrnmatrix(1:2*pts+1,0:maxrampdegree), Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz
  real*8                  :: newsortRg(1:4,lenmodelnewsort(3,1)+1:lenmodelnewsort(4,maxrgpa)), cpfmat(0:4) 
  real*8                  :: Fourgade, Tgmd, rQA2, totpart(1:9,0:maxrampdegree,1:25), cpfs, cpfs2, cpfp, cpfd
  real*8                  :: dangdth(1:25), dangdph(1:25), intt, x, y, temp(1:9)
  real*8                  :: d2angdth2(1:25), d2angdph2(1:25),d2angdthph(1:25)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8      :: inttlist(1:9,1:maxrgpa)
  real*8      :: acum(1:21), scale, conc0a(0:maxrampdegree), conc0b(0:maxrampdegree), conc1(0:maxrampdegree), conc2(0:maxrampdegree)
  integer :: lengthtwoelec, cnt, braketno, cnt2
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! genSSss calculates integrals of the form (SS|pp) where the |pp) 
  ! is not concentric and large exponent. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  warnswit = 0;    

  pi = dacos(-1d0);  sqrtpi = sqrt(pi); pi3r2 = pi**(3d0/2d0);   sq3r2 = sqrt(3d0)/2d0;
  
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
  !!                  MAIN LOOP OVER SHELL QUARTETS
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
  do j=1,noggpp
     bf3      = ggpp(1,j);   bf4      = ggpp(2,j);
     nc2      = ggpp(9,j)
     zeta     = ggpp(8,j);    sqrtzeta = sqrt(zeta);    z = 1d0/sqrtzeta;     rzeta    = 1/zeta;     
     Qx       = ggpp(3,j);    Qy       = ggpp(4,j);    Qz       = ggpp(5,j)
        Ax      = atoms(3,atid);    Ay      = atoms(4,atid);    Az      = atoms(5,atid) ;
        QAx     = Qx - Ax;          QAy     = Qy -Ay;           QAz     = Qz-Az;
        QA2 =(QAx**2+QAy**2+QAz**2);    
        
        T = zeta*QA2; 
        
        !*******!!!!!!*********!!!!!!!!**********!!!!!!!!!*********!!!!!!!!!***
        !   DECISION 2: T < 20
        !*******!!!!!!*********!!!!!!!!**********!!!!!!!!!*********!!!!!!!!!***
        if (T .lt. Tcutoff) then
           
           !*******!!!!!!*********!!!!!!!!**********!!!!!!!!!*********!!!!!!!!!***
           !   DECISION 3: Concentric or non-concentric
           !*******!!!!!!*********!!!!!!!!**********!!!!!!!!!*********!!!!!!!!!***
           
intcnt = intcnt + maxrgpa
     gamma    = ggpp(6,j);    delta    = ggpp(7,j)
     CDx       = ggpp(11,j);        CDy = ggpp(12,j);   CDz = ggpp(13,j)
     Fourgade = 4*gamma*delta;      Tgmd =2*(gamma-delta)
 QA = sqrt(QA2); rQA = 1d0/QA;  rQA2 = rQA**2;

           if (QA2 > 1d-15) then;
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
              !!              NON-CONCENTRIC CODE
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
              
              call calcder2thph(ang,dangdx,dangdy,dangdz,d2angdx2,d2angdxy,d2angdxz,d2angdy2,d2angdyz,&
                   d2angdz2,QAx,QAy,QAz,QA2,dangdth,dangdph,d2angdth2,d2angdph2,d2angdthph)

    erfQmr = derf(sqrtzeta*(QA-absc))
    erfQpr = derf(sqrtzeta*(QA+absc))
    expQmr = dexp(-zeta*(QA-absc)**2)
    expQpr = dexp(-zeta*(QA+absc)**2)

call calcssspppall(integrand,integrandsp,integrandpp(1:2*pts+1,0:4), z, &
QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, oddpow, evenpow)

              do v=1,2*pts+1;
                integrand(v,0:4)   = integrand(v,0:4)  *weights(v);
                integrandsp(v,0:4) = integrandsp(v,0:4)*weights(v);
                integrandpp(v,0:4) = integrandpp(v,0:4)*weights(v);
              end do

call dgemm('T','N',maxrampdegree+1,5,2*pts+1,1d0, Omrnmatrix,2*pts+1, integrand,2*pts+1,0d0,RIssmain,maxrampdegree+1)
call dgemm('T','N',maxrampdegree+1,5,2*pts+1,1d0,Omrnmatrix,2*pts+1, integrandsp,2*pts+1,0d0,RIspmain,maxrampdegree+1)
call dgemm('T','N',maxrampdegree+1,5,2*pts+1,1d0, Omrnmatrix,2*pts+1, integrandpp,2*pts+1,0d0,RIppmain,maxrampdegree+1)

do angL = 0,4
 scale = (157.91367041742973790d0/(2*angL+1d0))*0.25d0*rzeta**2*nc2;
     RIssmain(0:maxrampdegree,angL) = RIssmain(0:maxrampdegree,angL)*scale
     RIspmain(0:maxrampdegree,angL) = RIspmain(0:maxrampdegree,angL)* scale
     RIppmain(0:maxrampdegree,angL) = RIppmain(0:maxrampdegree,angL)* scale
end do


 angk = 0; do angL = 0,1; 
do angm=-angL,angL; angk=angk+1; 

 do ntot = 4,maxrampdegree

 totpart(1,ntot,angk) =   Tgmd*CDx*(ang(angk)*QAx*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*dangdx(angk)) + &
 (2*zeta - Fourgade*CDx**2)*ang(angk)*RIssmain(ntot,angL) + ang(angk)*RIspmain(ntot,angL)*(rQA-QAx**2*rQA**3) + &
 QAx**2*rQA2*ang(angk)*RIppmain(ntot,angL) + 2*dangdx(angk)*QAx*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*d2angdx2(angk)
                 
 totpart(2,ntot,angk) =   2*(gamma*CDy*QAx-delta*CDx*QAy)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDy*dangdx(angk) - delta*CDx*dangdy(angk)) &
 - Fourgade*CDx*CDy*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAy*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) +&
 RIspmain(ntot,angL)*rQA*(QAx*dangdy(angk)+QAy*dangdx(angk)) + &
 RIssmain(ntot,angL)*d2angdxy(angk)
                 
 totpart(3,ntot,angk) =   2*(gamma*CDz*QAx-delta*CDx*QAz)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDz*dangdx(angk) - delta*CDx*dangdz(angk)) &
 - Fourgade*CDx*CDz*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAz*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) +&
 RIspmain(ntot,angL)*rQA*(QAx*dangdz(angk)+QAz*dangdx(angk)) + &
 RIssmain(ntot,angL)*d2angdxz(angk)
                 
 totpart(4,ntot,angk) = 2*(gamma*CDx*QAy-delta*CDy*QAx)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDx*dangdy(angk) - delta*CDy*dangdx(angk)) &
 - Fourgade*CDx*CDy*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAy*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) + &
RIspmain(ntot,angL)*rQA*(QAy*dangdx(angk)+QAx*dangdy(angk)) + &
 RIssmain(ntot,angL)*d2angdxy(angk)
                 
 totpart(5,ntot,angk) =             Tgmd*CDy*(ang(angk)*QAy*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*dangdy(angk)) + &
 2*(zeta - 2*gamma*delta*CDy**2)*ang(angk)*RIssmain(ntot,angL) + ang(angk)*RIspmain(ntot,angL)*(rQA-QAy**2*rQA**3) + &
 QAy**2*rQA2*ang(angk)*RIppmain(ntot,angL) + 2*dangdy(angk)*QAy*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*d2angdy2(angk);
                 
 totpart(6,ntot,angk) =  +2*(gamma*CDz*QAy-delta*CDy*QAz)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDz*dangdy(angk) - delta*CDy*dangdz(angk))&
 - Fourgade*CDz*CDy*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAz*QAy*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) + &
RIspmain(ntot,angL)*rQA*(QAy*dangdz(angk)+QAz*dangdy(angk)) + &
 RIssmain(ntot,angL)*d2angdyz(angk)
                 
 totpart(7,ntot,angk) =  +2*(gamma*CDx*QAz-delta*CDz*QAx)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDx*dangdz(angk) - delta*CDz*dangdx(angk))&
 - Fourgade*CDx*CDz*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAz*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) + &
RIspmain(ntot,angL)*rQA*(QAx*dangdz(angk)+QAz*dangdx(angk)) + &
 RIssmain(ntot,angL)*d2angdxz(angk)
                 
 totpart(8,ntot,angk) =  +2*(gamma*CDy*QAz-delta*CDz*QAy)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDy*dangdz(angk) - delta*CDz*dangdy(angk)) &
 - Fourgade*CDy*CDz*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAy*QAz*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) &
+ RIspmain(ntot,angL)*rQA*(QAy*dangdz(angk)+QAz*dangdy(angk)) + &
 RIssmain(ntot,angL)*d2angdyz(angk)
                 
 totpart(9,ntot,angk) =         Tgmd*CDz*(ang(angk)*QAz*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*dangdz(angk)) + &
 2*(zeta - 2*gamma*delta*CDz**2)*ang(angk)*RIssmain(ntot,angL) + ang(angk)*RIspmain(ntot,angL)*(rQA-QAz**2*rQA**3) + &
 QAz**2*rQA2*ang(angk)*RIppmain(ntot,angL) + 2*dangdz(angk)*QAz*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*d2angdz2(angk);

end do; end do
end do

do angL = 2,4; 
do angm=-angL,angL; angk=angk+1; 

 do ntot = 4,30

 totpart(1,ntot,angk) =   Tgmd*CDx*(ang(angk)*QAx*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*dangdx(angk)) + &
 (2*zeta - Fourgade*CDx**2)*ang(angk)*RIssmain(ntot,angL) + ang(angk)*RIspmain(ntot,angL)*(rQA-QAx**2*rQA**3) + &
 QAx**2*rQA2*ang(angk)*RIppmain(ntot,angL) + 2*dangdx(angk)*QAx*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*d2angdx2(angk)
                 
 totpart(2,ntot,angk) =   2*(gamma*CDy*QAx-delta*CDx*QAy)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDy*dangdx(angk) - delta*CDx*dangdy(angk)) &
 - Fourgade*CDx*CDy*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAy*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) +&
 RIspmain(ntot,angL)*rQA*(QAx*dangdy(angk)+QAy*dangdx(angk)) + &
 RIssmain(ntot,angL)*d2angdxy(angk)
                 
 totpart(3,ntot,angk) =   2*(gamma*CDz*QAx-delta*CDx*QAz)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDz*dangdx(angk) - delta*CDx*dangdz(angk)) &
 - Fourgade*CDx*CDz*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAz*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) +&
 RIspmain(ntot,angL)*rQA*(QAx*dangdz(angk)+QAz*dangdx(angk)) + &
 RIssmain(ntot,angL)*d2angdxz(angk)
                 
 totpart(4,ntot,angk) = 2*(gamma*CDx*QAy-delta*CDy*QAx)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDx*dangdy(angk) - delta*CDy*dangdx(angk)) &
 - Fourgade*CDx*CDy*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAy*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) + &
RIspmain(ntot,angL)*rQA*(QAy*dangdx(angk)+QAx*dangdy(angk)) + &
 RIssmain(ntot,angL)*d2angdxy(angk)
                 
 totpart(5,ntot,angk) =             Tgmd*CDy*(ang(angk)*QAy*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*dangdy(angk)) + &
 2*(zeta - 2*gamma*delta*CDy**2)*ang(angk)*RIssmain(ntot,angL) + ang(angk)*RIspmain(ntot,angL)*(rQA-QAy**2*rQA**3) + &
 QAy**2*rQA2*ang(angk)*RIppmain(ntot,angL) + 2*dangdy(angk)*QAy*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*d2angdy2(angk);
                 
 totpart(6,ntot,angk) =  +2*(gamma*CDz*QAy-delta*CDy*QAz)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDz*dangdy(angk) - delta*CDy*dangdz(angk))&
 - Fourgade*CDz*CDy*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAz*QAy*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) + &
RIspmain(ntot,angL)*rQA*(QAy*dangdz(angk)+QAz*dangdy(angk)) + &
 RIssmain(ntot,angL)*d2angdyz(angk)
                 
 totpart(7,ntot,angk) =  +2*(gamma*CDx*QAz-delta*CDz*QAx)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDx*dangdz(angk) - delta*CDz*dangdx(angk))&
 - Fourgade*CDx*CDz*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAx*QAz*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) + &
RIspmain(ntot,angL)*rQA*(QAx*dangdz(angk)+QAz*dangdx(angk)) + &
 RIssmain(ntot,angL)*d2angdxz(angk)
                 
 totpart(8,ntot,angk) =  +2*(gamma*CDy*QAz-delta*CDz*QAy)*ang(angk)*RIspmain(ntot,angL)*rQA + &
 2*RIssmain(ntot,angL)*(gamma*CDy*dangdz(angk) - delta*CDz*dangdy(angk)) &
 - Fourgade*CDy*CDz*RIssmain(ntot,angL)*ang(angk) + &
 ang(angk)*QAy*QAz*rQA2*(RIppmain(ntot,angL)-RIspmain(ntot,angL)*rQA) &
+ RIspmain(ntot,angL)*rQA*(QAy*dangdz(angk)+QAz*dangdy(angk)) + &
 RIssmain(ntot,angL)*d2angdyz(angk)
                 
 totpart(9,ntot,angk) =         Tgmd*CDz*(ang(angk)*QAz*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*dangdz(angk)) + &
 2*(zeta - 2*gamma*delta*CDz**2)*ang(angk)*RIssmain(ntot,angL) + ang(angk)*RIspmain(ntot,angL)*(rQA-QAz**2*rQA**3) + &
 QAz**2*rQA2*ang(angk)*RIppmain(ntot,angL) + 2*dangdz(angk)*QAz*rQA*RIspmain(ntot,angL) + RIssmain(ntot,angL)*d2angdz2(angk);

end do; end do
end do

                bf1 = lenmodelnewsort(1,1); acum = 0d0;
          !*** Main loop over relevant ramp model components
              do i=1, maxrgpa
 temp = 0d0;
                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    temp =  temp + totpart(1:9,int(newsortRg(2,w)),int(newsortRg(3,w)))*newsortRg(1,w)
                 end do

inttlist(1:9,i) = temp;
end do
 do i=1,maxrgpa
temp = inttlist(1:9,i)
                bf2 = lenmodelnewsort(2,i)

if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if;
do cnt=1,3;
 do cnt2 = 1,3;
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3+cnt-1,bf4+cnt2-1,temp(3*cnt+cnt2-3),braketno)
 end do
end do
end if

     Gtot(bf2,bf1) = Gtot(bf2,bf1) + 2d0*sum(Ptot(bf4:bf4+2,bf3)*temp(1:3))  & 
          + 2d0*sum(Ptot(bf4:bf4+2,bf3+1)*temp(4:6)) + 2d0*sum(Ptot(bf4:bf4+2,bf3+2)*temp(7:9))

    acum(1) = acum(1) -sum(Palpha(bf2,bf4:bf4+2)*temp(1:3))
    acum(2) = acum(2) -sum(Pbeta(bf2,bf4:bf4+2) *temp(1:3))
    acum(3) = acum(3) -sum(Palpha(bf2,bf4:bf4+2)*temp(4:6))
    acum(4) = acum(4) -sum( Pbeta(bf2,bf4:bf4+2)*temp(4:6))
    acum(5) = acum(5) -sum(Palpha(bf2,bf4:bf4+2)*temp(7:9))
    acum(6) = acum(6) -sum( Pbeta(bf2,bf4:bf4+2)*temp(7:9))

     acum( 7: 9) = acum( 7: 9) + 2d0*Ptot(bf2,bf1)*temp(1:3)
     acum(10:12) = acum(10:12) + 2d0*Ptot(bf2,bf1)*temp(4:6)
     acum(13:15) = acum(13:15) + 2d0*Ptot(bf2,bf1)*temp(7:9)

     acum(16:18) = acum(16:18) -Palpha(bf2,bf3)*temp(1:3) -Palpha(bf2,bf3+1)*temp(4:6) -Palpha(bf2,bf3+2)*temp(7:9) 
     acum(19:21) = acum(19:21) - Pbeta(bf2,bf3)*temp(1:3) - Pbeta(bf2,bf3+1)*temp(4:6) - Pbeta(bf2,bf3+2)*temp(7:9)

     Galpha(bf3,bf2)   = Galpha(bf3,bf2)   -sum(Palpha(bf1,bf4:bf4+2)*temp(1:3))
     Gbeta(bf3,bf2)    = Gbeta(bf3,bf2)    -sum(Pbeta(bf1,bf4:bf4+2) *temp(1:3))
     Galpha(bf3+1,bf2) = Galpha(bf3+1,bf2) -sum(Palpha(bf1,bf4:bf4+2)*temp(4:6))
     Gbeta(bf3+1,bf2)  = Gbeta(bf3+1,bf2)  -sum(Pbeta(bf1,bf4:bf4+2) *temp(4:6))
     Galpha(bf3+2,bf2) = Galpha(bf3+2,bf2) -sum(Palpha(bf1,bf4:bf4+2)*temp(7:9))
     Gbeta(bf3+2,bf2)  = Gbeta(bf3+2,bf2)  -sum(Pbeta(bf1,bf4:bf4+2) *temp(7:9))


     Galpha(bf4:bf4+2,bf2) = Galpha(bf4:bf4+2,bf2) &
         -Palpha(bf1,bf3)*temp(1:3) -Palpha(bf1,bf3+1)*temp(4:6) -Palpha(bf1,bf3+2)*temp(7:9) 
     Gbeta(bf4:bf4+2,bf2)  = Gbeta(bf4:bf4+2,bf2)  &
         - Pbeta(bf1,bf3)*temp(1:3) - Pbeta(bf1,bf3+1)*temp(4:6) - Pbeta(bf1,bf3+2)*temp(7:9)


              end do
           
     Galpha(bf1,bf3)   = Galpha(bf1,bf3)   +acum(1)
     Gbeta(bf1,bf3)    = Gbeta(bf1,bf3)    +acum(2)
     Galpha(bf1,bf3+1) = Galpha(bf1,bf3+1) +acum(3)
     Gbeta(bf1,bf3+1)  = Gbeta(bf1,bf3+1)  +acum(4)
     Galpha(bf1,bf3+2) = Galpha(bf1,bf3+2) +acum(5)
     Gbeta(bf1,bf3+2)  = Gbeta(bf1,bf3+2)  +acum(6)

     Gtot(bf3  ,bf4:bf4+2) = Gtot(bf3  ,bf4:bf4+2) + acum(7:9)
     Gtot(bf3+1,bf4:bf4+2) = Gtot(bf3+1,bf4:bf4+2) + acum(10:12)
     Gtot(bf3+2,bf4:bf4+2) = Gtot(bf3+2,bf4:bf4+2) + acum(13:15)

     Galpha(bf1,bf4:bf4+2) = Galpha(bf1,bf4:bf4+2) + acum(16:18) 
     Gbeta(bf1,bf4:bf4+2)  = Gbeta(bf1,bf4:bf4+2)  + acum(19:21)


           
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
           !!                 CONCENTRIC CASE
!!!!!!!!!*************!!!!!!!!!!!!!********************!!!!!!!!!!!!!!!***************!!!!!!!!
        else if (QA2 .le. 1d-10) then;
 
           do v=1,2*pts+1
              integrand(v,0) = weights(v)*r1(v)*erf(r1(v)*sqrtzeta);
              integrand(v,1) = weights(v)*r1(v)*(-2*r1(v)*sqrtzeta*exp(-r2(v)*zeta) &
                          + 1.7724538509055160273d0*erf(r1(v)*sqrtzeta)   );
              integrand(v,2) = weights(v)*r1(v)*(3d0*sqrtpi*erf(r1(v)*sqrtzeta) &
                 -2*sqrtzeta*exp(-r2(v)*zeta)*r1(v)*(3+2*r2(v)*zeta));
              integrand(v,3) = weights(v)*r1(v)*(sqrtpi*erf(r1(v)*sqrtzeta) &
                 - (2d0/3d0)*r1(v)*sqrtzeta*exp(-r2(v)*zeta))
           end do

           cpfs = 2*pi**2*z**3*nc2*gamma*delta/zeta**2;   
           cpfs2 = pi**(3d0/2d0)*z**5*nc2;

           cpfp = 6.4297513358138322513d0*z**9*gamma*delta*nc2; 
           cpfd = 0.83007732812898906760d0*z**7*nc2;

do ntot = 0,maxrampdegree
 conc0a(ntot) = sum(integrand(1:2*pts+1,0)*Omrnmatrix(1:2*pts+1,ntot))*cpfs
 conc0b(ntot) = sum(integrand(1:2*pts+1,3)*Omrnmatrix(1:2*pts+1,ntot))*cpfs2
 conc1(ntot) = sum(integrand(1:2*pts+1,1)*Omrnmatrix(1:2*pts+1,ntot))*cpfp
 conc2(ntot)= sum(integrand(1:2*pts+1,2)*Omrnmatrix(1:2*pts+1,ntot))*cpfd
end do
             do i=1, maxrgpa
                 Rxx = 0d0; Rxy = 0d0; Rxz = 0d0;
                 Ryx = 0d0; Ryy = 0d0; Ryz = 0d0;
                 Rzx = 0d0; Rzy = 0d0; Rzz = 0d0;

                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    nc1  = newsortRg(1,w)
                    ntot = newsortRg(2,w)
                    angk  = newsortRg(3,w)
                    angl  = newsortRg(4,w)
  if (angk == 1) then; 
              tempSIss = conc0a(ntot)*nc1;
              tempRIpp = conc0b(ntot)*nc1;

                    Rxx = Rxx + tempRIpp - CDx**2*tempSIss; 
                    Rxy = Rxy - CDx*CDy*tempSIss
                    Rxz = Rxz - CDx*CDz*tempSIss
                    Ryx = Ryx - CDx*CDy*tempSIss
                    Ryy = Ryy + tempRIpp - CDy**2*tempSIss
                    Ryz = Ryz - CDz*CDy*tempSIss
                    Rzx = Rzx - CDz*CDx*tempSIss
                    Rzy = Rzy - CDz*CDy*tempSIss
                    Rzz = Rzz + tempRIpp - CDz**2*tempSIss
 else if (angk == 2) then;

                 tempPIsp = conc1(ntot)*nc1;
 
 Rxy = Rxy + CDx*tempPIsp
 Ryx = Ryx - CDx*tempPIsp
 Ryz = Ryz - CDz*tempPIsp
 Rzy = Rzy + CDz*tempPIsp

 else if (angk == 3) then;

                 tempPIsp = conc1(ntot)*nc1;
 
  Rzx = Rzx + CDx*tempPIsp;
  Rxz = Rxz - CDx*tempPIsp;
  Ryz = Ryz - CDy*tempPIsp;
  Rzy = Rzy + CDy*tempPIsp;


 else if (angk == 4) then; 
 
                 tempPIsp = conc1(ntot)*nc1;
 
 Rxy = Rxy - CDy*tempPIsp
 Rxz = Rxz - CDz*tempPIsp
 Ryx = Ryx + CDy*tempPIsp
 Rzx = Rzx + CDz*tempPIsp

else if (angk == 5) then

              D0Ipzpz = conc2(ntot)*nc1;

Rxy = Rxy + sq3r2* D0Ipzpz
Ryx = Ryx + sq3r2* D0Ipzpz

else if (angk == 6) then

              D0Ipzpz = conc2(ntot)*nc1;

Rzy = Rzy - sq3r2* D0Ipzpz
Ryz = Ryz - sq3r2* D0Ipzpz

else if (angk == 7) then

              D0Ipzpz = conc2(ntot)*nc1;

Rxx = Rxx - 0.5d0* D0Ipzpz
Ryy = Ryy - 0.5d0* D0Ipzpz
Rzz = Rzz +  D0Ipzpz


else if (angk == 8) then

              D0Ipzpz = conc2(ntot)*nc1;


Rzx = Rzx - sq3r2* D0Ipzpz
Rxz = Rxz - sq3r2* D0Ipzpz


else if (angk == 9) then

              D0Ipzpz = conc2(ntot)*nc1;
Rxx = Rxx + sq3r2* D0Ipzpz
Ryy = Ryy - sq3r2* D0Ipzpz
else
    exit

  end if
                 end do        
                        bf1 = lenmodelnewsort(1,i);
                       bf2 = lenmodelnewsort(2,i)

temp(1) = Rxx; temp(2) = Rxy; temp(3) = Rxz; temp(4) = Ryx; temp(5) = Ryy; temp(6) = Ryz; 
temp(7) = Rzx; temp(8) = Rzy; temp(9) = Rzz;


if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;        end if;
do cnt=1,3; do cnt2 = 1,3;
  call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3+cnt-1,bf4+cnt2-1,temp(3*cnt+cnt2-3),braketno)
end do; end do;
end if

     Gtot(bf1,bf2) = Gtot(bf1,bf2) + 2d0*sum(Ptot(bf3,bf4:bf4+2)*temp(1:3))  & 
          + 2d0*sum(Ptot(bf3+1,bf4:bf4+2)*temp(4:6)) + 2d0*sum(Ptot(bf3+2,bf4:bf4+2)*temp(7:9))

     Gtot(bf3  ,bf4:bf4+2) = Gtot(bf3  ,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp(1:3)
     Gtot(bf3+1,bf4:bf4+2) = Gtot(bf3+1,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp(4:6)
     Gtot(bf3+2,bf4:bf4+2) = Gtot(bf3+2,bf4:bf4+2) + 2d0*Ptot(bf1,bf2)*temp(7:9)

     Galpha(bf1,bf4:bf4+2) = Galpha(bf1,bf4:bf4+2) &
         -Palpha(bf3,bf2)*temp(1:3) -Palpha(bf3+1,bf2)*temp(4:6) -Palpha(bf3+2,bf2)*temp(7:9) 
     Gbeta(bf1,bf4:bf4+2)  = Gbeta(bf1,bf4:bf4+2)  &
         - Pbeta(bf3,bf2)*temp(1:3) - Pbeta(bf3+1,bf2)*temp(4:6) - Pbeta(bf3+2,bf2)*temp(7:9)

     Galpha(bf2,bf4:bf4+2) = Galpha(bf2,bf4:bf4+2) &
         -Palpha(bf1,bf3)*temp(1:3) -Palpha(bf1,bf3+1)*temp(4:6) -Palpha(bf1,bf3+2)*temp(7:9) 
     Gbeta(bf2,bf4:bf4+2)  = Gbeta(bf2,bf4:bf4+2)  &
         - Pbeta(bf1,bf3)*temp(1:3) - Pbeta(bf1,bf3+1)*temp(4:6) - Pbeta(bf1,bf3+2)*temp(7:9)

     Galpha(bf1,bf3)   = Galpha(bf1,bf3)   -sum(Palpha(bf2,bf4:bf4+2)*temp(1:3))
     Gbeta(bf1,bf3)    = Gbeta(bf1,bf3)    -sum(Pbeta(bf2,bf4:bf4+2) *temp(1:3))
     Galpha(bf1,bf3+1) = Galpha(bf1,bf3+1) -sum(Palpha(bf2,bf4:bf4+2)*temp(4:6))
     Gbeta(bf1,bf3+1)  = Gbeta(bf1,bf3+1)  -sum(Pbeta(bf2,bf4:bf4+2) *temp(4:6))
     Galpha(bf1,bf3+2) = Galpha(bf1,bf3+2) -sum(Palpha(bf2,bf4:bf4+2)*temp(7:9))
     Gbeta(bf1,bf3+2)  = Gbeta(bf1,bf3+2)  -sum(Pbeta(bf2,bf4:bf4+2) *temp(7:9))

     Galpha(bf2,bf3)   = Galpha(bf2,bf3)   -sum(Palpha(bf1,bf4:bf4+2)*temp(1:3))
     Gbeta(bf2,bf3)    = Gbeta(bf2,bf3)    -sum(Pbeta(bf1,bf4:bf4+2) *temp(1:3))
     Galpha(bf2,bf3+1) = Galpha(bf2,bf3+1) -sum(Palpha(bf1,bf4:bf4+2)*temp(4:6))
     Gbeta(bf2,bf3+1)  = Gbeta(bf2,bf3+1)  -sum(Pbeta(bf1,bf4:bf4+2) *temp(4:6))
     Galpha(bf2,bf3+2) = Galpha(bf2,bf3+2) -sum(Palpha(bf1,bf4:bf4+2)*temp(7:9))
     Gbeta(bf2,bf3+2)  = Gbeta(bf2,bf3+2)  -sum(Pbeta(bf1,bf4:bf4+2) *temp(7:9))

              end do;           !** Ends the loop over the models with K angular momentum

        end if;  !*** End selection concentric vs non-concentric
     end if;  !*** End selection on short-range vs long-range terms
end do;  !*** Ends loop over pp shell pairs


return
end subroutine calcSRgenRIppdigestoneRperatom



!************************************
!************************************
!************************************
!************************************



!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine printint(intt,bf1,bf2,bf3,bf4,printopt)
 implicit none
 integer, intent(in) :: bf1, bf2, bf3, bf4, printopt
 real*8, intent(in)  :: intt

if (printopt .ne. 1) then; 
 return
end if

if (abs(intt) .gt. 1d-10) then;
if (bf1 == bf2 .and. bf3 == bf4) then;
print *, "(",bf1, bf2, "|", bf3, bf4, ") = ", 4*intt
else if (bf1 == bf2 .or. bf3 == bf4) then;
print *, "(",bf1, bf2, "|", bf3, bf4, ") = ", 2*intt
else;
print *, "(",bf1, bf2, "|", bf3, bf4, ") = ", intt
end if
end if

end subroutine printint

