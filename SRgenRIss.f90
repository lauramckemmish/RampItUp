




!**********************************

subroutine calcSRgenRIssdigestoneRperatom(noggss, ggss, noatoms, atoms, maxrampdegree, Tcutoff,&
     pts,a, newsortRg,lenmodelnewsort,maxrgpa, Omrnmatrix, RIssmain, &
      Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,integrand,intcnt,Gtot,absc,weights, &
     r1, r2, r3, r4, r6, r8,r10,r12, expQmr,expQpr,erfQmr,erfQpr,evenpow,oddpow,ang,  newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     :: noatoms, noggss, maxrampdegree, pts, a, nobasisfun 
  real*8, intent(in)      :: atoms(1:5,1:noatoms), Tcutoff, ggss(1:20,1:noggss)
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  integer                 :: i,j, ntot, k, atid, m, angl, angm, angk, w, warnswit(1:10), v, intcnt, printopt
  integer                 :: lsp, ksp,  lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, bf1, bf2, bf3, bf4
  real*8                  :: nc1, nc2, z, Ax, Ay, Az, Qx, Qy, Qz, QAx, QAy, QAz, QA2, QA, sqrtpi
  real*8                  :: zeta, sqrtzeta, cf, cf2, T, pi, ang(1:25), Ptotbf3bf4, inttlist(1:maxrgpa)
  real*8                  :: cpf,cpfmat(1:25), Omrnmatrix(1:2*pts+1,0:maxrampdegree)
  real*8                  :: absc(1:2*pts+1), weights(1:2*pts+1),RIssmain(0:maxrampdegree,0:4)
  real*8                  :: r1(1:2*pts+1), r2(1:2*pts+1), r3(1:2*pts+1), r4(1:2*pts+1)
  real*8                  :: r6(1:2*pts+1),  r8(1:2*pts+1),  r10(1:2*pts+1),  r12(1:2*pts+1)
  real*8                  :: expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1)
  real*8                  :: integrand(1:2*pts+1,0:4), evenpow(1:2*pts+1), oddpow(1:2*pts+1), eeint2
  real*8                  :: newsortRg(1:4,lenmodelnewsort(3, 1)+1:lenmodelnewsort(4,maxrgpa)), eeint,eeintpart(1:5), x, y
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8                  :: acum(1:4), t11, stor(1:4), tempint
  integer :: lengthtwoelec, braketno
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SRgenRIss calculates integrals of the form (SS|ss) where the |ss) is not concentric and large exponent. 
  !
  ! All integrals are correct. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "DEBUG.h"
printopt = 0
#ifdef DEBUG_PRINT
print *, "in SR (RI|ss)"
printopt = 1;
#endif

  warnswit = 0;  
  pi = dacos(-1d0);  sqrtpi = sqrt(pi);
  
   bf1 = int(lenmodelnewsort(1,1))
  !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
  !   MAIN LOOP
  !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
        Ax      = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)
  do j=1,noggss
     !*****!!!!!!!!! LOOP OVER ss SHELL PAIRS   !!!!***********!!!!!
     Qx       = ggss(4,j);        Qy       = ggss(5,j);       Qz       = ggss(6,j)
     zeta     = ggss(8,j);        
    !*****!!!!!!!!! LOOP OVER atom centers   !!!!***********!!!!!
        QAx     = Qx - Ax;          QAy     = Qy - Ay;          QAz     = Qz - Az; 
        QA2     = sum((ggss(4:6,j)-atoms(3:5,a))**2);     
        T       = zeta*QA2;
        
        if (T < Tcutoff) then
x = 0d0;
     bf3      = ggss(1,j);        bf4      = ggss(2,j);       nc2      = ggss(9,j)
 sqrtzeta = sqrt(zeta);      z        = 1d0/sqrtzeta;
     Ptotbf3bf4 = 2d0* Ptot(bf3, bf4);
  QA = sqrt(QA2);

#ifdef INTCNTPRT
intcnt = intcnt + maxrgpa
#endif
           !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
           !   SELECT SHORT-RANGE SHELL QUARTETS ONLY
           !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
           
          if (QA2 .gt. 1d-10 ) then;
              !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              !   NON-CONCENTRIC CASE
              !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              !*** Preliminaries
              
              call calcfullang(ang,QAx,QAy,QAz);

    erfQmr = derf(sqrtzeta*(QA-absc))
    erfQpr = derf(sqrtzeta*(QA+absc))
    expQmr = dexp(-zeta*(QA-absc)**2)
    expQpr = dexp(-zeta*(QA+absc)**2)

              call calcr2lrall(integrand, z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, &
                  oddpow,evenpow)

              do v=1,2*pts+1;   integrand(v,0:4) = integrand(v,0:4)*weights(v);   end do

    call dgemm('T','N',maxrampdegree+1,5,2*pts+1,1d0,Omrnmatrix,2*pts+1,integrand,2*pts+1,0d0,RIssmain,maxrampdegree+1)

                cpfmat(1) = ang(1)*157.91367041742973790d0*nc2
                cpfmat(2:4) = ang(2:4)*157.91367041742973790d0/(3d0)*nc2
                cpfmat(5:9) = ang(5:9)*157.91367041742973790d0/(5d0)*nc2
                cpfmat(10:16) = ang(10:16)*157.91367041742973790d0/(7d0)*nc2
                cpfmat(17:25) = ang(17:25)*157.91367041742973790d0/(9d0)*nc2


              !*** Main loop over relevant ramp model components

 acum = 0d0; 
stor(1:4) = (/-Palpha(bf3,bf1),-Pbeta(bf3,bf1),-Palpha(bf4,bf1),-Pbeta(bf4,bf1)/)

              do i=1, maxrgpa
                 eeintpart = 0d0;
                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)-6,5
                    eeintpart(1) = eeintpart(1) + RIssmain(int(newsortRg(2,w)),int(newsortRg(4,w)))&
                                     * cpfmat(int(newsortRg(3,w)))*newsortRg(1,w)
                    eeintpart(2) = eeintpart(2) + RIssmain(int(newsortRg(2,w+1)),int(newsortRg(4,w+1)))&
                                     * cpfmat(int(newsortRg(3,w+1)))*newsortRg(1,w+1)
                    eeintpart(3) = eeintpart(3) + RIssmain(int(newsortRg(2,w+2)),int(newsortRg(4,w+2)))&
                                     * cpfmat(int(newsortRg(3,w+2)))*newsortRg(1,w+2)
                    eeintpart(4) = eeintpart(4) + RIssmain(int(newsortRg(2,w+3)),int(newsortRg(4,w+3)))&
                                     * cpfmat(int(newsortRg(3,w+3)))*newsortRg(1,w+3)
                    eeintpart(5) = eeintpart(5) + RIssmain(int(newsortRg(2,w+4)),int(newsortRg(4,w+4)))&
                                     * cpfmat(int(newsortRg(3,w+4)))*newsortRg(1,w+4)
                 end do
!finishing the ones not divisible by 5
v = w;
eeint = sum(eeintpart)
do w= v, lenmodelnewsort(4,i)
  eeint = eeint + RIssmain(int(newsortRg(2,w)),int(newsortRg(4,w))) * cpfmat(int(newsortRg(3,w)))*newsortRg(1,w)
end do   
inttlist(i) = eeint
end do

do i=1,maxrgpa
    eeint = inttlist(i)
     bf2 = int(lenmodelnewsort(2,i))

     Gtot(bf2,bf1) = Gtot(bf2,bf1) + Ptotbf3bf4* eeint 

     x = x+ 2d0*Ptot(bf2,bf1)*eeint;
     
     acum(1)  = acum(1) -Palpha(bf2,bf3)*eeint
     acum(2)  = acum(2)  - Pbeta(bf2,bf3)*eeint
     acum(3)  = acum(3) -Palpha(bf2,bf4)*eeint
     acum(4)  = acum(4)  - Pbeta(bf2,bf4)*eeint
     Galpha(bf2,bf4) = Galpha(bf2,bf4) +stor(1)*eeint
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  +stor(2)*eeint
     Galpha(bf2,bf3) = Galpha(bf2,bf3) +stor(3)*eeint
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  +stor(4)*eeint

#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if;  

if (braketno == 1 .and. i .ne. j) then;
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,2*eeint,braketno)
else
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,eeint,braketno)
end if
end if
#endif


              end do;           !** Ends the loop over the models with K angular momentum

     Galpha(bf4,bf1) = Galpha(bf4,bf1) +acum(1) 
     Gbeta(bf4,bf1)  = Gbeta(bf4,bf1)  +acum(2)
     Galpha(bf3,bf1) = Galpha(bf3,bf1) +acum(3)
     Gbeta(bf3,bf1)  = Gbeta(bf3,bf1)  +acum(4)

          else
!*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              !   CONCENTRIC CASE
              !  (only SI|ss) is non-zero
              !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              do v=1,2*pts+1;  integrand(v,0) = weights(v)*r1(v)*erf(r1(v)*sqrtzeta);       end do

              !*** Main loop over relevant ramp model components
              cpf = 19.739208802178717238d0*z**3*nc2;
do ntot=0,maxrampdegree
RIssmain(ntot,0) = sum(integrand(1:2*pts+1,0)*Omrnmatrix(1:2*pts+1,ntot))*cpf
end do
                    bf1 = int(lenmodelnewsort(1,1))
              do i=1, maxrgpa
                 eeint = 0d0;
                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    nc1  = newsortRg(1,w)
                    ntot = newsortRg(2,w)
                    lsp  = newsortRg(4,w)
  if (lsp == 0) then; 
                    eeint = eeint + RIssmain(ntot,0)*nc1
  else 
    exit
  end if
                 end do        
                    bf2 = int(lenmodelnewsort(2,i))

     Gtot(bf2,bf1) = Gtot(bf2,bf1) + Ptotbf3bf4* eeint 

     x = x+ 2d0*Ptot(bf2,bf1)*eeint;


#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if; 
if (braketno == 1 .and. bf1 == 1 .and. bf2 == 2) then
if (i .ne. j) then
 print *, i, j, bf1,bf2,bf3,bf4, 2*eeint
else
 print *, i, j, bf1,bf2,bf3,bf4, eeint
end if
else if (braketno == 1 .and. bf1 == 2 .and. bf2 == 1) then
if (i .ne. j) then
 print *, bf1,bf2,bf3,bf4, 2*eeint
else
 print *, bf1,bf2,bf3,bf4, eeint
end if
end if

!if (braketno == 1 .and. i .ne. j) then;
!    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,2*eeint,braketno)
!else
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,eeint,2)
!end if
end if
#endif
     
     Galpha(bf1,bf4) = Galpha(bf1,bf4) -Palpha(bf2,bf3)*eeint 
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - Pbeta(bf2,bf3)*eeint
     Galpha(bf1,bf3) = Galpha(bf1,bf3) -Palpha(bf2,bf4)*eeint
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - Pbeta(bf2,bf4)*eeint
     Galpha(bf2,bf4) = Galpha(bf2,bf4) -Palpha(bf1,bf3)*eeint
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - Pbeta(bf1,bf3)*eeint
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -Palpha(bf1,bf4)*eeint
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - Pbeta(bf1,bf4)*eeint
                   
              end do;           !** Ends the loop over the models with K angular momentum
           end if             !** Ends concentric vs non-concentric decision
     Gtot(bf3,bf4) = Gtot(bf3,bf4) + x
        end if;            !** Ends the decision on whether to use SR or LR code. 
  end do;           !** Ends loop over Gaussian ss shell pairs

  return
end subroutine calcSRgenRIssdigestoneRperatom



!**********************************

subroutine calcSRgenRIssdigest(noggss, ggss, noatoms, atoms, maxrampdegree, Tcutoff,&
     pts,a, newsortRg,lenmodelnewsort,maxrgpa, Omrnmatrix, RIssmain, &
      Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,integrand, ang,intcnt,Gtot, & 
           newtwo,lengthtwoelec)
  implicit none
  integer, intent(in)     :: noatoms, noggss, maxrampdegree, pts, a, nobasisfun 
  real*8, intent(in)      :: atoms(1:5,1:noatoms), Tcutoff, ggss(1:20,1:noggss)
  real*8                  :: Gtot(1:nobasisfun,1:nobasisfun)
  integer                 :: i,j, ntot, k, atid, m, angl, angm, angk, w, warnswit(1:10), v, intcnt, printopt
  integer                 :: lsp, ksp,  lenmodelnewsort(1:4,1:maxrgpa), maxrgpa, bf1, bf2, bf3, bf4
  real*8                  :: nc1, nc2, z, Ax, Ay, Az, Qx, Qy, Qz, QAx, QAy, QAz, QA2, QA, sqrtpi
  real*8                  :: zeta, sqrtzeta, cf, cf2, T, pi, ang(1:25), Ptotbf3bf4
  real*8                  :: cpf,cpfmat(1:25), Omrnmatrix(1:2*pts+1,0:maxrampdegree)
  real*8                  :: absc(1:2*pts+1), weights(1:2*pts+1),RIssmain(0:maxrampdegree,0:4)
  real*8                  :: r1(1:2*pts+1), r2(1:2*pts+1), r3(1:2*pts+1), r4(1:2*pts+1),  tempint
  real*8                  :: r6(1:2*pts+1),  r8(1:2*pts+1),  r10(1:2*pts+1),  r12(1:2*pts+1)
  real*8                  :: expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1)
  real*8                  :: integrand(1:2*pts+1,0:4), evenpow(1:2*pts+1), oddpow(1:2*pts+1)
  real*8                  :: newsortRg(1:4,lenmodelnewsort(3, 1)+1:lenmodelnewsort(4,maxrgpa)), eeint, x, y
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  integer :: lengthtwoelec, braketno
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1) 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! SRgenRIss calculates integrals of the form (SS|ss) where the |ss) is not concentric and large exponent. 
  !
  ! All integrals are correct. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "DEBUG.h"
printopt = 0
#ifdef DEBUG_PRINT
print *, "in SR (RI|ss)"
printopt = 1;
#endif
  warnswit = 0;  
  pi = dacos(-1d0);  sqrtpi = sqrt(pi);

  call getquad(weights,absc,pts);  call evalallrk(r1,r2,r3,r4, r6, r8,absc,pts);
  r10 = r8*r2;  r12 = r8*r4;  
  
  do ntot = 0,maxrampdegree;    call evalOmrn(Omrnmatrix(1:2*pts+1,ntot),absc,pts,ntot);  end do

  !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
  !   MAIN LOOP
  !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
  do j=1,noggss
     !*****!!!!!!!!! LOOP OVER ss SHELL PAIRS   !!!!***********!!!!!
     bf3      = ggss(1,j);        bf4      = ggss(2,j);       nc2      = ggss(9,j)
     Qx       = ggss(4,j);        Qy       = ggss(5,j);       Qz       = ggss(6,j)
     zeta     = ggss(8,j);        sqrtzeta = sqrt(zeta);      z        = 1d0/sqrtzeta;
     Ptotbf3bf4 = 2d0* Ptot(bf3, bf4)
     !*****!!!!!!!!! LOOP OVER atom centers   !!!!***********!!!!!
        Ax      = atoms(3,a);    Ay      = atoms(4,a);    Az      = atoms(5,a)
        QAx     = Qx - Ax;          QAy     = Qy - Ay;          QAz     = Qz - Az; 
        QA2     = QAx**2+QAy**2+QAz**2;       QA = sqrt(QA2);  
        T       = zeta*QA2;
        
        if (T < Tcutoff) then
x = 0d0;

intcnt = intcnt + maxrgpa
           !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
           !   SELECT SHORT-RANGE SHELL QUARTETS ONLY
           !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
           
          if (QA2 .gt. 1d-10 ) then;
              !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              !   NON-CONCENTRIC CASE
              !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              !*** Preliminaries
              
              call calcfullang(ang,QAx,QAy,QAz);

              call allerfexp(erfQmr,erfQpr,expQmr,expQpr,absc,pts,QA, sqrtzeta,zeta); 

              call calcr2lrall(integrand(1:2*pts+1,0:4), z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, &
                  oddpow,evenpow)

              do v=1,2*pts+1;   integrand(v,0:4) = integrand(v,0:4)*weights(v);   end do

    call dgemm('T','N',maxrampdegree+1,5,2*pts+1,1d0,Omrnmatrix,2*pts+1,integrand,2*pts+1,0d0,RIssmain,maxrampdegree+1)

              angk = 0;
              do angL = 0,4; do angm=-angL,angL; angk = angk + 1;
                cpfmat(angk) = 157.91367041742973790d0*ang(angk)/(2d0*angL+1d0)*nc2;
              end do; end do;


              !*** Main loop over relevant ramp model components
              do i=1, maxrgpa
                 eeint = 0d0;
                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    nc1  = newsortRg(1,w)
                    ntot = newsortRg(2,w)
                    ksp  = newsortRg(3,w)
                    lsp  = newsortRg(4,w)
                    eeint = eeint + RIssmain(ntot,lsp)* cpfmat(ksp)*nc1
                 end do         

                    bf1 = int(lenmodelnewsort(1,i))
                    bf2 = int(lenmodelnewsort(2,i))

#ifdef DEBUG_PRINT
printopt = 1;
 call printint(eeint,bf1,bf2,bf3,bf4, printopt)
#endif

     Gtot(bf1,bf2) = Gtot(bf1,bf2) + Ptotbf3bf4* eeint 

     x = x+ 2d0*Ptot(bf1,bf2)*eeint;
     
     Galpha(bf1,bf4) = Galpha(bf1,bf4) -Palpha(bf2,bf3)*eeint 
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - Pbeta(bf2,bf3)*eeint
     Galpha(bf1,bf3) = Galpha(bf1,bf3) -Palpha(bf2,bf4)*eeint
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - Pbeta(bf2,bf4)*eeint
     Galpha(bf2,bf4) = Galpha(bf2,bf4) -Palpha(bf1,bf3)*eeint
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - Pbeta(bf1,bf3)*eeint
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -Palpha(bf1,bf4)*eeint
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - Pbeta(bf1,bf4)*eeint


#ifdef INTPRINT
if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if; 
if (braketno == 1 .and. i .ne. j) then;
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,2*eeint,braketno)
else
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,eeint,braketno)
end if
end if
#endif
              end do;           !** Ends the loop over the models with K angular momentum
          else
!*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              !   CONCENTRIC CASE
              !  (only SI|ss) is non-zero
              !*****!!!!!!!!!!!!!***********!!!!!!!!!!!!!***********!!!!!!!!!!***********!!!!!!!!!
              do v=1,2*pts+1;  integrand(v,0) = weights(v)*r1(v)*erf(r1(v)*sqrtzeta);       end do

              !*** Main loop over relevant ramp model components
              cpf = 19.739208802178717238d0*z**3;
              do i=1, maxrgpa
                 eeint = 0d0;
                 do w= lenmodelnewsort(3,i)+1, lenmodelnewsort(4,i)
                    nc1  = newsortRg(1,w)
                    ntot = newsortRg(2,w)
                    lsp  = newsortRg(4,w)
  if (lsp == 0) then; 
                    eeint = eeint + sum(integrand(1:2*pts+1,0)*Omrnmatrix(1:2*pts+1,ntot))*cpf*nc1*nc2
  end if
                 end do        
                    bf1 = int(lenmodelnewsort(1,i))
                    bf2 = int(lenmodelnewsort(2,i))

#ifdef DEBUG_PRINT
printopt = 1;
 call printint(eeint,bf1,bf2,bf3,bf4, printopt)
#endif

#ifdef INTPRINT
if (abs(Palpha(1,1)) .lt. 1d-10) then
        if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
        else; braketno = 2;
        end if; 
if (braketno == 1 .and. i .ne. j) then;
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,2*eeint,braketno)
else
    call makenewtwo(newtwo,nobasisfun,bf1,bf2,bf3,bf4,eeint,braketno)
end if
end if
#endif

     Gtot(bf2,bf1) = Gtot(bf2,bf1) + Ptotbf3bf4* eeint 

     x = x+ 2d0*Ptot(bf2,bf1)*eeint;
     
     Galpha(bf1,bf4) = Galpha(bf1,bf4) -Palpha(bf2,bf3)*eeint 
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - Pbeta(bf2,bf3)*eeint
     Galpha(bf1,bf3) = Galpha(bf1,bf3) -Palpha(bf2,bf4)*eeint
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - Pbeta(bf2,bf4)*eeint
     Galpha(bf2,bf4) = Galpha(bf2,bf4) -Palpha(bf1,bf3)*eeint
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - Pbeta(bf1,bf3)*eeint
     Galpha(bf2,bf3) = Galpha(bf2,bf3) -Palpha(bf1,bf4)*eeint
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - Pbeta(bf1,bf4)*eeint
                   

              end do;           !** Ends the loop over the models with K angular momentum

           end if             !** Ends concentric vs non-concentric decision
     Gtot(bf3,bf4) = Gtot(bf3,bf4) + x
        end if;            !** Ends the decision on whether to use SR or LR code. 
  end do;           !** Ends loop over Gaussian ss shell pairs
  


  return
end subroutine calcSRgenRIssdigest
