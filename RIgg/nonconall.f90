subroutine calcssspppall(integrand,integrandsp,integrandpp, z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, &
         oddpow, evenpow)
  implicit none
  integer, intent(in)   :: pts
  real*8, intent(in)    :: z,  QA, r1(1:2*pts+1), r2(1:2*pts+1),r4(1:2*pts+1),r6(1:2*pts+1)
  real*8, intent(in)    :: r8(1:2*pts+1),r10(1:2*pts+1),r12(1:2*pts+1), r3(1:2*pts+1)
  real*8, intent(in)    :: expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1)
  real*8, intent(out)   :: integrand(1:2*pts+1,0:4),integrandsp(1:2*pts+1,0:4),integrandpp(1:2*pts+1,0:4)
  real*8                :: pfa, pfb, pfc, pf, pf1, pf3, pf2, pf4, pf5, pf6, pf7, pf8,pf9,pf10,pf11
  real*8                :: oddpow(1:2*pts+1), evenpow(1:2*pts+1), powQA(-7:6), rQA, pfd, pfe, sqrtpi
  integer               :: l
  real*8                :: z2, z4, z6, z8,z10,z12,z14

rQA = 1d0/QA;
z2 = z**2;
z4 = z2**2
z6 = z2*z4
z8 = z4**2
z10 = z4*z6;
z12 = z6**2
z14 = z6*z8

sqrtpi = 1.77245385090551602729816748334d0! 
rQA = 1d0/QA;
powQA(-1) = rQA;
do l=2,7
powQA(-l) = powQA(-l+1)*rQA;
end do
powQA(2) = QA**2;
do l=3,6
powQA(l) = QA*powQA(l-1);
end do


!*** r*(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
  pfa   = -sqrtpi*z**3*0.125d0;
!*** r*(exp(Qm1) - exp(Qp1))
  pfb   = -1d0*z4*rQA*0.125d0;
  !*** r2*(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pfc   = sqrtpi*z**3*0.125d0*rQA;
  integrand(1:2*pts+1,0) = r1*(pfa*(erfQmr-erfQpr) + pfb*(expQmr-expQpr))+ pfc*r2*(erfQmr + erfQpr);


!!!!!!!!!!!!!!!!!* r1l1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
  pfa   = -sqrtpi*z**3*0.125d0*QA;
  
  !*** r*(exp(Qm1) - exp(Qp1))
  pfb   = 0.125d0*z4*(0.5d0*z2*powQA(-2) - 1);

  !*** r^2(exp(Qm1) + exp(Qp1))
  pfc   = -0.125d0*z4*rQA;

! r2 gt r1,  l = 1

  !*** r^4 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pfd   = sqrtpi*z**3*0.125d0*powQA(-2)

  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pfe   = -z4*0.125d0*powQA(-2);
  
  integrand(1:2*pts+1,1) =  pfa*r1*(erfQmr-erfQpr)+ (pfb*r1+pfe*r3)*(expQmr-expQpr)+ pfc*r2*(expQmr+expQpr) &
                + r4*pfd*(erfQmr+erfQpr)


! r^2

  !*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
  pfa   = -sqrtpi*z**3*0.125d0*powQA(2);
  !*** r^(exp(Qm1) - exp(Qp1))
  pf1   = -3d0*z8*powQA(-3)*0.03125d0+z6*rQA*0.0625d0 -z4*0.125d0*QA;
  pf3   = -0.125d0*z4*rQA;
  oddpow = r1*(pf1 + pf3*r2);
  
  !*** r^2(exp(Qm1) + exp(Qp1))
  pf   = 3d0*z6*powQA(-2)*0.0625d0 -0.125d0*z4;
  integrand(1:2*pts+1,2) = pfa*r1*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr)+ pf*(expQmr+expQpr)*r2
  

! r2 gt r1,  l = 2


  !*** r^6 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pfa   = sqrtpi*z**3*0.125d0*powQA(-3)

  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pf3   = z6*0.0625d0*powQA(-3);
  pf5   = -z4*0.125d0*powQA(-3);
  oddpow = r3*(pf3+r2*pf5)
  
  !*** r^4 *(exp(Qm1) + exp(Qp1)
  pf   = pf5*QA! -z4*0.125d0*powQA(-2);
  integrand(1:2*pts+1,2) = integrand(1:2*pts+1,2) + r6*pfa*(erfQmr+erfQpr)+ r4*pf*(expQmr+expQpr) + oddpow*(expQmr-expQpr)



! r^3


  !*** r*(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1))*powQA(3);
  pf   = -sqrtpi*z**3*0.125d0*powQA(3);
  
  !*** r^odd*(exp(Qm1) - exp(Qp1))
  pf1   = 15d0*z10*0.015625d0*powQA(-4)-3d0*z8*0.03125d0*powQA(-2)+ z6*0.0625d0 -z4*0.125d0*powQA(2);
  pf3   = 3*z6*powQA(-2)*0.125d0 -0.125d0*z4;
  oddpow = r1*(pf1+r2*pf3)
 
  !*** r^even(exp(Qm1) + exp(Qp1))
  pf2   = -15d0*z8*powQA(-3)*0.03125d0 +3d0*z6*0.0625d0*rQA -0.125d0*z4*QA  ;
  pf4   =  -z4*rQA*0.125d0;
  evenpow = r2*(pf2+r2*pf4);
  integrand(1:2*pts+1,3) = r1*pf*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)



! r2 gt r1,  l = 3

  
  !*** r^8(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pf   = sqrtpi*z**3*0.125d0*powQA(-4);
  
  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pf3   = -3*z8*0.03125d0*powQA(-4);
  pf5   = z6*0.0625d0*powQA(-4) - z4*0.125d0*powQA(-2);
  pf7   = -z4*0.125d0*powQA(-4);
  oddpow = r3*(pf3 + r2*(pf5+r2*pf7))
  
  !*** r^4 *(exp(Qm1) + exp(Qp1)
  pf4   = 3*z6*0.0625d0*powQA(-3);
  pf6   = pf7*QA!  -z4*0.125d0*powQA(-3);
  evenpow = r4*(pf4+r2*pf6)
  integrand(1:2*pts+1,3) = integrand(1:2*pts+1,3)+ r8*pf*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)




! r^4


  
  !*** r*(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1))*powQA(3);
  pf   = -sqrtpi*z**3*0.125d0*QA**4;
  
  !*** r^odd*(exp(Qm1) - exp(Qp1))
  pf1   = -(z4*(105*z8+QA**2*(-30*z6+QA**2*(12*z4+QA**2*(16*QA**2-8*z2)))))*powQA(-5)*0.0078125d0;
  pf3   = -(z4*(45*z4+QA**2*(4*QA**2-12*z2)))*powQA(-3)*0.03125d0;
  pf5   =  -z4*0.125d0*rQA;
  oddpow = r1*(pf1+r2*(pf3+r2*pf5));
  
  !*** r^even(exp(Qm1) + exp(Qp1))
  pf2   = (z4*(105*z6+QA**2*(-30*z4+QA**2*(-8*QA**2+12*z2))))*powQA(-4)*0.015625d0;
  pf4   =  -(z4*(QA**2-5*z2))*powQA(-2)*0.125d0;
  evenpow = r2*(pf2+r2*pf4)
  integrand(1:2*pts+1,4) = r1*pf*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr) + evenpow*(expQmr+expQpr)  

! r2 gt r1,  l = 4
  
  !*** r^10(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pf   = sqrtpi*z**3*0.125d0*powQA(-5);
  
  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pf3   = 15*z10*0.015625d0*powQA(-5);
  pf5   = -3*z8*0.03125d0*powQA(-5) + 3*z6*0.125d0*powQA(-3);
  pf7   = z6*0.0625d0*powQA(-5) - z4*0.125d0*powQA(-3);
  pf9   = -z4*0.125d0*powQA(-5);
  oddpow = r3*(pf3+r2*(pf5+r2*(pf7+r2*pf9)))
  
  !*** r^4 *(exp(Qm1) + exp(Qp1)
  pf4   = -15*z8*0.03125d0*powQA(-4);
  pf6   = 3*z6*0.0625d0*powQA(-4) - z4*0.125d0*powQA(-2);
  pf8   = pf9*QA;! -z4*0.125d0*powQA(-4);
  evenpow = r4*(pf4 + r2*(pf6+r2*pf8))
  integrand(1:2*pts+1,4) = integrand(1:2*pts+1,4) + r10*pf*(erfQmr+erfQpr)+ oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)



!*** r*(exp(Qm1) - exp(Qp1))
! pfa   = 0.125d0*z4*powQA(-2);

!*** r^2*(exp(Qm1) + exp(Qp1))
! pfb   = -z2*0.25d0*rQA;
 integrandsp(1:2*pts+1,0) = r1*0.125d0*z4*powQA(-2)*(expQmr-expQpr)  -z2*0.25d0*rQA*r2*(expQmr+expQpr)



!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pfa   = -1.77245385090551602729816748334d0*z**3*0.125d0;

!*** r*(exp(Qm1) - exp(Qp1))
 pf1   = -z6*0.125d0*powQA(-3) - z4*0.125d0*rQA;
 pf3   = -0.25d0*z2*rQA;
 oddpow = r1*(pf1+r2*pf3)

!*** r^2(exp(Qm1) + exp(Qp1))
 pfb   = 0.25d0*z4*powQA(-2);
 integrandsp(1:2*pts+1,1) = r1*pfa*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr) + r2*pfb*(expQmr+expQpr)


!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pf   = -1.77245385090551602729816748334d0*z**3*QA*0.25d0;

!*** r^odd*(exp(Qm1) - exp(Qp1))
 pf1   = 9*z8*0.03125d0*powQA(-4) + z6*0.125d0*powQA(-2) - z4*0.25d0;
 pf3   = z4*0.5d0*powQA(-2);
 oddpow = r1*(pf1+r2*pf3) 

!*** r^even*(exp(Qm1) + exp(Qp1))
 pf2   = -9*z6*0.0625d0*powQA(-3) - z4*0.25d0*rQA;
 pf4   = -z2*0.25d0*rQA;
 evenpow = r2*(pf2+r2*pf4)
 integrandsp(1:2*pts+1,2) =  r1*pf*(erfQmr-erfQpr)+ oddpow*(expQmr-expQpr) + evenpow*(expQmr+expQpr)


!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pf   = -3*1.77245385090551602729816748334d0*z**3*0.125d0*QA**2;

!*** r^(exp(Qm1) - exp(Qp1))
 pf1   = -15*z10*0.0625d0*powQA(-5) - 9*z8*0.03125d0*powQA(-3) + 3*z6*0.0625d0*rQA - 3*QA*z4*0.125d0;
 pf3   = -27*z6*0.0625d0*powQA(-3) - 3*z4*0.125d0*rQA;
 pf5   = -z2*0.25d0*rQA;
 oddpow = r1*(pf1+r2*(pf3+r2*pf5));

!*** r^2(exp(Qm1) + exp(Qp1))
 pf2   = 15*z8*0.125d0*powQA(-4) + 9*z6*0.0625d0*powQA(-2) - 3*z4*0.125d0;
 pf4   = 7*z4*0.125d0*powQA(-2);
 evenpow = r2*(pf2+r2*pf4)
 integrandsp(1:2*pts+1,3) = r1*pf*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)


!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pf   = -1.77245385090551602729816748334d0*z**3*powQA(3)*0.5d0;

!*** r^(exp(Qm1) - exp(Qp1))
 pf1   = 525*z12*0.0078125d0*powQA(-6) + 15*z10*0.0625d0*powQA(-4) -3*z8*0.125d0*powQA(-2) + z6*0.25d0 &
       - QA**2*z4*0.5d0;
 pf3   = 15*z8*0.5d0*powQA(-4) + 3*z6*0.5d0*powQA(-2) - z4*0.5d0;
 pf5   = 11*z4*0.125d0*powQA(-2) ;
 oddpow = r1*(pf1+r2*(pf3+r2*pf5));

!*** r^2(exp(Qm1) + exp(Qp1))
 pf2   = -525*z10*0.015625d0*powQA(-5) - 15*z8*0.125d0*powQA(-3) + z6*0.75d0*rQA - QA*z4*0.5d0;
 pf4   = -65*z6*0.0625d0*powQA(-3) - z4*0.5d0*rQA;
 pf6   = -z2*0.25d0*rQA ;
 evenpow = r2*(pf2+r2*(pf4+r2*pf6));
 integrandsp(1:2*pts+1,4) = r1*pf*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)

!*** r^2(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pfa   = -1.77245385090551602729816748334d0*z**3*0.125d0*powQA(-2);

!*** r^2(exp(zeta*Qm1^2) + exp(zeta*Qp1^2));
 pfb   = z2*0.25d0*rQA;
 integrandsp(1:2*pts+1,0) = integrandsp(1:2*pts+1,0)+  r2*(pfa*(erfQmr+erfQpr) + pfb*(expQmr+expQpr))


!*** r^4 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pfa   = -1.77245385090551602729816748334d0*z**3*0.25d0*powQA(-3);

!*** r^3 *(exp(Qm1) - exp(Qp1)
 pfb   = z4*0.25d0*powQA(-3)+ z2*0.25d0*rQA;
 integrandsp(1:2*pts+1,1) = integrandsp(1:2*pts+1,1)+ pfa*r4*(erfQmr+erfQpr) + pfb*r3*(expQmr-expQpr);



!*** r^6 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pfa   = -3*1.77245385090551602729816748334d0*z**3*0.125d0*powQA(-4);

!*** r^3 *(exp(Qm1) - exp(Qp1)
 pf3   = -3*z6*0.0625d0*powQA(-4) - z4*0.125d0*powQA(-2);
 pf5   = 3*z4*0.125d0*powQA(-4);
 oddpow = r3*(pf3+r2*pf5)

!*** r^4 *(exp(Qm1) + exp(Qp1)
 pfb   = 3*z4*0.125d0*powQA(-3) + z2*0.25d0*rQA ;
 integrandsp(1:2*pts+1,2) = integrandsp(1:2*pts+1,2) + pfa*r6*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr)+ r4*pfb*(expQmr+expQpr)


!*** r^8(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pf   = -1.77245385090551602729816748334d0*z**3*0.5d0*powQA(-5);

!*** r^3 *(exp(Qm1) - exp(Qp1)
 pf3   = 3*z8*0.125d0*powQA(-5)+3*z6*0.0625d0*powQA(-3) ;
 pf5   = -z6*0.25d0*powQA(-5) + z4*0.5d0*powQA(-3) + z2*0.25d0*rQA;
 pf7   = z4*0.5d0*powQA(-5);
 oddpow = r3*(pf3+r2*(pf5+r2*pf7));

!*** r^4 *(exp(Qm1) + exp(Qp1)
 pf4   = - 3*z6*0.25d0*powQA(-4) - 3*z4*0.125d0*powQA(-2);
 pf6   = pf7*QA!  z4*0.5d0*powQA(-4);
 evenpow = r4*(pf4+r2*pf6)
 integrandsp(1:2*pts+1,3) = integrandsp(1:2*pts+1,3) + pf*r8*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)

!*** r^10(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pf   = -5*1.77245385090551602729816748334d0*z**3*0.125d0*powQA(-6);

!*** r^pow *(exp(Qm1) - exp(Qp1)
 pf3   = -75*z10*0.015625d0*powQA(-6) - 15*z8*0.03125d0*powQA(-4) ;
 pf5   = 15*z8*0.03125d0*powQA(-6) - 15*z6*0.125d0*powQA(-4) - 3*z4*0.25d0*powQA(-2);
 pf7   = -5*z6*0.0625d0*powQA(-6) + 5*z4*0.125d0*powQA(-4);
 pf9   = 5*z4*0.125d0*powQA(-6)
 oddpow = r3*(pf3+r2*(pf5+r2*(pf7+r2*pf9))); 

!*** r^even *(exp(Qm1) + exp(Qp1)
 pf4   = 75*z8*0.03125d0*powQA(-5) + 15*z6*0.0625d0*powQA(-3) ;
 pf6   = -15*z6*0.0625d0*powQA(-5) + 5*z4*0.125d0*powQA(-3) + z2*0.25d0*rQA;
 pf8   = pf9*QA ! 5*z4*0.125d0*powQA(-5);
 evenpow = r4*(pf4+r2*(pf6+r2*pf8));
 integrandsp(1:2*pts+1,4) = integrandsp(1:2*pts+1,4) + pf*r10*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr) + evenpow*(expQmr+expQpr)


  !*** r*(exp(Qm1) - exp(Qp1))
  pf1  = -z4*0.25d0*powQA(-3) - z2*0.25d0*rQA;
  pf3   = -1d0*0.5d0*rQA;
  oddpow = r1*(pf1+r2*pf3);
  
  !*** r^2*(exp(Qm1) + exp(Qp1))
  pf   = 0.5d0 + z2*0.5d0*powQA(-2);
  integrandpp(1:2*pts+1,0) = oddpow*(expQmr-expQpr) + pf*r2*(expQmr+expQpr)


  !*** r*(exp(Qm1) - exp(Qp1))
  pf1   = 3*z6*0.125d0*powQA(-4) + 3*z4*0.125d0*powQA(-2)
  pf3   = 0.5d0 + 3*z2*0.25d0*powQA(-2);
  oddpow = r1*(pf1 + r2*pf3);
  
  !*** r^2(exp(Qm1) + exp(Qp1))
  pf2   = -3*z4*0.25d0*powQA(-3) - 3*z2*0.25d0*rQA;
  pf4   = -1d0*0.5d0*rQA;
  evenpow = r2*(pf2+r2*pf4)
  integrandpp(1:2*pts+1,1) =  oddpow*(expQmr-expQpr) + evenpow*(expQmr+expQpr)


  !*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
  pf   = -1.77245385090551602729816748334d0*z**3*0.25d0;
  
  !*** r^(exp(Qm1) - exp(Qp1))
  pf1   = -9*z8*0.125d0*powQA(-5) -13*z6*0.0625d0*powQA(-3) - z4*0.25d0*rQA;
  pf3   = -17*z4*0.125d0*powQA(-3) - 3*z2*0.5d0*rQA;
  pf5   =-1*0.5d0*rQA;
  oddpow = r1*(pf1+r2*(pf3+r2*pf5));
  
  !*** r^2(exp(Qm1) + exp(Qp1))
  pf2   = 9*z6*0.25d0*powQA(-4)+13*z4*0.125d0*powQA(-2);
  pf4   =0.5d0 + 5*z2*0.25d0*powQA(-2);
  evenpow = r2*(pf2+r2*pf4)
  integrandpp(1:2*pts+1,2) = + r1*pf*(erfQmr-erfQpr) + evenpow*(expQmr+expQpr)+ oddpow*(expQmr-expQpr)


  !*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
  pf   = -3*1.77245385090551602729816748334d0*z**3*QA*0.25d0

  !*** r^(exp(Qm1) - exp(Qp1))
  pf1   = 75*z10*0.0625d0*powQA(-6) +87*z8*0.03125d0*powQA(-4) + 3*z6*0.125d0*powQA(-2) - 3*z4*0.25d0;
  pf3   = 141*z6*0.0625d0*powQA(-4) +39*z4*0.125d0*powQA(-2);
  pf5   = 0.5d0+2*z2*powQA(-2);
  oddpow = r1*(pf1+r2*(pf3+r2*pf5));

  
  !*** r^2(exp(Qm1) + exp(Qp1))
  pf2   = -75*z8*0.125d0*powQA(-5) - 87*z6*0.0625d0*powQA(-3) - 3*z4*0.25d0*rQA;
  pf4   = -41*z4*0.125d0*powQA(-3) - 5*z2*0.5d0*rQA;
  pf6   = -1*0.5d0*rQA;
  evenpow = r2*(pf2+r2*(pf4+r2*pf6));
  integrandpp(1:2*pts+1,3) = + r1*pf*(erfQmr-erfQpr) + evenpow*(expQmr+expQpr)+ oddpow*(expQmr-expQpr)

  !*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
  pf   = -3*1.77245385090551602729816748334d0*z**3*QA**2*0.5d0

  !*** r^(exp(Qm1) - exp(Qp1))
  pf1   = -1575*z12*0.015625d0*powQA(-7) -765*z10*0.015625d0*powQA(-5)  -9*z8*0.125d0*powQA(-3) + 3*z6*0.25d0*rQA&
             -3*QA*z4*0.5d0;
  pf3   = -1485*z8*0.03125d0*powQA(-5) -87*z6*0.25d0*powQA(-3) - 3*z4*0.5d0*rQA;
  pf5   = -87*z4*0.125d0*powQA(-3) - 15*z2*0.25d0*rQA ;
  pf7   = -1*0.5d0*rQA;
  oddpow = r1*(pf1+r2*(pf3+r2*(pf5+r2*pf7)));
  
  !*** r^2(exp(Qm1) + exp(Qp1))
  pf2   = 1575*z10*0.03125d0*powQA(-6) +765*z8*0.03125d0*powQA(-4) + 9*z6*0.25d0*powQA(-2) - 3*z4*0.5d0;
  pf4   = 435*z6*0.0625d0*powQA(-4) +93*z4*0.125d0*powQA(-2);
  pf6   = 0.5d0 + 3*z2*powQA(-2);
  evenpow = r2*(pf2+r2*(pf4+r2*pf6));
  integrandpp(1:2*pts+1,4) = + r1*pf*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)



  
  !*** r^2(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pfa   = 1.77245385090551602729816748334d0*z**3*0.25d0*powQA(-3);
  
  !*** r^2(exp(zeta*Qm1^2) + exp(zeta*Qp1^2));
  pfb   = -0.5d0 -z2*0.5d0*powQA(-2);
  
  !*** r^3(exp(zeta*Qm1^2) - exp(zeta*Qp1^2));
  pfc   = 0.5d0*rQA;
  integrandpp(1:2*pts+1,0) =  integrandpp(1:2*pts+1,0) + r2*(pfa*(erfQmr+erfQpr) +pfb*(expQmr+expQpr))+ pfc*r3*(expQmr-expQpr)


  !*** r^4 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pfa   = 3*1.77245385090551602729816748334d0*z**3*0.25d0*powQA(-4);
  
  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pfb   = -0.5d0 - 3*z4*0.25d0*powQA(-4) - 3*z2*0.25d0*powQA(-2);
  
  !*** r^4 *(exp(Qm1) + exp(Qp1)
!  pfc   = 1d0*0.5d0*rQA;
  integrandpp(1:2*pts+1,1) =  integrandpp(1:2*pts+1,1) + pfa*r4*(erfQmr+erfQpr)+ pfb*r3*(expQmr-expQpr) + r4*pfc*(expQmr+expQpr)

 

  !*** r^6 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pfa   = 3*1.77245385090551602729816748334d0*z**3*0.5d0*powQA(-5);
  
  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pf3   = 3*z6*0.25d0*powQA(-5) + 5*z4*0.125d0*powQA(-3)+z2*0.25d0*rQA;
  pf5   = 1*0.5d0*rQA - 3*z4*0.5d0*powQA(-5);
  oddpow = r3*(pf3+r2*pf5)
  
  !*** r^4 *(exp(Qm1) + exp(Qp1)
  pfb   = -0.5d0 - 3*z4*0.5d0*powQA(-4) - 5*z2*0.25d0*powQA(-2) ;
  integrandpp(1:2*pts+1,2) =  integrandpp(1:2*pts+1,2) + pfa*r6*(erfQmr+erfQpr)  + oddpow*(expQmr-expQpr) + pfb*r4*(expQmr+expQpr)



  !*** r^8(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pf   = 5*1.77245385090551602729816748334d0*z**3*0.5d0*powQA(-6);

  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pf3   = -15*z8*0.125d0*powQA(-6)-21*z6*0.0625d0*powQA(-4) - 3*z4*0.125d0*powQA(-2);
  pf5   = - 0.5d0 + 5*z6*0.25d0*powQA(-6) - 5*z4*0.5d0*powQA(-4) - 2*z2*powQA(-2);
  pf7   = -5*z4*0.5d0*powQA(-6);
  oddpow = r3*(pf3+r2*(pf5+r2*pf7));

  !*** r^4 *(exp(Qm1) + exp(Qp1)
  pf4   = 15*z6*0.25d0*powQA(-5) + 21*z4*0.125d0*powQA(-3) + 3*z2*0.25d0*rQA;
  pf6   = 1*0.5d0*rQA - 5*z4*0.5d0*powQA(-5);
  evenpow = r4*(pf4+r2*pf6);
  integrandpp(1:2*pts+1,3) =  integrandpp(1:2*pts+1,3) +  pf*r8*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)


  !*** r^10(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
  pf   = 15*1.77245385090551602729816748334d0*z**3*0.25d0*powQA(-7);
  
  !*** r^3 *(exp(Qm1) - exp(Qp1)
  pf3   = 225*z10*0.03125d0*powQA(-7) + 135*z8*0.03125d0*powQA(-5)  + 15*z6*0.0625d0*powQA(-3);
  pf5   =  -45*z8*0.0625d0*powQA(-7)  + 45*z6*0.25d0*powQA(-5) + 57*z4*0.125d0*powQA(-3)+ 3*z2*0.5d0*rQA;
  pf7   = 1*0.5d0*rQA + 15*z6*0.125d0*powQA(-7) - 15*z4*0.25d0*powQA(-5);
  pf9   = -15*z4*0.25d0*powQA(-7);
  oddpow = r3*(pf3+r2*(pf5+r2*(pf7+r2*pf9)))
   
  !*** r^4 *(exp(Qm1) + exp(Qp1)
  pf4   = -225*z8*0.0625d0*powQA(-6) -135*z6*0.0625d0*powQA(-4) -15*z4*0.125d0*powQA(-2);
  pf6   = 45*z6*0.125d0*powQA(-6) -15*z4*0.25d0*powQA(-4) - 3*z2*powQA(-2) - 0.5d0 ;
  pf8   = -15*z4*0.25d0*powQA(-6);
  evenpow = r4*(pf4+r2*(pf6+r2*pf8));
  integrandpp(1:2*pts+1,4) =  integrandpp(1:2*pts+1,4) +pf*r10*(erfQmr+erfQpr)+ oddpow*(expQmr-expQpr) + evenpow*(expQmr+expQpr)

return
end subroutine calcssspppall


!***********************************************************