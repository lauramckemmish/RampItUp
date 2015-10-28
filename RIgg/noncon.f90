subroutine calcr2lrall(integrand, z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, &
         oddpow, evenpow)
  implicit none
  integer, intent(in)   :: pts
  real*8, intent(in)    :: z,  QA, r1(1:2*pts+1), r2(1:2*pts+1),r4(1:2*pts+1),r6(1:2*pts+1)
  real*8, intent(in)    :: r8(1:2*pts+1),r10(1:2*pts+1),r12(1:2*pts+1), r3(1:2*pts+1)
  real*8, intent(in)    :: expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1)
  real*8, intent(out)   :: integrand(1:2*pts+1,0:4)
  real*8                :: pfa, pfb, pfc, pf, pf1, pf3, pf2, pf4, pf5, pf6, pf7, pf8,pf9,pf10,pf11
  real*8                :: oddpow(1:2*pts+1), evenpow(1:2*pts+1), powQA(-6:6), rQA, pfd, pfe, sqrtpi
  real*8                :: z2, z4, z6, z8,z10,z12,z14
  integer               :: l
  
!!!!!!!!!!!!!!!!!* r1l0!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sqrtpi = 1.77245385090551602729816748334d0! 
rQA = 1d0/QA;
powQA(-1) = rQA;
do l=2,6
powQA(-l) = powQA(-l+1)*rQA;
end do
powQA(2) = QA**2;
do l=3,6
powQA(l) = QA*powQA(l-1);
end do

z2 = z**2;
z4 = z2**2
z6 = z2*z4
z8 = z4**2
z10 = z4*z6;
z12 = z6**2
z14 = z6*z8

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


  !*** r*(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1))*QA**3;
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


  
  !*** r*(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1))*QA**3;
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

  return
end subroutine calcr2lrall


!**********************************************************************
