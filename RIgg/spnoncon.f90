subroutine calcspall(integrandsp, z, QA,erfQmr,erfQpr,expQmr,expQpr,r1,r2,r3,r4,r6,r8,r10,r12,pts, &
         oddpow, evenpow)
  implicit none
  integer, intent(in)   :: pts
  real*8, intent(in)    :: z,  QA, r1(1:2*pts+1), r2(1:2*pts+1),r4(1:2*pts+1),r6(1:2*pts+1)
  real*8, intent(in)    :: r8(1:2*pts+1),r10(1:2*pts+1),r12(1:2*pts+1), r3(1:2*pts+1)
  real*8, intent(in)    :: expQmr(1:2*pts+1), expQpr(1:2*pts+1),erfQmr(1:2*pts+1), erfQpr(1:2*pts+1)
  real*8, intent(out)   :: integrandsp(1:2*pts+1,0:4)
  real*8                :: pfa, pfb, pfc, pf, pf1, pf3, pf2, pf4, pf5, pf6, pf7, pf8,pf9,pf10,pf11
  real*8                :: oddpow(1:2*pts+1), evenpow(1:2*pts+1), powQA(-6:6), rQA, pfd, pfe, sqrtpi
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

!*** r*(exp(Qm1) - exp(Qp1))
 pfa   = 0.125d0*z4*rQA**2;

!*** r^2*(exp(Qm1) + exp(Qp1))
 pfb   = -z2*0.25d0*rQA;
 integrandsp(1:2*pts+1,0) = r1*pfa*(expQmr-expQpr) + pfb*r2*(expQmr+expQpr)



!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pfa   = -1.77245385090551602729816748334d0*z**3*0.125d0;

!*** r*(exp(Qm1) - exp(Qp1))
 pf1   = -z6*0.125d0*rQA**3 - z4*0.125d0*rQA;
 pf3   = -0.25d0*z2*rQA;
 oddpow = r1*(pf1+r2*pf3)

!*** r^2(exp(Qm1) + exp(Qp1))
 pfb   = 0.25d0*z4*rQA**2;
 integrandsp(1:2*pts+1,1) = r1*pfa*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr) + r2*pfb*(expQmr+expQpr)


!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pf   = -1.77245385090551602729816748334d0*z**3*QA*0.25d0;

!*** r^odd*(exp(Qm1) - exp(Qp1))
 pf1   = 9*z8*0.03125d0*rQA**4 + z6*0.125d0*rQA**2 - z4*0.25d0;
 pf3   = z4*0.5d0*rQA**2;
 oddpow = r1*(pf1+r2*pf3) 

!*** r^even*(exp(Qm1) + exp(Qp1))
 pf2   = -9*z6*0.0625d0*rQA**3 - z4*0.25d0*rQA;
 pf4   = -z2*0.25d0*rQA;
 evenpow = r2*(pf2+r2*pf4)
 integrandsp(1:2*pts+1,2) =  r1*pf*(erfQmr-erfQpr)+ oddpow*(expQmr-expQpr) + evenpow*(expQmr+expQpr)


!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pf   = -3*1.77245385090551602729816748334d0*z**3*0.125d0*QA**2;

!*** r^(exp(Qm1) - exp(Qp1))
 pf1   = -15*z10*0.0625d0*rQA**5 - 9*z8*0.03125d0*rQA**3 + 3*z6*0.0625d0*rQA - 3*QA*z4*0.125d0;
 pf3   = -27*z6*0.0625d0*rQA**3 - 3*z4*0.125d0*rQA;
 pf5   = -z2*0.25d0*rQA;
 oddpow = r1*(pf1+r2*(pf3+r2*pf5));

!*** r^2(exp(Qm1) + exp(Qp1))
 pf2   = 15*z8*0.125d0*rQA**4 + 9*z6*0.0625d0*rQA**2 - 3*z4*0.125d0;
 pf4   = 7*z4*0.125d0*rQA**2;
 evenpow = r2*(pf2+r2*pf4)
 integrandsp(1:2*pts+1,3) = r1*pf*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)


!*** r^(erf(sqrtzeta*Qm1) - erf(sqrtzeta*Qp1));
 pf   = -1.77245385090551602729816748334d0*z**3*QA**3*0.5d0;

!*** r^(exp(Qm1) - exp(Qp1))
 pf1   = 525*z12*0.0078125d0*rQA**6 + 15*z10*0.0625d0*rQA**4 -3*z8*0.125d0*rQA**2 + z6*0.25d0 &
       - QA**2*z4*0.5d0;
 pf3   = 15*z8*0.5d0*rQA**4 + 3*z6*0.5d0*rQA**2 - z4*0.5d0;
 pf5   = 11*z4*0.125d0*rQA**2 ;
 oddpow = r1*(pf1+r2*(pf3+r2*pf5));

!*** r^2(exp(Qm1) + exp(Qp1))
 pf2   = -525*z10*0.015625d0*rQA**5 - 15*z8*0.125d0*rQA**3 + z6*0.75d0*rQA - QA*z4*0.5d0;
 pf4   = -65*z6*0.0625d0*rQA**3 - z4*0.5d0*rQA;
 pf6   = -z2*0.25d0*rQA ;
 evenpow = r2*(pf2+r2*(pf4+r2*pf6));
 integrandsp(1:2*pts+1,4) = r1*pf*(erfQmr-erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)

!*** r^2(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pfa   = -1.77245385090551602729816748334d0*z**3*0.125d0*rQA**2;

!*** r^2(exp(zeta*Qm1^2) + exp(zeta*Qp1^2));
 pfb   = z2*0.25d0*rQA;
 integrandsp(1:2*pts+1,0) = integrandsp(1:2*pts+1,0)+  r2*(pfa*(erfQmr+erfQpr) + pfb*(expQmr+expQpr))


!*** r^4 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pfa   = -1.77245385090551602729816748334d0*z**3*0.25d0*rQA**3;

!*** r^3 *(exp(Qm1) - exp(Qp1)
 pfb   = z4*0.25d0*rQA**3+ z2*0.25d0*rQA;
 integrandsp(1:2*pts+1,1) = integrandsp(1:2*pts+1,1)+ pfa*r4*(erfQmr+erfQpr) + pfb*r3*(expQmr-expQpr);



!*** r^6 (erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pfa   = -3*1.77245385090551602729816748334d0*z**3*0.125d0*rQA**4;

!*** r^3 *(exp(Qm1) - exp(Qp1)
 pf3   = -3*z6*0.0625d0*rQA**4 - z4*0.125d0*rQA**2;
 pf5   = 3*z4*0.125d0*rQA**4;
 oddpow = r3*(pf3+r2*pf5)

!*** r^4 *(exp(Qm1) + exp(Qp1)
 pfb   = 3*z4*0.125d0*rQA**3 + z2*0.25d0*rQA ;
 integrandsp(1:2*pts+1,2) = integrandsp(1:2*pts+1,2) + pfa*r6*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr)+ r4*pfb*(expQmr+expQpr)


!*** r^8(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pf   = -1.77245385090551602729816748334d0*z**3*0.5d0*rQA**5;

!*** r^3 *(exp(Qm1) - exp(Qp1)
 pf3   = 3*z8*0.125d0*rQA**5+3*z6*0.0625d0*rQA**3 ;
 pf5   = -z6*0.25d0*rQA**5 + z4*0.5d0*rQA**3 + z2*0.25d0*rQA;
 pf7   = z4*0.5d0*rQA**5;
 oddpow = r3*(pf3+r2*(pf5+r2*pf7));

!*** r^4 *(exp(Qm1) + exp(Qp1)
 pf4   = - 3*z6*0.25d0*rQA**4 - 3*z4*0.125d0*rQA**2;
 pf6   = pf7*QA!  z4*0.5d0*rQA**4;
 evenpow = r4*(pf4+r2*pf6)
 integrandsp(1:2*pts+1,3) = integrandsp(1:2*pts+1,3) + pf*r8*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr)+ evenpow*(expQmr+expQpr)

!*** r^10(erf(sqrtzeta*Qm1) + erf(sqrtzeta*Qp1));
 pf   = -5*1.77245385090551602729816748334d0*z**3*0.125d0*rQA**6;

!*** r^pow *(exp(Qm1) - exp(Qp1)
 pf3   = -75*z10*0.015625d0*rQA**6 - 15*z8*0.03125d0*rQA**4 ;
 pf5   = 15*z8*0.03125d0*rQA**6 - 15*z6*0.125d0*rQA**4 - 3*z4*0.25d0*rQA**2;
 pf7   = -5*z6*0.0625d0*rQA**6 + 5*z4*0.125d0*rQA**4;
 pf9   = 5*z4*0.125d0*rQA**6
 oddpow = r3*(pf3+r2*(pf5+r2*(pf7+r2*pf9))); 

!*** r^even *(exp(Qm1) + exp(Qp1)
 pf4   = 75*z8*0.03125d0*rQA**5 + 15*z6*0.0625d0*rQA**3 ;
 pf6   = -15*z6*0.0625d0*rQA**5 + 5*z4*0.125d0*rQA**3 + z2*0.25d0*rQA;
 pf8   = pf9*QA ! 5*z4*0.125d0*rQA**5;
 evenpow = r4*(pf4+r2*(pf6+r2*pf8));
 integrandsp(1:2*pts+1,4) = integrandsp(1:2*pts+1,4) + pf*r10*(erfQmr+erfQpr) + oddpow*(expQmr-expQpr) + evenpow*(expQmr+expQpr)


end subroutine calcspall





