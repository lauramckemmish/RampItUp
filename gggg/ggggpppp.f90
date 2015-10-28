subroutine calcggggpppp(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggpp,noggpp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
  implicit none
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun), basis(1:nobasisfun,1:30)
  real*8, intent(in)      :: ggpp(1:20,1:noggpp)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: noggpp, nobasisfun,out,compare
  character(len=30), intent(in)           :: moleculename
  integer  :: i, j, k, bf1, bf2, bf3, bf4, printopt, braketno, i1, i2, i3, i4, lengthtwoelec
  real*8   :: zeta, eta, pf, f0eval, alpha,beta,gamma,delta, expeval,RPQ2, nc1, nc2, gpf
  real*8  :: t1, tempint(1:3,1:3,1:3,1:3),  uaabb1, ubrabraket(1:3,1:3,1:3) 
  real*8  :: ubra2ketket(1:3,1:3,1:3),  ubrabraketket(1:3,1:3,1:3,1:3), PQ(1:3), AB(1:3), CD(1:3)
  real*8  :: u, uij, uijk, uiij, uijkl, uiijj, uiipart,uketket(1:3,1:3), RAB2, RCD2
  real*8  :: ubrabra(1:3,1:3), ubra1(1:3),ubra2(1:3), uket1(1:3),uket2(1:3), ubraint, uketint
  real*8  :: ubra1ketket(1:3,1:3,1:3), ubrabraket1(1:3,1:3,1:3), ubrabraket2(1:3,1:3,1:3)
  real*8  :: ubra1ket1(1:3,1:3),ubra1ket2(1:3,1:3),ubra2ket1(1:3,1:3),ubra2ket2(1:3,1:3)
  real*8  :: normtwoelecmat2(0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)
  real*8  :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)

#include "DEBUG.h"

#ifdef PSEUDOCHECK
normtwoelecmat2 = 0d0;
if (abs(Palpha(1,1)) .lt. 1d-10) then
call readQchemtwoeecompare(normtwoelecmat2,nobasisfun,basis,moleculename)
print *, "in ggggpppp";   compare = 0;   printopt = 0;
else;    compare = 0;   printopt = 0;   end if
#endif

  pf = (dacos(-1d0))**(5d0/2d0);

  do i=1,noggpp
    alpha = ggpp(6,i); beta = ggpp(7,i);    zeta = ggpp(8,i);
    bf1 = int(ggpp(1,i)); bf2 = int(ggpp(2,i));     nc1 = ggpp(9,i);   
    AB = basis(bf1,4:6)-basis(bf2,4:6);  RAB2 = AB(1)**2+AB(2)**2+AB(3)**2;
    do j=i,noggpp
      gamma = ggpp(6,j); delta = ggpp(7,j);eta = ggpp(8,j);
      bf3 = int(ggpp(1,j)); bf4 = int(ggpp(2,j));      nc2 = ggpp(9,j);
      PQ = ggpp(3:5,i)-ggpp(3:5,j);         RPQ2 = PQ(1)**2+PQ(2)**2+PQ(3)**2;
      CD = basis(bf3,4:6)-basis(bf4,4:6);   RCD2 = CD(1)**2+CD(2)**2+CD(3)**2;

      call calcf0(f0eval,zeta*eta/(zeta+eta)*RPQ2);      expeval = exp(-RPQ2*zeta*eta/(zeta+eta));

 tempint = 0d0;

  u = 2*pf/(zeta*eta*sqrt(zeta+eta))*f0eval;

  ubraint = 2*pf/(RPQ2*eta*zeta**2*sqrt(zeta+eta))*(expeval-f0eval);
  uketint = 2*pf/(RPQ2*zeta*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);
  ubra1 = ubraint*PQ*alpha;   ubra2 = ubraint*PQ*beta; 
  uket1 = uketint*PQ*gamma;   uket2 = uketint*PQ*delta; 

  uij = 2*pf*gamma*delta/(RPQ2**2*zeta*eta**3*sqrt(zeta+eta))*&
            (-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta)+ 3*f0eval);
  uiipart = -2*pf*gamma*delta/(RPQ2*zeta*eta**3*sqrt(zeta+eta))*(-expeval+f0eval);
  do i1=1,3; do i2=1,3;
   uketket(i1,i2) = uij*PQ(i1)*PQ(i2)
   if (i1 == i2) then;     uketket(i1,i2) = uketket(i1,i2) + uiipart;   end if
  end do; end do

  uij = -2*pf/(RPQ2**2*zeta**2*eta**2*sqrt(zeta+eta))*(-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta)+ 3*f0eval);
  uiipart = 2*pf/(RPQ2*zeta**2*eta**2*sqrt(zeta+eta))*(-expeval+f0eval);
  do i1=1,3; do i2=1,3;
   ubra1ket1(i1,i2) = alpha*gamma*uij*PQ(i1)*PQ(i2);   ubra2ket1(i1,i2) = beta*gamma*uij*PQ(i1)*PQ(i2);
   ubra1ket2(i1,i2) = alpha*delta*uij*PQ(i1)*PQ(i2);   ubra2ket2(i1,i2) = beta*delta*uij*PQ(i1)*PQ(i2);
   if (i1 == i2) then;
      ubra1ket1(i1,i2) = ubra1ket1(i1,i2) + uiipart*alpha*gamma
      ubra2ket1(i1,i2) = ubra2ket1(i1,i2) + uiipart*beta*gamma
      ubra1ket2(i1,i2) = ubra1ket2(i1,i2) + uiipart*alpha*delta
      ubra2ket2(i1,i2) = ubra2ket2(i1,i2) + uiipart*beta*delta
   end if;  end do; end do

  uij = 2*pf*alpha*beta/(RPQ2**2*eta*zeta**3*sqrt(zeta+eta))* (-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta)+ 3*f0eval);
  uiipart = -2*pf*alpha*beta/(RPQ2*eta*zeta**3*sqrt(zeta+eta))*(-expeval+f0eval);

  do i1=1,3; do i2=1,3;
   ubrabra(i1,i2) = uij*PQ(i1)*PQ(i2)
   if (i1 == i2) then;      ubrabra(i1,i2) = ubrabra(i1,i2) + uiipart;   end if; 
  end do; end do

  uijk = 2*pf*gamma*delta/(RPQ2**3*zeta**2*eta**3*sqrt(zeta+eta))*(-15*f0eval + (expeval/(zeta+eta)**2)*&
         (15*eta**2 + 10*eta*(3 + eta*RPQ2)*zeta + (15 + 2*eta*RPQ2*(5 + 2*eta*RPQ2))*zeta**2))
  uiij = 2*pf*gamma*delta/(RPQ2**2*zeta**2*eta**3*sqrt(zeta+eta))* (3*f0eval-expeval/(zeta+eta)*(3*(eta+zeta)+2*zeta*eta*RPQ2));
 do i1=1,3; do i2=1,3; do i3=1,3;
  ubra1ketket(i1,i2,i3) = uijk*PQ(i1)*PQ(i2)*PQ(i3)*alpha;
  ubra2ketket(i1,i2,i3) = uijk*PQ(i1)*PQ(i2)*PQ(i3)*beta; 
  if (i1 == i2) then
     ubra1ketket(i1,i2,i3) = ubra1ketket(i1,i2,i3) + uiij*PQ(i3)*alpha
     ubra2ketket(i1,i2,i3) = ubra2ketket(i1,i2,i3) + uiij*PQ(i3)*beta
  end if;
  if (i1 == i3) then;
     ubra1ketket(i1,i2,i3) = ubra1ketket(i1,i2,i3) + uiij*PQ(i2)*alpha
     ubra2ketket(i1,i2,i3) = ubra2ketket(i1,i2,i3) + uiij*PQ(i2)*beta
  end if;
  if (i2 == i3) then
     ubra1ketket(i1,i2,i3) = ubra1ketket(i1,i2,i3) + uiij*PQ(i1)*alpha
     ubra2ketket(i1,i2,i3) = ubra2ketket(i1,i2,i3) + uiij*PQ(i1)*beta
  end if; end do ; end do; end do

  uijk = -2*pf*alpha*beta/(RPQ2**3*zeta**3*eta**2*sqrt(zeta+eta))*(-15*f0eval + (expeval/(zeta+eta)**2)*&
         (15*(zeta+eta)**2 + 10*RPQ2*zeta*eta*(zeta+eta) + 4*RPQ2**2*zeta**2*eta**2))
  uiij = -2*pf*alpha*beta/(RPQ2**2*zeta**3*eta**2*sqrt(zeta+eta))*(3*f0eval-expeval/(zeta+eta)*(3*(eta+zeta)+2*zeta*eta*RPQ2));
 do i1=1,3; do i2=1,3; do i3=1,3;
  ubrabraket1(i1,i2,i3) = uijk*PQ(i1)*PQ(i2)*PQ(i3)*gamma;
  ubrabraket2(i1,i2,i3) = uijk*PQ(i1)*PQ(i2)*PQ(i3)*delta; 
  if (i1 == i2) then
     ubrabraket1(i1,i2,i3) = ubrabraket1(i1,i2,i3) + uiij*PQ(i3)*gamma
     ubrabraket2(i1,i2,i3) = ubrabraket2(i1,i2,i3) + uiij*PQ(i3)*delta
  end if
  if (i1 == i3) then;
     ubrabraket1(i1,i2,i3) = ubrabraket1(i1,i2,i3) + uiij*PQ(i2)*gamma
     ubrabraket2(i1,i2,i3) = ubrabraket2(i1,i2,i3) + uiij*PQ(i2)*delta
  end if
  if (i2 == i3) then
     ubrabraket1(i1,i2,i3) = ubrabraket1(i1,i2,i3) + uiij*PQ(i1)*gamma
     ubrabraket2(i1,i2,i3) = ubrabraket2(i1,i2,i3) + uiij*PQ(i1)*delta
  end if; end do ; end do; end do


 uiijj = 2*pf/(RPQ2**3*eta**3*zeta**3*sqrt(zeta+eta))*&
                  (-15*f0eval+(expeval/(zeta+eta)**2)*(15*(zeta+eta)**2+2*RPQ2*zeta*eta*(5*(zeta+eta)+2*zeta*eta*RPQ2)));
 t1 = 105+70*RPQ2*zeta*eta/(zeta+eta)+28*RPQ2**2*zeta**2*eta**2/(zeta+eta)**2+8*zeta**3*eta**3*RPQ2**3/(zeta+eta)**3;
 uijkl = 2*pf/(RPQ2**4*eta**3*zeta**3*sqrt(zeta+eta))*(105*f0eval-(expeval)*t1);
uaabb1 = 2*pf/(RPQ2**2*zeta**3*eta**3*sqrt(zeta+eta))*(3*f0eval-expeval/(zeta+eta)*(3*(zeta+eta)+2*RPQ2*zeta*eta));
if (RPQ2 .lt. 3d-5 .and. RPQ2 .gt. 1d-10) then
  uaabb1 = 16d0*pf*(1d0/(10d0*zeta*eta*(zeta+eta)**(2.5d0)) - RPQ2/(14*(zeta+eta)**(3.5d0))&
                      +RPQ2**2*zeta*eta/(36d0*(zeta+eta)**4.5d0) - RPQ2**3*zeta**2*eta**2/(132d0*(zeta+eta)**5.5d0))
  uiijj = -16d0*pf*(-1d0/(7d0*(zeta+eta)**3.5d0)+ zeta*eta*RPQ2/(9*(zeta+eta)**4.5d0)&
                      - zeta**2*eta**2*RPQ2**2/(22d0*(zeta+eta)**5.5d0) +&
                        RPQ2**3*zeta**3*eta**3/(78d0*(zeta+eta)**6.5d0))
  uijkl = 16d0*pf*(2d0*zeta*eta/(9*(zeta+eta**4.5d0))- 2*zeta**2*eta**2*RPQ2/(11*(zeta+eta)**5.5d0)&
                      + RPQ2**2*zeta**3*eta**3/(13d0*(zeta+eta)**7.5d0) &
                     - RPQ2**3*zeta**4*eta**4/(45d0*(zeta+eta)**8.5d0))
end if
 ubrabraketket = 0d0; 
 do i1=1,3; do i2=1,3;
     ubrabraketket(i1,i1,i2,i2) = ubrabraketket(i1,i1,i2,i2) + uaabb1
     ubrabraketket(i1,i2,i1,i2) = ubrabraketket(i1,i2,i1,i2) + uaabb1
     ubrabraketket(i1,i2,i2,i1) = ubrabraketket(i1,i2,i2,i1) + uaabb1
 end do; end do;


 do i1=1,3; do i2=1,3; do i3=1,3; do i4=1,3;
   ubrabraketket(i1,i2,i3,i4) = ubrabraketket(i1,i2,i3,i4) + PQ(i1)*PQ(i2)*PQ(i3)*PQ(i4)*uijkl
 end do; end do; end do; end do;


 do i1=1,3; do i2=1,3; do i3=1,3; do i4=1,3;
  if (i1 == i2) then;      ubrabraketket(i1,i2,i3,i4) = ubrabraketket(i1,i2,i3,i4) + PQ(i3)*PQ(i4)*uiijj;  end if
  if (i1 == i3) then;      ubrabraketket(i1,i2,i3,i4) = ubrabraketket(i1,i2,i3,i4) + PQ(i2)*PQ(i4)*uiijj;  end if
  if (i1 == i4) then;      ubrabraketket(i1,i2,i3,i4) = ubrabraketket(i1,i2,i3,i4) + PQ(i2)*PQ(i3)*uiijj;  end if
  if (i2 == i3) then;      ubrabraketket(i1,i2,i3,i4) = ubrabraketket(i1,i2,i3,i4) + PQ(i1)*PQ(i4)*uiijj;  end if
  if (i2 == i4) then;      ubrabraketket(i1,i2,i3,i4) = ubrabraketket(i1,i2,i3,i4) + PQ(i1)*PQ(i3)*uiijj;  end if
  if (i3 == i4) then;      ubrabraketket(i1,i2,i3,i4) = ubrabraketket(i1,i2,i3,i4) + PQ(i1)*PQ(i2)*uiijj;  end if
 end do; end do; end do; end do;


 if (zeta*eta/(zeta+eta)*RPQ2 .lt. 1d-8) then
    do i1=1,3; do i2=1,3; do i3=1,3; do i4=1,3;
      tempint(i1,i2,i3,i4) = 2*AB(i1)*AB(i2)*CD(i3)*CD(i4)*pf*alpha*beta*gamma*delta/(zeta**3*eta**3*sqrt(zeta+eta));
    if (i1 == i2) then;   
        tempint(i1,i2,i3,i4)=tempint(i1,i2,i3,i4)- CD(i3)*CD(i4)*pf*gamma*delta*(3*zeta+2*eta)/(3*zeta**2*eta**3*(zeta+eta)**1.5d0); 
    end if
    if (i3 == i4) then;   
        tempint(i1,i2,i3,i4)=tempint(i1,i2,i3,i4)- AB(i1)*AB(i2)*pf*alpha*beta*(2*zeta+3*eta)/(3*zeta**3*eta**2*(zeta+eta)**1.5d0); 
    end if
    if (i1 == i3) then;   
        tempint(i1,i2,i3,i4)=tempint(i1,i2,i3,i4)+AB(i2)*CD(i4)*pf*alpha*gamma/(3*zeta**2*eta**2*(zeta+eta)**1.5d0);    end if
    if (i1 == i4) then; 
        tempint(i1,i2,i3,i4)=tempint(i1,i2,i3,i4)-AB(i2)*CD(i3)*pf*alpha*delta/(3*zeta**2*eta**2*(zeta+eta)**1.5d0);    end if
    if (i2 == i3) then;   
        tempint(i1,i2,i3,i4)=tempint(i1,i2,i3,i4)-AB(i1)*CD(i4)*pf*beta*gamma/(3*zeta**2*eta**2*(zeta+eta)**1.5d0);    end if;
    if (i2 == i4) then;   
        tempint(i1,i2,i3,i4)=tempint(i1,i2,i3,i4)+AB(i1)*CD(i3)*pf*beta*delta/(3*zeta**2*eta**2*(zeta+eta)**1.5d0);    end if;
    end do; end do; end do; end do;

     do i1=1,3; do i2=1,3
       tempint(i1,i1,i2,i2) = pf*(10*zeta**2+10*eta**2+23*zeta*eta)/(30*zeta**2*eta**2*(zeta+eta)**(2.5d0)) &
                   -AB(i1)**2*pf*alpha*beta*(2*zeta+3*eta)/(3*zeta**3*eta**2*(zeta+eta)**1.5d0) &
                   -CD(i2)**2*pf*gamma*delta*(3*zeta+2*eta)/(3*zeta**2*eta**3*(zeta+eta)**1.5d0) &
                   + 2*AB(i1)**2*CD(i2)**2*alpha*beta*gamma*delta*pf/(zeta**3*eta**3*sqrt(zeta+eta)) ;
       tempint(i1,i2,i1,i2) =  pf/(10*zeta*eta*(zeta+eta)**(2.5d0)) + &
                                (AB(i1)*CD(i1)*beta*delta+AB(i2)*CD(i2)*alpha*gamma)*pf/(3*zeta**2*eta**2*(zeta+eta)**(1.5d0)) & 
                                +2*AB(i1)*AB(i2)*CD(i1)*CD(i2)*alpha*beta*gamma*delta*pf/(zeta**3*eta**3*sqrt(zeta+eta));
       tempint(i1,i2,i2,i1) =  pf/(10*zeta*eta*(zeta+eta)**(2.5d0)) - &
                                (AB(i1)*CD(i1)*beta*gamma+AB(i2)*CD(i2)*alpha*delta)*pf/(3*zeta**2*eta**2*(zeta+eta)**(1.5d0)) & 
                                +2*AB(i1)*AB(i2)*CD(i1)*CD(i2)*alpha*beta*gamma*delta*pf/(zeta**3*eta**3*sqrt(zeta+eta));
     end do; end do;

     do i1=1,3;
     tempint(i1,i1,i1,i1) = (5*zeta+2*eta)*(5*eta+2*zeta)*pf/(30d0*zeta**2*eta**2*(zeta+eta)**(2.5d0))  &
                 -AB(i1)**2*pf*alpha*beta*(2*zeta+3*eta)/(3d0*zeta**3*eta**2*(zeta+eta)**(1.5d0)) & 
                 -CD(i1)**2*pf*gamma*delta*(3*zeta+2*eta)/(3d0*zeta**2*eta**3*(zeta+eta)**(1.5d0)) &
                 + 2*AB(i1)**2*CD(i1)**2*alpha*beta*gamma*delta*pf/(zeta**3*eta**3*sqrt(zeta+eta)) &
                 + pf*AB(i1)*CD(i1)*(alpha-beta)*(gamma-delta)/(3d0*zeta**2*eta**2*(zeta+eta)**(1.5d0)) 
     end do;

 else ! if (RPQ2 .gt. 1d-10) 


    do i1=1,3; do i2=1,3;      tempint(i1,i1,i2,i2) = tempint(i1,i1,i2,i2)+ u/(4*zeta*eta);    end do; end do;

    do i1=1,3; do i2=1,3; do i3=1,3; do i4=1,3;
   if (i3==i4) then
      tempint(i1,i2,i3,i4) = tempint(i1,i2,i3,i4) +  ubrabra(i1,i2)/(8d0*alpha*beta*eta)  &
               -AB(i1)*AB(i2)*alpha*beta*u/(2d0*zeta**2*eta) +(AB(i2)*ubra1(i1)-ubra2(i2)*AB(i1))/(4*zeta*eta) ;
   end if;
   if (i1==i2) then
      tempint(i1,i2,i3,i4) = tempint(i1,i2,i3,i4) + uketket(i3,i4)/(8d0*gamma*delta*zeta) &
                 - CD(i3)*CD(i4)*gamma*delta*u/(2d0*zeta*eta**2)  +(CD(i4)*uket1(i3)-uket2(i4)*CD(i3))/(4d0*zeta*eta) 
   end if
   end do; end do; end do; end do;

  do i1=1,3; do i2=1,3; do i3=1,3; do i4=1,3; 
    tempint(i1,i2,i3,i4) = tempint(i1,i2,i3,i4) + ubrabraketket(i1,i2,i3,i4)/(16d0) + &
          AB(i1)*AB(i2)*CD(i3)*CD(i4)*u*alpha*beta*gamma*delta/(zeta**2*eta**2) + &
          AB(i1)*AB(i2)*alpha*beta/(2d0*zeta**2*eta)*(CD(i3)*uket2(i4)-CD(i4)*uket1(i3)) + &
          CD(i3)*CD(i4)*gamma*delta/(2d0*eta**2*zeta)*(AB(i1)*ubra2(i2)-AB(i2)*ubra1(i1))  &
          -alpha*beta/(4d0*gamma*delta*zeta**2)*AB(i1)*AB(i2)*uketket(i3,i4) &
          -gamma*delta/(4d0*alpha*beta*eta**2)*CD(i3)*CD(i4)*ubrabra(i1,i2) &   
          + 1d0/(4d0*zeta*eta)*( AB(i1)*CD(i3)*ubra2ket2(i2,i4)-AB(i1)*CD(i4)*ubra2ket1(i2,i3) & 
            - AB(i2)*CD(i3)*ubra1ket2(i1,i4) + AB(i2)*CD(i4)*ubra1ket1(i1,i3)) +&
          (AB(i2)*ubra1ketket(i1,i3,i4) -AB(i1)*ubra2ketket(i2,i3,i4))/(8d0*gamma*delta*zeta) +&
          (CD(i4)*ubrabraket1(i1,i2,i3) -CD(i3)*ubrabraket2(i1,i2,i4))/(8d0*alpha*beta*eta)

   end do; end do; end do; end do;


 end if

if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; 
   if (i .ne. j) then;     nc2 = nc2*2d0;   end if; braketno = 1;        else; braketno = 2;        end if;

tempint = tempint*nc1*nc2

 do i1=1,3; do i2=1,3; do i3=1,3; do i4=1,3;
#ifdef PSEUDOCHECK
if (abs(Palpha(1,1)) .lt. 1d-10) then
  call makenewtwo(newtwo,nobasisfun,bf1+i1-1,bf2+i2-1,bf3+i3-1,bf4+i4-1,tempint(i1,i2,i3,i4),braketno)
end if
#endif

 call digestintsnewcompare(tempint(i1,i2,i3,i4),bf1+i1-1,bf2+i2-1,bf3+i3-1,bf4+i4-1,& 
         Galpha,Gbeta,Ptot,Palpha,Pbeta, nobasisfun,printopt,braketno, normtwoelecmat2,out,compare)

 end do; end do; end do; end do;  


    end do
  end do

#ifdef INTPRINT
if (abs(Palpha(1,1)) .lt. 1d-10) then
!print *, "TESTING DIFFERENCES"
do bf1=1,nobasisfun; do bf2=bf1,nobasisfun; do bf3=bf1,nobasisfun; do bf4=bf2,nobasisfun;
   if (abs(newtwo(bf1,bf2,bf3,bf4,1)) .gt. 1d-10) then;
 if ( abs(newtwo(bf1,bf2,bf3,bf4,1)- normtwoelecmat2(bf1,bf2,bf3,bf4,1)) .gt. 1d-9) then
!    print *, bf1,bf2,bf3,bf4, newtwo(bf1,bf2,bf3,bf4,1), normtwoelecmat2(bf1,bf2,bf3,bf4,1)
 end if
   end if 
end do; end do; end do; end do;
!print *, "END OF TESTING DIFFERENCES"
end if
#endif

end subroutine calcggggpppp





!********

subroutine digestintsnewcompare(intt, bf1, bf2, bf3, bf4, Galpha,Gbeta,Ptot,Palpha,Pbeta,&
      nobasisfun,printopt,braketno, normtwoelecmat2, out,compare)
 implicit none
 integer, intent(in)   :: bf1, bf2, bf3, bf4, nobasisfun,printopt, compare
 real*8, intent(in)    :: intt
 real*8                :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
 real*8                :: Ptot(1:nobasisfun,1:nobasisfun), Galpha(1:nobasisfun,1:nobasisfun)
 real*8                :: Gbeta(1:nobasisfun,1:nobasisfun), x, y
 real*8  :: normtwoelecmat2(0:nobasisfun,0: nobasisfun,0: nobasisfun,0: nobasisfun,1:1)
 integer :: braketno,out

     Galpha(bf1,bf2) = Galpha(bf1,bf2) + 2d0*Ptot(bf3,bf4)* intt*0.5d0*braketno
     Gbeta(bf1,bf2) = Gbeta(bf1,bf2) + 2d0*Ptot(bf3,bf4)* intt*0.5d0*braketno

     Galpha(bf1,bf4) = Galpha(bf1,bf4) - 0.5d0*braketno*Palpha(bf2,bf3)* intt 
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - 0.5d0*braketno*Pbeta(bf2,bf3)* intt
     Galpha(bf1,bf3) = Galpha(bf1,bf3) - 0.5d0*braketno*Palpha(bf2,bf4)* intt
     Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - 0.5d0*braketno*Pbeta(bf2,bf4)* intt
     Galpha(bf2,bf4) = Galpha(bf2,bf4) - 0.5d0*braketno*Palpha(bf1,bf3)* intt
     Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - 0.5d0*braketno* Pbeta(bf1,bf3)* intt
     Galpha(bf2,bf3) = Galpha(bf2,bf3) - 0.5d0*braketno*Palpha(bf1,bf4)* intt
     Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - 0.5d0*braketno*Pbeta(bf1,bf4)* intt

     Galpha(bf3,bf4) = Galpha(bf3,bf4) + 2d0*Ptot(bf1,bf2)* intt*0.5d0*braketno
     Gbeta(bf3,bf4) = Gbeta(bf3,bf4) + 2d0*Ptot(bf1,bf2)* intt*0.5d0*braketno


end subroutine digestintsnewcompare




