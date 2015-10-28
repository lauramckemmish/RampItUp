
!************************************
subroutine calcfullang(ang,QAx,QAy,QAz)
  implicit none
  real*8, intent(in)  :: QAx, QAy, QAz
  real*8, intent(out) :: ang(1:25)
  real*8              :: pi, QA
  integer             :: l
  real*8                :: costh, sinth(0:4), cosph, sinph, cosmph(0:4), sinmph(0:4)
 
  sinth(0) = 1d0; 
  !***** Calculating the long-range angular component!!
!  theta = datan2(sqrt(QAx**2+QAy**2),QAz);       phi   = datan2(QAy , QAx);
  QA = sqrt(QAx**2+QAy**2+QAz**2);
  costh = QAz/QA;
  sinth(1) = sqrt(QAx**2+QAy**2)/QA;
  cosph = QAx/(QA*sinth(1));
  sinph = QAy/(QA*sinth(1));
!  costh = cos(theta); sinth(1) = sin(theta); cosph = cos(phi); sinph = sin(phi)
  do l=2,4; sinth(l) = sinth(l-1)*sinth(1); end do;
  cosmph(0) = 1;
  cosmph(1) = cosph;
  sinmph(0) = 0; 
  sinmph(1) = sinph;

  do l=2,4;
   cosmph(l) = 2*cosph*cosmph(l-1)-cosmph(l-2); 
   sinmph(l) = 2*cosph*sinmph(l-1)-sinmph(l-2); 
  end do

  call calcang(ang,costh,sinth,cosph,sinph, cosmph,sinmph)
  
  if (abs(QAx) .lt. 1d-10 .and. abs(QAy) .lt. 1d-10) then;
    call QAxyzeroang(ang,QAz,4)
  end if;  

  return
end subroutine calcfullang
!*******************************************************************

subroutine calcderthph(dangdx,dangdy,dangdz,ang,QAx,QAy,QAz,QA,QA2, dangdth,dangdph)
  real*8, intent(in)        :: QAx, QAy, QAz, QA, QA2
  real*8, intent(out)       :: ang(1:25), dangdx(1:25)
  real*8, intent(out)       :: dangdy(1:25), dangdz(1:25)
  real*8                    :: theta,phi, dangdth(1:25), dangdph(1:25)
  real*8                    :: dphdx, dphdy, dthdx,dthdy,dthdz, pi
  real*8                    :: costh, sinth(0:4), cosph, sinph, v0, v1, v2, cosmph(0:4), sinmph(0:4)
  integer                   :: l
 
  costh = QAz/QA;
  sinth(1) = sqrt(QAx**2+QAy**2)/QA;
  cosph = QAx/(QA*sinth(1));
  sinph = QAy/(QA*sinth(1));

  do l=2,4; sinth(l) = sinth(l-1)*sinth(1); end do;
  cosmph(0) = 1;
  cosmph(1) = cosph;
  sinmph(0) = 0; 
  sinmph(1) = sinph;

  do l=2,4;
   cosmph(l) = 2*cosph*cosmph(l-1)-cosmph(l-2); 
   sinmph(l) = 2*cosph*sinmph(l-1)-sinmph(l-2); 
  end do

  call calcang(ang,costh,sinth,cosph,sinph, cosmph,sinmph)
   

  sinth(0) = 1d0;
v0 = sqrt(QAx**2+QAy**2)
v1 = 1d0/v0
v2 = 1d0/QA2;
!*** Calculates the orientation prefactors
  
  call calcdangdth(dangdth,costh,sinth,cosph,sinph,cosmph,sinmph)
  call calcdangdph(dangdph,costh,sinth,cosph,sinph, cosmph,sinmph)
  
  dthdx = (QAx*QAz)*v1*v2
  dthdy = (QAy*QAz)*v1*v2
  dthdz = -v0*v2
  
  dphdx = -QAy*v1**2;
  dphdy =  QAx*v1**2;

  dangdx = dangdth*dthdx + dangdph*dphdx
  dangdy = dangdth*dthdy + dangdph*dphdy
  dangdz = dangdth*dthdz 

  if (abs(QAx) .lt. 1d-10 .and. abs(QAy) .lt. 1d-10) then;
   call QAxyzeroang(ang,QAz,4)
    call QAxyzerodang(dangdx,dangdy,dangdz,QAz,4)
  end if;  

  return
end subroutine calcderthph

!********************************************************
subroutine calcder2thph(ang,dangdx,dangdy,dangdz,d2angdx2,d2angdxy,d2angdxz,d2angdy2,d2angdyz,d2angdz2,QAx,QAy,QAz,QA2, &
        dangdth,dangdph,d2angdth2,d2angdph2,d2angdthph)
  implicit none
  real*8, intent(in)    :: QAx, QAy, QAz, QA2
  real*8, intent(out)   :: ang(1:25), dangdx(1:25), dangdy(1:25),dangdz(1:25)
  real*8, intent(out)   :: d2angdx2(1:25), d2angdxy(1:25), d2angdxz(1:25),d2angdy2(1:25)
  real*8, intent(out)   :: d2angdyz(1:25), d2angdz2(1:25)
  real*8                :: theta, phi, pi
  real*8                :: d2angdth2(1:25), d2angdph2(1:25),d2angdthph(1:25)
  real*8                :: dphdx, dphdy, dthdx,dthdy,dthdz
  real*8                :: dangdth(1:25), dangdph(1:25)
  real*8                :: d2thdx2, d2thdy2, d2thdz2, d2thdxy, d2thdxz, d2thdyz
  real*8                :: d2phdx2, d2phdy2, d2phdxy, QA
  real*8                :: costh, sinth(0:4), cosph, sinph, v0, v1, v2, cosmph(0:4), sinmph(0:4)
  integer               :: l
  

v0 = sqrt(QAx**2+QAy**2)
v1 = 1d0/v0
v2 = 1d0/QA2;
  sinth(0) = 1d0;
    QA = sqrt(QA2)
  costh = QAz/QA;
  sinth(1) = sqrt(QAx**2+QAy**2)/QA;
  cosph = QAx/(QA*sinth(1));
  sinph = QAy/(QA*sinth(1));

  do l=2,4; sinth(l) = sinth(l-1)*sinth(1); end do;
  cosmph(0) = 1;
  cosmph(1) = cosph;
  sinmph(0) = 0; 
  sinmph(1) = sinph;

  do l=2,4;
   cosmph(l) = 2*cosph*cosmph(l-1)-cosmph(l-2); 
   sinmph(l) = 2*cosph*sinmph(l-1)-sinmph(l-2); 
  end do

  call calcang(ang,costh,sinth,cosph,sinph, cosmph,sinmph)

  call calcdangdth(dangdth,costh,sinth,cosph,sinph,cosmph,sinmph)
  call calcdangdph(dangdph,costh,sinth,cosph,sinph, cosmph,sinmph)

  dthdx = (QAx*QAz)*v1*v2
  dthdy = (QAy*QAz)*v1*v2
  dthdz = -v0*v2
  
  dphdx = -QAy*v1**2;
  dphdy =  QAx*v1**2;

  dangdx = dangdth*dthdx + dangdph*dphdx
  dangdy = dangdth*dthdy + dangdph*dphdy
  dangdz = dangdth*dthdz 

  !*** Calculates the orientation prefactors
  call calcd2angdth2(d2angdth2,costh,sinth,cosph,sinph,cosmph,sinmph)
  call calcd2angdph2(d2angdph2,costh,sinth,cosph,sinph,cosmph,sinmph)
  call calcd2angdthph(d2angdthph,costh, sinth, cosph, sinph,cosmph,sinmph)
  
  d2thdx2 = -QAz*(2*QAx**4+QAx**2*QAy**2-QAy**2*( QAy**2+QAz**2))*v1**3*v2**2
  d2thdy2 = -QAz*(2*QAy**4-QAx**4       -QAx**2*(-QAy**2+QAz**2))*v1**3*v2**2
  d2thdz2 = 2*v0*QAz*v2**2
  
  d2thdxy = -QAx*QAy*QAz*(3*(QAx**2+QAy**2)+QAz**2)*v1**3*v2**2
  d2thdxz = QAx*(QAx**2+QAy**2-QAz**2)*v1*v2**2
  d2thdyz = QAy*(QAx**2+QAy**2-QAz**2)*v1*v2**2
  
  d2phdx2 = 2*QAx*QAy*v1**4
  d2phdy2 = -d2phdx2
  
  d2phdxy = -(QAx**2-QAy**2)*v1**4
 
  d2angdx2 = d2phdx2*dangdph + (dphdx)**2*d2angdph2 + d2thdx2*dangdth + (dthdx)**2*d2angdth2 + &
       2*dphdx*dthdx*d2angdthph
  
  d2angdy2 = d2phdy2*dangdph + (dphdy)**2*d2angdph2 + d2thdy2*dangdth + (dthdy)**2*d2angdth2 + &
       2*dphdy*dthdy*d2angdthph
  
  d2angdz2 = d2thdz2*dangdth + (dthdz)**2*d2angdth2 
  
  d2angdxy = d2phdxy*dangdph + dphdx*dphdy*d2angdph2 + d2thdxy*dangdth + dthdx*dthdy*d2angdth2 + &
       dphdx*dthdy*d2angdthph +  dphdy*dthdx*d2angdthph
  
  d2angdxz = d2thdxz*dangdth + dthdx*dthdz*d2angdth2 +  dphdx*dthdz*d2angdthph 
  
  d2angdyz =  d2thdyz*dangdth + dthdz*(dthdy*d2angdth2 + dphdy*d2angdthph)
  
  
  if (abs(QAx) .lt. 1d-10 .and. abs(QAy) .lt. 1d-10) then;
    call QAxyzeroang(ang,QAz,4)
    call QAxyzerodang(dangdx,dangdy,dangdz,QAz,4)
    call QAxyzerod2ang(d2angdx2,d2angdy2,d2angdz2,d2angdxy,d2angdxz,d2angdyz,QAz,4)
  end if;  

  return
end subroutine calcder2thph


!*********************************************************

!***************************************************************

subroutine calcang(ang,costh, sinth,cosph,sinph, cosmph,sinmph) 
  implicit none
  real*8,intent(out)           :: ang(1:25)
  integer :: l
  real*8                       :: costh, sinth(0:4), cosph,sinph, cosmph(0:4), sinmph(0:4), pf(0:5)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calking calculates the orientation dependence of the integral
  ! = RY_{l,m} [theta, phi]
  !
  ! Speed increases could be made Qy precomputing powers of sin/ cos th/phi etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 pf(1) = -4.886025119029199d-1; pf(2) = 0.54627421529603953527d0;
 pf(3) = -0.59004358992664351035d0; pf(4) = 0.62583573544917613459d0; 
 pf(5) = -0.65638205684017010281d0;

  do L=1,4
     pf(L)  = pf(L)*sinth(L)
     ang((L+1)**2) = pf(L)*cosmph(L)
     ang(L**2+1)   = pf(L)*sinmph(L)
  end do

  ang(1)  = 2.8209479177387813d-1
    ang(3)  = 4.886025119029199d-1*costh

  ang(6)  = -1.0925484305920792d0*costh*sinph*sinth(1)
  ang(7)  = -3.1539156525251997d-1+9.461746957575599d-1*costh**2
  ang(8)  = -1.0925484305920792d0*cosph*costh*sinth(1)
  
  ang(11)  = 1.4453057213202770277d0*costh*sinth(2)*sinmph(2)
  ang(12)  = 0.45704579946446573616d0*sinph*sinth(1)*(1-5*costh**2)
  ang(13)  = 0.37317633259011539141d0*costh*(-3+5*costh**2)
  ang(14)  = 0.45704579946446573616d0*cosph*sinth(1)*(1-5*costh**2)
  ang(15)  = 1.4453057213202770277d0*costh*sinth(2)*cosmph(2)
  
  ang(18)  = -1.7701307697799305310d0*costh*sinth(3)*sinmph(3)
  ang(19)  = 0.47308734787878000905d0*sinth(2)*sinmph(2)*(-1+7*costh**2)
  ang(20)  = 2.0071396306718675039d0*sinph*sinth(1)*costh*(8*sinth(2)-7)
  ang(21)  = 0.10578554691520430380d0*(3+costh**2*(-30+35*costh**2))
  ang(22)  = 2.0071396306718675039d0*cosph*sinth(1)*costh*(1-7*costh**2)
  ang(23)  = 0.47308734787878000905d0*sinth(2)*cosmph(2)*(-1+7*costh**2)
  ang(24)  = -1.7701307697799305310d0*costh*sinth(3)*cosmph(3)


 return
end subroutine calcang





!***************************************************************
subroutine calcdangdth(dangdth,costh,sinth,cosph,sinph,cosmph,sinmph) 
  implicit none
  real*8,intent(out)           :: dangdth(1:25)
  integer                      :: l
  real*8, intent(in)           :: costh, sinth(0:4), cosph,sinph
  real*8                       :: pf(1:5), sinmph(0:4), cosmph(0:4)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calking calculates the orientation dependence of the integral
! = RY_{l,m} [theta, phi]
!
! Speed increases could be made by precomputing powers of sin/ cos th/phi etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 pf(1) = -4.886025119029199d-1; pf(2) = 0.54627421529603953527d0;
 pf(3) = -0.59004358992664351035d0; pf(4) = 0.62583573544917613459d0; 
 pf(5) = -0.65638205684017010281d0;

  dangdth(1)  = 0

  do L=1,4
     pf(L)  = pf(L)*L*sinth(L-1)*costh
     dangdth((L+1)**2) = pf(L)*cosmph(L)
     dangdth(L**2+1)   = pf(L)*sinmph(L)
  end do

  dangdth(3)  = -4.886025119029199d-1*sinth(1)
  
  dangdth(6)  = 1.0925484305920790705d0*sinph*(-1+2*sinth(2))
  dangdth(7)  = -1.8923493915151202d0*costh*sinth(1)
  dangdth(8)  = 1.0925484305920790705d0*cosph*(-1+2*sinth(2))
  
  dangdth(11)  =  1.4453057213202770277d0*sinth(1)*sinph*(2-3*sinth(2))
  dangdth(12)  = 0.45704579946446573616d0*costh*sinph*(-4+15*sinth(2))
  dangdth(13)  = 0.37317633259011539141d0*(-12*sinth(1)+15*sinth(3))
  dangdth(14)  = 0.45704579946446573616d0*costh*cosph*(-4+15*sinth(2))
  dangdth(15)  = 1.4453057213202770277d0*sinth(1)*cosph*(2-3*sinth(2))

  
  dangdth(18)  = 7.0805230791197221241d0*sinmph(3)*(sinth(4)-3*sinth(2))
  dangdth(19)  = -1.8923493915151202d0*cosph*costh*sinph*sinth(1)+1.324644574060584d1*cosph*costh**3*sinph*sinth(1)-  &
       1.324644574060584d1*cosph*costh*sinph*sinth(3)
  dangdth(20)  = -3.345232717786446d-1*costh**2*sinph-2.341662902450512d0*costh**4*sinph+  &
       3.345232717786446d-1*sinph*sinth(2)+1.4049977414703072d1*costh**2*sinph*sinth(2)-  &
       2.341662902450512d0*sinph*sinth(4)
  dangdth(21)  = -1.057855469152043d0*costh*sinth(1)-7.404988284064301d0*costh**3*sinth(1)+  &
       7.404988284064301d0*costh*sinth(3)
  dangdth(22)  = -3.345232717786446d-1*cosph*costh**2-2.341662902450512d0*cosph*costh**4+  &
       3.345232717786446d-1*cosph*sinth(2)+1.4049977414703072d1*cosph*costh**2*sinth(2)-  &
       2.341662902450512d0*cosph*sinth(4)
  dangdth(23)  = -9.461746957575599d-1*cosph**2*costh*sinth(1)+6.623222870302921d0*cosph**2*costh**3*sinth(1)+  &
       9.461746957575599d-1*costh*sinph**2*sinth(1)-6.623222870302921d0*costh**3*sinph**2*sinth(1)-  &
       6.623222870302921d0*cosph**2*costh*sinth(3)+6.623222870302921d0*costh*sinph**2*sinth(3)
  dangdth(24)  = 7.0805230791197221241d0*cosmph(3)*(sinth(4)-3*sinth(2))
  return
end subroutine calcdangdth

!***************************************************************
subroutine calcdangdph(dangdph,costh,sinth,cosph,sinph,cosmph,sinmph) 
  implicit none
  real*8,intent(out)           :: dangdph(1:25)
  integer                      :: l
  real*8, intent(in)           :: costh, sinth(0:4), cosph,sinph,cosmph(0:4), sinmph(0:4)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calking calculates the orientation dependence of the integral
  ! = RY_{l,m} [theta, phi]
  !
  ! Speed increases could be made by precomputing powers of sin/ cos th/phi etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   

  do l=0,4
    dangdph(L**2+L+1) = 0d0;
  end do

  dangdph(2)  = -4.886025119029199d-1*cosph*sinth(1)
  dangdph(4)  = 4.886025119029199d-1*sinph*sinth(1)
  
  dangdph(5)  = 1.0925484305920790705d0*cosmph(2)*sinth(2)
  dangdph(6)  = -1.0925484305920792d0*cosph*costh*sinth(1)
  dangdph(8)  = 1.0925484305920792d0*costh*sinph*sinth(1)
  dangdph(9)  = -1.0925484305920792d0*sinmph(2)*sinth(2)
  
  dangdph(10)  = -1.7701307697799305310d0*cosmph(3)*sinth(3)
  dangdph(11)  = 2.8906114426405540554d0*costh*sinth(2)*cosmph(2)
  dangdph(12)  = 0.45704579946446573616d0*cosph*(-4*sinth(1)+5*sinth(3))
  dangdph(14)  = 0.45704579946446573616d0*sinph*( 4*sinth(1)-5*sinth(3))
  dangdph(15)  = -2.8906114426405540554d0*costh*sinth(2)*sinmph(2)
  dangdph(16)  = 1.7701307697799305310d0*sinmph(3)*sinth(3)
  
  dangdph(17)  = 2.5033429417967045383d0*cosmph(4)*sinth(4)
  dangdph(18)  = -5.3103923093397915931d0*costh*cosmph(3)*sinth(3)
  dangdph(19)  = 0.31539156525252000603d0*cosmph(2)*(18*sinth(2)-21*sinth(4))
  dangdph(20)  = 0.22301551451909638932d0*costh*cosph*(-12*sinth(1)+21*sinth(3))
  dangdph(22)  = -0.22301551451909638932d0*costh*sinph*(-12*sinth(1)+21*sinth(3))
  dangdph(23)  = -0.31539156525252000603d0*sinmph(2)*(18*sinth(2)-21*sinth(4))
  dangdph(24)  = 5.3103923093397915931d0*costh*sinmph(3)*sinth(3)
  dangdph(25)  = -2.5033429417967045383d0*sinth(4)*sinmph(4)

   return
end subroutine calcdangdph




!******************************************************************



!***************************************************************
subroutine calcd2angdth2(d2angdth2,costh,sinth,cosph,sinph,cosmph,sinmph) 
  implicit none
  real*8,intent(out)           :: d2angdth2(1:25)
  integer                      :: l
  real*8, intent(in)           :: costh, sinth(0:4), cosph,sinph,cosmph(0:4),sinmph(0:4)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calking calculates the orientation dependence of the integral
  ! = RY_{l,m} [theta, phi]
  !
  ! Speed increases could be made by precomputing powers of sin/ cos th/phi etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  d2angdth2(1)  = 0
  
  d2angdth2(2)  = 4.886025119029199d-1*sinph*sinth(1)
  d2angdth2(3)  = -4.886025119029199d-1*costh
  d2angdth2(4)  = 4.886025119029199d-1*cosph*sinth(1)
  
  d2angdth2(5)  = 1.0925484305920790705d0*sinmph(2)*(1-2*sinth(2))
  d2angdth2(6)  = 4.370193722368317d0*costh*sinph*sinth(1)
  d2angdth2(7)  = 1.8923493915151200362d0*(-1+2*sinth(2))
  d2angdth2(8)  = 4.370193722368317d0*cosph*costh*sinth(1)
  d2angdth2(9)  = 1.0925484305920790705d0*cosmph(2)*(1-2*sinth(2))
  
  d2angdth2(10)  = 0.59004358992664351035d0*sinmph(3)*(-6*sinth(1)+9*sinth(3))
  d2angdth2(11)  = -7.226528606601386d-1*cosph*costh*sinph+6.503875745941246d0*cosph*costh**3*sinph-  &
       1.9511627237823739d1*cosph*costh*sinph*sinth(2)
  d2angdth2(12)  = 1.1426144986611644d-1*sinph*sinth(1)+1.5425295731925717d1*costh**2*sinph*sinth(1)-  &
       5.1417652439752395d0*sinph*sinth(3)
  d2angdth2(13)  = -2.7988224944258655d-1*costh-4.1982337416387985d0*costh**3+1.2594701224916394d1*costh*sinth(2)
  d2angdth2(14)  = 1.1426144986611644d-1*cosph*sinth(1)+1.5425295731925717d1*cosph*costh**2*sinth(1)-  &
       5.1417652439752395d0*cosph*sinth(3)
  d2angdth2(15)  = -3.613264303300693d-1*cosph**2*costh+3.251937872970623d0*cosph**2*costh**3+  &
       3.613264303300693d-1*costh*sinph**2-3.251937872970623d0*costh**3*sinph**2-  &
       9.755813618911871d0*cosph**2*costh*sinth(2)+9.755813618911871d0*costh*sinph**2*sinth(2)
  d2angdth2(16)  = 0.59004358992664351035d0*cosmph(3)*(-6*sinth(1)+9*sinth(3))
  
  d2angdth2(17)  = 5.006685883593409d0*cosph**3*costh**2*sinph-5.006685883593409d0*cosph**3*costh**4*sinph-  &
       5.006685883593409d0*cosph*costh**2*sinph**3+5.006685883593409d0*cosph*costh**4*sinph**3-  &
       5.006685883593409d0*cosph**3*sinph*sinth(2)+3.0040115301560455d1*cosph**3*costh**2*sinph*sinth(2)+  &
       5.006685883593409d0*cosph*sinph**3*sinth(2)-3.0040115301560455d1*cosph*costh**2*sinph**3*sinth(2)-  &
       5.006685883593409d0*cosph**3*sinph*sinth(4)+5.006685883593409d0*cosph*sinph**3*sinth(4)
  d2angdth2(18)  = 1.0620784618679584d1*cosph**2*costh*sinph*sinth(1)-4.248313847471834d1*cosph**2*costh**3*sinph*sinth(1)-  &
       3.5402615395598613d0*costh*sinph**3*sinth(1)+1.4161046158239443d1*costh**3*sinph**3*sinth(1)+  &
       4.248313847471834d1*cosph**2*costh*sinph*sinth(3)-1.4161046158239443d1*costh*sinph**3*sinth(3)
  d2angdth2(19)  = -1.8923493915151202d0*cosph*costh**2*sinph+1.324644574060584d1*cosph*costh**4*sinph+  &
       1.8923493915151202d0*cosph*sinph*sinth(2)-7.947867444363505d1*cosph*costh**2*sinph*sinth(2)+  &
       1.324644574060584d1*cosph*sinph*sinth(4)
  d2angdth2(20)  = 1.3380930871145782d0*costh*sinph*sinth(1)+3.746660643920819d1*costh**3*sinph*sinth(1)-  &
       3.746660643920819d1*costh*sinph*sinth(3)
  d2angdth2(21)  = -1.057855469152043d0*costh**2-7.404988284064301d0*costh**4+1.057855469152043d0*sinth(2)+  &
       4.44299297043858d1*costh**2*sinth(2)-7.404988284064301d0*sinth(4)
  d2angdth2(22)  = 1.3380930871145782d0*cosph*costh*sinth(1)+3.746660643920819d1*cosph*costh**3*sinth(1)-  &
       3.746660643920819d1*cosph*costh*sinth(3)
  d2angdth2(23)  = -9.461746957575599d-1*cosph**2*costh**2+6.623222870302921d0*cosph**2*costh**4+  &
       9.461746957575599d-1*costh**2*sinph**2-6.623222870302921d0*costh**4*sinph**2+  &
       9.461746957575599d-1*cosph**2*sinth(2)-3.9739337221817523d1*cosph**2*costh**2*sinth(2)-  &
       9.461746957575599d-1*sinph**2*sinth(2)+3.9739337221817523d1*costh**2*sinph**2*sinth(2)+  &
       6.623222870302921d0*cosph**2*sinth(4)-6.623222870302921d0*sinph**2*sinth(4)
  d2angdth2(24)  = 3.5402615395598613d0*cosph**3*costh*sinth(1)-1.4161046158239443d1*cosph**3*costh**3*sinth(1)-  &
       1.0620784618679584d1*cosph*costh*sinph**2*sinth(1)+4.248313847471834d1*cosph*costh**3*sinph**2*sinth(1)+  &
       1.4161046158239443d1*cosph**3*costh*sinth(3)-4.248313847471834d1*cosph*costh*sinph**2*sinth(3)
  d2angdth2(25)  = 1.2516714708983523d0*cosph**4*costh**2-1.2516714708983523d0*cosph**4*costh**4-  &
       7.510028825390114d0*cosph**2*costh**2*sinph**2+7.510028825390114d0*cosph**2*costh**4*sinph**2+  &
       1.2516714708983523d0*costh**2*sinph**4-1.2516714708983523d0*costh**4*sinph**4-  &
       1.2516714708983523d0*cosph**4*sinth(2)+7.510028825390114d0*cosph**4*costh**2*sinth(2)+  &
       7.510028825390114d0*cosph**2*sinph**2*sinth(2)-  &
       4.506017295234068d1*cosph**2*costh**2*sinph**2*sinth(2)-1.2516714708983523d0*sinph**4*sinth(2)+  &
       7.510028825390114d0*costh**2*sinph**4*sinth(2)-1.2516714708983523d0*cosph**4*sinth(4)+  &
       7.510028825390114d0*cosph**2*sinph**2*sinth(4)-1.2516714708983523d0*sinph**4*sinth(4)
  
  
  return
end subroutine calcd2angdth2

!***************************************************************


subroutine calcd2angdph2(d2angdph2,costh,sinth,cosph,sinph,cosmph,sinmph) 
 implicit none
 real*8,intent(out)           :: d2angdph2(1:25)
 integer                      :: l
 real*8, intent(in)           :: costh, sinth(0:4), cosph,sinph,cosmph(0:4),sinmph(0:4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calking calculates the orientation dependence of the integral
! = RY_{l,m} [theta, phi]
!
! Speed increases could be made by precomputing powers of sin/ cos th/phi etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do l=0,4
    d2angdph2(L**2+L+1) = 0d0;
  end do

    d2angdph2(2)  = 4.886025119029199d-1*sinph*sinth(1)
    d2angdph2(4)  = 4.886025119029199d-1*cosph*sinth(1)
 
    d2angdph2(5)  = -2.1850968611841581411d0*sinth(2)*sinmph(2)
    d2angdph2(6)  = 1.0925484305920792d0*costh*sinph*sinth(1)
    d2angdph2(8)  = 1.0925484305920792d0*cosph*costh*sinth(1)
    d2angdph2(9)  = -2.1850968611841581411d0*sinth(2)*cosmph(2)

    d2angdph2(10)  = 5.3103923093397915931d0*sinth(3)*sinmph(3)
    d2angdph2(11)  = -5.7812228852811081108d0*costh*sinth(2)*sinmph(2)
    d2angdph2(12)  = 1.1426144986611644d-1*sinph*sinth(1)+1.7139217479917466d0*costh**2*sinph*sinth(1)-  &
5.7130724933058215d-1*sinph*sinth(3)
    d2angdph2(14)  = 1.1426144986611644d-1*cosph*sinth(1)+1.7139217479917466d0*cosph*costh**2*sinth(1)-  &
5.7130724933058215d-1*cosph*sinth(3)
    d2angdph2(15)  = -5.7812228852811081108d0*costh*sinth(2)*cosmph(2)
    d2angdph2(16)  = 5.3103923093397915931d0*sinth(3)*cosmph(3)
 
    d2angdph2(17)  = -10.013371767186818153d0*sinth(4)*sinmph(4)
    d2angdph2(18)  = 15.931176928019374779d0*sinth(3)*costh*sinmph(3)
    d2angdph2(19)  = -1.4192620436363401d0*cosph*sinph-1.8923493915151202d0*cosph*costh**2*sinph+  &
3.3116114351514603d0*cosph*costh**4*sinph+1.8923493915151202d0*cosph*sinph*sinth(2)-  &
1.9869668610908762d1*cosph*costh**2*sinph*sinth(2)+3.3116114351514603d0*cosph*sinph*sinth(4)
    d2angdph2(20)  = 3.345232717786446d-1*costh*sinph*sinth(1)+2.341662902450512d0*costh**3*sinph*sinth(1)-  &
2.341662902450512d0*costh*sinph*sinth(3)
    d2angdph2(22)  = 3.345232717786446d-1*cosph*costh*sinth(1)+2.341662902450512d0*cosph*costh**3*sinth(1)-  &
2.341662902450512d0*cosph*costh*sinth(3)
    d2angdph2(23)  = -7.096310218181701d-1*cosph**2-9.461746957575599d-1*cosph**2*costh**2+  &
1.6558057175757301d0*cosph**2*costh**4+7.096310218181701d-1*sinph**2+  &
9.461746957575599d-1*costh**2*sinph**2-1.6558057175757301d0*costh**4*sinph**2+  &
9.461746957575599d-1*cosph**2*sinth(2)-9.93483430545438d0*cosph**2*costh**2*sinth(2)-  &
9.461746957575599d-1*sinph**2*sinth(2)+9.93483430545438d0*costh**2*sinph**2*sinth(2)+  &
1.6558057175757301d0*cosph**2*sinth(4)-1.6558057175757301d0*sinph**2*sinth(4)
    d2angdph2(24)  = 15.931176928019374779d0*sinth(3)*costh*cosmph(3)
    d2angdph2(25)  = -10.013371767186818153d0*sinth(4)*cosmph(4)

 
return
end subroutine calcd2angdph2

!***************************************************************



subroutine calcd2angdthph(d2angdthph,costh, sinth, cosph, sinph, cosmph,sinmph)
 implicit none
 real*8,intent(out)           :: d2angdthph(1:25)
 real*8,intent(in)            :: costh, cosph, sinph, cosmph(0:4), sinmph(0:4)
 integer                      :: l
 real*8                       ::sinth(0:4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calking calculates the orientation dependence of the integral
! = RY_{l,m} [theta, phi]
!
! Speed increases could be made by precomputing powers of sin/ cos th/phi etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do l=0,4
    d2angdthph(L**2+L+1) = 0d0;
  end do
 
    d2angdthph(2)  = -4.886025119029199d-1*cosph*costh
    d2angdthph(4)  = 4.886025119029199d-1*costh*sinph
 
    d2angdthph(5)  = 2.1850968611841581411d0*costh*cosmph(2)*sinth(1)
    d2angdthph(6)  = -1.0925484305920792d0*cosph*costh**2+1.0925484305920792d0*cosph*sinth(2)
    d2angdthph(8)  = 1.0925484305920792d0*costh**2*sinph-1.0925484305920792d0*sinph*sinth(2)
    d2angdthph(9)  = -2.1850968611841581411d0*costh*sinmph(2)*sinth(1)
 
    d2angdthph(10)  = -5.3103923093397915931d0*costh*cosmph(3)*sinth(2)
    d2angdthph(11)  = -7.226528606601386d-1*cosph**2*sinth(1)+6.503875745941246d0*cosph**2*costh**2*sinth(1)+  &
7.226528606601386d-1*sinph**2*sinth(1)-6.503875745941246d0*costh**2*sinph**2*sinth(1)-  &
2.1679585819804155d0*cosph**2*sinth(3)+2.1679585819804155d0*sinph**2*sinth(3)
    d2angdthph(12)  = -1.1426144986611644d-1*cosph*costh-1.7139217479917466d0*cosph*costh**3+  &
5.1417652439752395d0*cosph*costh*sinth(2)
    d2angdthph(14)  = 1.1426144986611644d-1*costh*sinph+1.7139217479917466d0*costh**3*sinph-  &
5.1417652439752395d0*costh*sinph*sinth(2)
    d2angdthph(15)  = 1.445305721320277d0*cosph*sinph*sinth(1)-1.3007751491882493d1*cosph*costh**2*sinph*sinth(1)+  &
4.335917163960831d0*cosph*sinph*sinth(3)
    d2angdthph(16)  = 5.3103923093397915931d0*costh*sinmph(3)*sinth(2)
 
    d2angdthph(17)  = 10.013371767186818153d0*costh*sinth(3)*cosmph(4)
    d2angdthph(18)  = -2.655196154669896d0*cosph**3*costh**2+2.655196154669896d0*cosph**3*costh**4+  &
7.965588464009688d0*cosph*costh**2*sinph**2-7.965588464009688d0*cosph*costh**4*sinph**2+  &
2.655196154669896d0*cosph**3*sinth(2)-1.5931176928019375d1*cosph**3*costh**2*sinth(2)-  &
7.965588464009688d0*cosph*sinph**2*sinth(2)+4.779353078405812d1*cosph*costh**2*sinph**2*sinth(2)+  &
2.655196154669896d0*cosph**3*sinth(4)-7.965588464009688d0*cosph*sinph**2*sinth(4)
    d2angdthph(19)  = -1.8923493915151202d0*cosph**2*costh*sinth(1)+1.324644574060584d1*cosph**2*costh**3*sinth(1)+  &
1.8923493915151202d0*costh*sinph**2*sinth(1)-1.324644574060584d1*costh**3*sinph**2*sinth(1)-  &
1.324644574060584d1*cosph**2*costh*sinth(3)+1.324644574060584d1*costh*sinph**2*sinth(3)
    d2angdthph(20)  = -3.345232717786446d-1*cosph*costh**2-2.341662902450512d0*cosph*costh**4+  &
3.345232717786446d-1*cosph*sinth(2)+1.4049977414703072d1*cosph*costh**2*sinth(2)-  &
2.341662902450512d0*cosph*sinth(4)
    d2angdthph(22)  = 3.345232717786446d-1*costh**2*sinph+2.341662902450512d0*costh**4*sinph-  &
3.345232717786446d-1*sinph*sinth(2)-1.4049977414703072d1*costh**2*sinph*sinth(2)+  &
2.341662902450512d0*sinph*sinth(4)
    d2angdthph(23)  = 3.7846987830302403d0*cosph*costh*sinph*sinth(1)-2.649289148121168d1*cosph*costh**3*sinph*sinth(1)+  &
2.649289148121168d1*cosph*costh*sinph*sinth(3)
    d2angdthph(24)  = 7.965588464009688d0*cosph**2*costh**2*sinph-7.965588464009688d0*cosph**2*costh**4*sinph-  &
2.655196154669896d0*costh**2*sinph**3+2.655196154669896d0*costh**4*sinph**3-  &
7.965588464009688d0*cosph**2*sinph*sinth(2)+4.779353078405812d1*cosph**2*costh**2*sinph*sinth(2)+  &
2.655196154669896d0*sinph**3*sinth(2)-1.5931176928019375d1*costh**2*sinph**3*sinth(2)-  &
7.965588464009688d0*cosph**2*sinph*sinth(4)+2.655196154669896d0*sinph**3*sinth(4)
    d2angdthph(25)  = -10.013371767186818153d0*costh*sinth(3)*sinmph(4)


return
end subroutine calcd2angdthph

!***************************************************************


