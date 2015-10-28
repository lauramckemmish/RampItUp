subroutine QAxyzeroang(ang,QAz,maxL)
 implicit none
 integer, intent(in)       :: maxL
 real*8, intent(in)        :: QAz
 real*8, intent(out)       :: ang(1:(1+maxL)**2)
 real*8                    :: pf
 integer                   :: L

 ang = 0d0;

 pf = 0.28209479177387814347d0;
 do L = 0,maxL
    ang(L**2+L+1) = sign(1d0,QAz)**(L)*sqrt(2d0*L+1d0)*pf
 end do;

return
end subroutine QAxyzeroang

!****************************************************************************************************

subroutine QAxyzerodang(dangdx,dangdy,dangdz,QAz,maxL)
 implicit none
 integer, intent(in)       :: maxL
 real*8, intent(in)        :: QAz
 real*8, intent(out)       :: dangdx(1:(1+maxL)**2), dangdy(1:(1+maxL)**2), dangdz(1:(1+maxL)**2)

 dangdx = 0d0; dangdy = 0d0; dangdz = 0d0;  

  if (maxL == 0) then; return;
  end if
  dangdx(4) =  -0.48860251190291992159d0/abs(QAz)
  dangdy(2) =  -0.48860251190291992159d0/abs(QAz)

  if (maxL == 1) then; return;
  end if
  dangdx(8) =  -1.0925484305920790705d0/QAz
  dangdy(6) =  -1.0925484305920790705d0/QAz

  if (maxL == 2) then; return;
  end if
  dangdx(14) =  -1.8281831978578629446d0/abs(QAz)
  dangdy(12) =  -1.8281831978578629446d0/abs(QAz)

  if (maxL == 3) then; return;
  end if

  dangdx(22) =  -2.6761861742291566718d0/QAz
  dangdy(20) =  -2.6761861742291566718d0/QAz

  if (maxL == 4) then; return;
  end if

  dangdx(32) =  -3.6235732095655753704d0/abs(QAz)
  dangdy(30) =  -3.6235732095655753704d0/abs(QAz)

  if (maxL == 5) then; return;
  end if

  dangdx(44) = -4.6609709001498511107d0/QAz
  dangdy(42) = -4.6609709001498511107d0/QAz

  if (maxL == 6) then; return;
  end if

  dangdx(58) = -5.7812228852811081108d0/abs(QAz)
  dangdy(56) = -5.7812228852811081108d0/abs(QAz)

  if (maxL == 7) then; return;
  end if
 
  if (maxL > 7) then; 
     print *, "This isn't coded in QAxyzerodang"; call exit
   end if


return
end subroutine QAxyzerodang



!****************************************************************************************************

subroutine QAxyzerod2ang(d2angdx2,d2angdy2,d2angdz2,d2angdxy,d2angdxz,d2angdyz,QAz,maxL)
 implicit none
 integer, intent(in)       :: maxL
 real*8, intent(in)        :: QAz
 real*8, intent(out)       :: d2angdx2(1:(1+maxL)**2),d2angdy2(1:(1+maxL)**2),d2angdz2(1:(1+maxL)**2)
 real*8, intent(out)       :: d2angdxy(1:(1+maxL)**2),d2angdxz(1:(1+maxL)**2),d2angdyz(1:(1+maxL)**2)

 d2angdx2 = 0d0; d2angdy2 = 0d0; d2angdz2 = 0d0; d2angdxy = 0d0; d2angdxz = 0d0; d2angdyz = 0d0;
 
  if (maxL == 0) then; return;
  end if
  
  d2angdx2(3) = -0.48860251190291992159d0/(abs(QAz)*QAz)
  d2angdy2(3) = -0.48860251190291992159d0/(abs(QAz)*QAz)
 
  d2angdxz(4) = 0.48860251190291992159d0/(abs(QAz)*QAz)
  d2angdyz(2) = 0.48860251190291992159d0/(abs(QAz)*QAz)

  if (maxL == 1) then; return;
  end if
  d2angdx2(7) = -1.8923493915151200362d0/(QAz**2)
  d2angdx2(9) =  1.0925484305920790705d0/(QAz**2)
  d2angdy2(7) = -1.8923493915151200362d0/(QAz**2)
  d2angdy2(9) = -1.0925484305920790705d0/(QAz**2)
  
  d2angdxy(5) =  1.0925484305920790705d0/(QAz**2);
  d2angdxz(8) =  1.0925484305920790705d0/(QAz**2);
  d2angdyz(6) =  1.0925484305920790705d0/(QAz**2);

  if (maxL == 2) then; return;
  end if

  d2angdx2(13) = -4.4781159910813846970d0/(abs(QAz)*QAz)
  d2angdx2(15) =  2.8906114426405540554d0/(abs(QAz)*QAz)
  d2angdy2(13) = -4.4781159910813846970d0/(abs(QAz)*QAz)
  d2angdy2(15) = -2.8906114426405540554d0/(abs(QAz)*QAz)

  d2angdxy(11) =  2.8906114426405540554d0/(abs(QAz)*QAz)
  d2angdxz(14) =  1.8281831978578629446d0/(abs(QAz)*QAz)
  d2angdyz(12) =  1.8281831978578629446d0/(abs(QAz)*QAz)

  if (maxL == 3) then; return;
  end if
  d2angdx2(21) = -8.4628437532163443042d0/(QAz**2)
  d2angdx2(23) =  5.6770481745453601086d0/(QAz**2)
  d2angdy2(21) = -8.4628437532163443042d0/(QAz**2)
  d2angdy2(23) = -5.6770481745453601086d0/(QAz**2)
  
  d2angdxy(19) =  5.6770481745453601086d0/(QAz**2)
  d2angdxz(22) =  2.6761861742291566718d0/(QAz**2)
  d2angdyz(20) =  2.6761861742291566718d0/(QAz**2)

  if (maxL == 4) then; return;
  end if

  d2angdx2(31) = -14.034038694410831573d0/(abs(QAz)*QAz)
  d2angdx2(33) =  9.5870735699466475100d0/(abs(QAz)*QAz)
  d2angdy2(31) = -14.034038694410831573d0/(abs(QAz)*QAz)
  d2angdy2(33) = -9.5870735699466475100d0/(abs(QAz)*QAz)

  d2angdxy(29) =  9.5870735699466475100d0/(abs(QAz)*QAz)
  d2angdxz(32) =  3.6235732095655753704d0/(abs(QAz)*QAz)
  d2angdyz(30) =  3.6235732095655753704d0/(abs(QAz)*QAz)

  if (maxL == 5) then; return;
  end if

  d2angdx2(43) = -21.359251961923151113d0/(QAz**2)
  d2angdx2(45) =  14.739284152238775986d0/(QAz**2)
  d2angdy2(43) = -21.359251961923151113d0/(QAz**2)
  d2angdy2(45) = -14.739284152238775986d0/(QAz**2)

  d2angdxy(41) =  14.739284152238775986d0/(QAz**2)
  d2angdxz(44) =  4.6609709001498511107d0/(QAz**2)
  d2angdyz(42) =  4.6609709001498511107d0/(QAz**2)

  if (maxL == 6) then; return;
  end if

  d2angdx2(57) = -30.591356056578213975d0/(abs(QAz)*QAz)
  d2angdx2(59) =  21.241569237359166372d0/(abs(QAz)*QAz)
  d2angdy2(57) = -30.591356056578213975d0/(abs(QAz)*QAz)
  d2angdy2(59) = -21.241569237359166372d0/(abs(QAz)*QAz)

  d2angdxy(55) =  21.241569237359166372d0/(abs(QAz)*QAz)
  d2angdxz(58) =  5.7812228852811081108d0/(abs(QAz)*QAz) 
  d2angdyz(56) =  5.7812228852811081108d0/(abs(QAz)*QAz) 

  if (maxL == 7) then; return;
  end if
 
  if (maxL > 7) then; 
     print *, "This isn't coded in QAxyzerod2ang"; call exit
   end if


return
end subroutine QAxyzerod2ang