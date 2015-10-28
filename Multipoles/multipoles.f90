subroutine multipoles(flatdigmp,calcmp, L,nLpure,noRR, noRgSs, noRgSp, RR, RgSs,RgSp, nobasisfun,noatoms,&
          atoms, mpprint,kmax, algchange, moleculename,pts,type)
  implicit none
  integer, intent(in)                   :: L, nLpure, noatoms, nobasisfun , kmax, algchange(1:30)
  integer, intent(in)                   :: noRR, noRgSs, noRgSp, mpprint,pts,type
  real*8, intent(in)   :: atoms(1:5,1:noatoms),RR(1:30,1:noRR), RgSs(1:30,1:noRgSs), RgSp(1:30,1:noRgSp)
  character(len=30), intent(in)         :: moleculename
  integer, intent(out)                  :: calcmp(1:nobasisfun,1:nobasisfun,1:noatoms)
  real*8, intent(out)                   :: flatdigmp(1:nobasisfun*(1+nobasisfun)/2,1:nLpure+2)
  real*8         :: mpRgSs(1:noRgSs,1:nLpure+4), mpRgSp(1:3*noRgSp,1:nLpure+4), mpRRSS(1:noRR,1:nLpure+4)
  integer        :: i, j, cnt, nLCart


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This code controls the calculations of all multipole moments and the digestion of these
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  calcmp = 0;  flatdigmp = 0d0; nLCart =  (L+1)*(L+2)*(L+3)/6

  call calcmpRRSS(mpRRSS,  noRR, RR,noatoms,atoms, nLpure)

  call calcgenmpRgSs(mpRgSs(1:noRgSs,1:nLpure+4),noRgSs,RgSs(1:30,1:noRgSs),&
       L, nLCart,nLpure,kmax) 
  
  mpRgSs(1:noRgSs,nLPure+3) =  mpRgSs(1:noRgSs,nLPure+3)/(2d0*sqrt(dacos(-1d0)));

!DELETE ALL Ss shell pairs formed from two non-concentric Ss shells 
  do i=1,noRgSs
    if (int(RgSs(14,i)) == 111) then;
print *, "IMPORTANT -", int(RgSs(1:2,i)), "multipoles clobbered!! - charge originally",  mpRgSs(i,1)
        mpRgSs(i,1:nLpure) = 0d0;
    end if
  end do

  call calcgenmpRgSp(mpRgSp(1:3*noRGSp,1:nLpure+4),noRgSp,RgSp(1:30,1:noRGSp),&
      L, nLCart,nLpure,kmax)
  mpRgSp(1:3*noRgSp,nLPure+3) =  mpRgSp(1:3*noRgSp,nLPure+3)/(2d0*sqrt(dacos(-1d0)));

  call mpdigestion(flatdigmp,calcmp, mpRRSS ,mpRgSs,mpRgSp,nobasisfun,noRR, noRgSs,noRgSp,nLpure, noatoms)


!DELETE ALL Ss shell pairs formed from two non-concentric Ss shells 
  do i=1,noRgSs
    if (int(RgSs(14,i)) == 111) then;
       calcmp(RgSs(1,i),RgSs(2,i),1:noatoms) = 0;
       calcmp(RgSs(2,i),RgSs(1,i),1:noatoms) = 0;
    end if
  end do

  if (mpprint >=1) then
     !!************ Print out Multipole moments ******************
     print *, "These are digested monopole moments"; cnt = 0;
     do i=1,nobasisfun;   do j=i,nobasisfun;   cnt = cnt + 1;
!        if (flatdigmp(cnt,nLpure+2) == 1) then
           print "(a1,i3,a1,i3,a3,f20.15)", "(", i, ",", j, ") = ",flatdigmp(cnt,1)
!        end if;  
     end do; end do
     cnt = 0; do i=1,nobasisfun;   do j=i,nobasisfun;      cnt = cnt + 1;
!     if (flatdigmp(cnt,nLpure+2) == 1) then
        print "(a1,i3,a1,i3,a3,f20.15,f20.15,f20.15)", "(", i, ",", j, ") = ",flatdigmp(cnt,2:4)
!     end if;   
  end do; end do;
end if

if (mpprint > 1) then; print *, "These are digested (2,m)-pole moments"; cnt = 0;
 do i=1,nobasisfun;   do j=i,nobasisfun;      cnt = cnt + 1;
      if (flatdigmp(cnt,nLpure+2) == 1) then
        print "(a1,i3,a1,i3,a3,f20.15,f20.15,f20.15,f20.15,f20.15)", "(", i, ",", j, ") = ",flatdigmp(cnt,5:9)
      end if;   end do; end do;end if

if (mpprint > 2 .and. L > 2) then; print *, "These are digested (3,m)-pole moments"; cnt = 0;
 do i=1,nobasisfun;   do j=i,nobasisfun;      cnt = cnt + 1;
      if (flatdigmp(cnt,nLpure+2) == 1) then
        print "(a1,i3,a1,i3,a3,7f20.15)", "(", i, ",", j, ") = ",flatdigmp(cnt,10:16)
      end if;   end do; end do; end if

if (mpprint > 3 .and. L > 3) then; print *, "These are digested (4,m)-pole moments"; cnt = 0;
 do i=1,nobasisfun;   do j=i,nobasisfun;      cnt = cnt + 1;
      if (flatdigmp(cnt,nLpure+2) == 1) then
        print "(a1,i3,a1,i3,a3,9f15.10)", "(", i, ",", j, ") = ",flatdigmp(cnt,17:25)
      end if;   end do; end do; end if


if (mpprint > 4 .and. L > 4) then; print *, "These are digested (5,m)-pole moments * 1d10"; cnt = 0;
 do i=1,nobasisfun;   do j=i,nobasisfun;      cnt = cnt + 1;
      if (flatdigmp(cnt,nLpure+2) == 1) then
        print "(a1,i3,a1,i3,a3,11f15.10)", "(", i, ",", j, ") = ",1d10*flatdigmp(cnt,26:36)
      end if;   end do; end do; end if

if (mpprint > 5 .and. L > 5) then; print *, "These are digested (6,m)-pole moments * 1d10"; cnt = 0;
 do i=1,nobasisfun;   do j=i,nobasisfun;      cnt = cnt + 1;
      if (flatdigmp(cnt,nLpure+2) == 1) then
        print "(a1,i3,a1,i3,a3,13f15.10)", "(", i, ",", j, ") = ",1d10*flatdigmp(cnt,37:49)
      end if;   end do; end do; end if

if (mpprint > 6 .and. L > 6) then; print *, "These are digested (7,m)-pole moments * 1d10"; cnt = 0;
 do i=1,nobasisfun;   do j=i,nobasisfun;      cnt = cnt + 1;
      if (flatdigmp(cnt,nLpure+2) == 1) then
        print "(a1,i3,a1,i3,a3,15f15.10)", "(", i, ",", j, ") = ",1d10*flatdigmp(cnt,50:64)
      end if;   end do; end do; end if

!print *, "finish multipoles"
return
end subroutine multipoles
