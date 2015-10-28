


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine readQchemoverlap(qchemol,nobasisfun,basis, olreadprint,totolread,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun, olreadprint
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(out)                     :: qchemol(1:nobasisfun**2,1:4)
  integer, intent(out)                    :: totolread
  real*8                                  :: x, alpha, nc, pi
  integer                                 :: i, index2, bf, isinput(1:nobasisfun,1:nobasisfun)
  character(len=30)                       :: moleculename
  character(len=40)                       :: filename
  
  qchemol = 0d0;
  totolread = 0;
  
  isinput = 0;
  filename = 'Data/'//trim(moleculename)
  filename = trim(filename) // "_ol_int.csv"
  open(unit=15, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  
  do i=1,nobasisfun**2
     
     ! This if statement ensures we only read unique ket values. 
     if (index2 > 0) then
        if (isinput(int(qchemol(index2,2)), int(qchemol(index2,3))) == 0) then
           isinput(int(qchemol(index2,2)),int( qchemol(index2,3))) = 1;
           isinput(int(qchemol(index2,3)),int( qchemol(index2,2))) = 1;
           index2 = index2+1
        end if
     else 
        index2 = index2 + 1
     end if
     read(15, *, end=200) qchemol(index2,2)
     read(15, *, end=200) qchemol(index2,3)
     read(15, *, end=200) qchemol(index2,1)
     qchemol(index2,4) = 1
     bf = int(qchemol(index2,2))
     if (basis(bf,3) ==5) then
        alpha = basis(bf,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        qchemol(index2,4) = qchemol(index2,4)*nc*basis(bf,22)
     elseif (basis(bf,3) ==0) then
        qchemol(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 6 .and. basis(bf,3) .le. 8) then;
        qchemol(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 9 .and. basis(bf,3) .le. 11) then;
        alpha = basis(bf,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        qchemol(index2,4) = qchemol(index2,4)*nc*basis(bf,22)
     end if
     bf = int(qchemol(index2,3))
     if (basis(bf,3) ==5) then
        alpha = basis(bf,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        qchemol(index2,4) = qchemol(index2,4)*nc*basis(bf,22)
     elseif (basis(bf,3) ==0) then
        qchemol(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 6 .and. basis(bf,3) .le. 8) then;
        qchemol(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 9 .and. basis(bf,3) .le. 11) then;
        alpha = basis(bf,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        qchemol(index2,4) = qchemol(index2,4)*nc*basis(bf,22)
     end if
     
  end do
200 totolread = index2
  close(15)
  
  if (olreadprint ==1) then
     print *, "OL read print integrals"
     do i=1,totolread
        print "(f5.0,f5.0,f20.15,f20.15)", qchemol(i,2:3), qchemol(i,1), qchemol(i,4) 
     end do
  end if
  return
end subroutine readQchemoverlap


!*************************************************

subroutine readQchemkinetic(qchemke,nobasisfun,basis, kereadprint,totkeread,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun, kereadprint
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(out)                     :: qchemke(1:nobasisfun**2,1:4)
  integer, intent(out)                    :: totkeread
  integer                                 :: i, index2, bf, isinput(1:nobasisfun,1:nobasisfun)
  real*8                                  :: x, alpha, nc, pi
  character(len=30)                       :: moleculename
  character(len=40)                       :: filename
  
  qchemke = 0d0;
  isinput = 0;
  
  filename = 'Data/'//trim(moleculename)
  filename = trim(filename) // "_ke_int.csv"
  open(unit=16, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  do i=1,nobasisfun**2 
     if (index2 > 0) then
        if (isinput(int(qchemke(index2,2)), int(qchemke(index2,3))) == 0) then
           isinput(int(qchemke(index2,2)),int( qchemke(index2,3))) = 1;
           isinput(int(qchemke(index2,3)),int( qchemke(index2,2))) = 1;
           index2 = index2+1
        end if
     else 
        index2 = index2 + 1
     end if
     read(16,*, end=250) qchemke(index2,2)
     read(16,*,end=250) qchemke(index2,3)
     read(16,*,end=250) qchemke(index2,1)
     qchemke(index2,4) = 1
     bf = int(qchemke(index2,2))
     if (basis(bf,3) ==5) then
        alpha = basis(bf,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        qchemke(index2,4) = qchemke(index2,4)*nc*basis(bf,22)
     elseif (basis(bf,3) ==0) then
        qchemke(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 6 .and. basis(bf,3) .le. 8) then;
        qchemke(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 9 .and. basis(bf,3) .le. 11) then;
        alpha = basis(bf,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        qchemke(index2,4) = qchemke(index2,4)*nc*basis(bf,22)
     end if
     bf = int(qchemke(index2,3))
     if (basis(bf,3) ==5) then
        alpha = basis(bf,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        qchemke(index2,4) = qchemke(index2,4)*nc*basis(bf,22)
     elseif (basis(bf,3) ==0) then
        qchemke(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 6 .and. basis(bf,3) .le. 8) then;
        qchemke(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 9 .and. basis(bf,3) .le. 11) then;
        alpha = basis(bf,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        qchemke(index2,4) = qchemke(index2,4)*nc*basis(bf,22)
     end if
  end do
250 totkeread = index2
  
  close(15)
  
  if (kereadprint ==1) then
     print *, "KE read print integrals"
     do i=1,totkeread
        print "(f5.0,f5.0,f20.15,f20.15)", qchemke(i,2:3), qchemke(i,1), qchemke(i,4) 
     end do
  end if
  return
end subroutine readQchemkinetic





!**********************************************************
subroutine readQchemna(qchemna,nobasisfun,basis, noatoms, nareadprint,totnaread,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun, noatoms, nareadprint
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(out)                     :: qchemna(1:nobasisfun**2,1:4)
  integer, intent(out)                    :: totnaread
  real*8                                  :: x, alpha, nc, pi
  integer                                 :: i, index2, bf, atid, isinput(0:nobasisfun,0:nobasisfun)
  character(len=30)                       :: moleculename
  character(len=40)                       :: filename
  print *, "starting"
  qchemna = 0d0;
  
  isinput = 0;
  
  filename = 'Data/'//trim(moleculename)
  filename = trim(filename) // "_na_int.csv"
  open(unit=15, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  do i=1,nobasisfun**2
     if (index2 > 0) then
        if (isinput( int( qchemna(index2,2) ), int( qchemna(index2,3) ) ) == 0) then
           isinput(int(qchemna(index2,2)),int( qchemna(index2,3))) = 1;
           isinput(int(qchemna(index2,3)),int( qchemna(index2,2))) = 1;
           index2 = index2+1
        end if
     else 
        index2 = index2 + 1
     end if
     read(15, *, end=200) qchemna(index2,2)
     read(15, *, end=200) qchemna(index2,3)
     read(15, *, end=200) qchemna(index2,1)
     qchemna(index2,4) = 1d0
     bf = int(qchemna(index2,2))
     if (basis(bf,3) ==5) then
        alpha = basis(bf,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        qchemna(index2,4) = qchemna(index2,4)*nc*basis(bf,22)
     elseif (basis(bf,3) ==0) then
        qchemna(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 6 .and. basis(bf,3) .le. 8) then;
        qchemna(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 9 .and. basis(bf,3) .le. 11) then;
        alpha = basis(bf,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        qchemna(index2,4) = qchemna(index2,4)*nc*basis(bf,22)
     end if
     bf = int(qchemna(index2,3))
     if (basis(bf,3) ==5) then
        alpha = basis(bf,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        qchemna(index2,4) = qchemna(index2,4)*nc*basis(bf,22)
     elseif (basis(bf,3) ==0) then
        qchemna(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 6 .and. basis(bf,3) .le. 8) then;
        qchemna(index2,4) = 0d0;
     elseif (basis(bf,3) .ge. 9 .and. basis(bf,3) .le. 11) then;
        alpha = basis(bf,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        qchemna(index2,4) = qchemna(index2,4)*nc*basis(bf,22)
     end if

  end do
200 totnaread = i
  close(15)

  
  if (nareadprint ==1) then
     print *, "Read NA integrals from Qchem"
     do i=1,nobasisfun**2
        print "(f5.0,f5.0,f20.15,f20.15)", qchemna(i,2:4), qchemna(i,1)   
     end do
  end if
  print *, "finished"
  return
end subroutine readQchemna


!*****************************************************
subroutine readQchemtwoee(Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,basis,moleculename,compareprint)
  implicit none
  integer, intent(in)                     :: nobasisfun,compareprint
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(in)                      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)                      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)                   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8                                  :: normtwoelecmat2(0:nobasisfun,0: nobasisfun,0: nobasisfun,0: nobasisfun,1:2)
  integer                                 :: i, index2, bf1, bf2, bf3, bf4, m, n, l, s
  real*8                                  :: x, alpha, nc, pi, intval
  character(len=30), intent(in)           :: moleculename
  character(len=100)                      :: filename

  normtwoelecmat2 = 0d0;
  
  filename = "Data/"
  filename = trim(filename)//trim(moleculename)
  filename = trim(filename) // "_twoelec_int.csv"
  open(unit=15, file = filename)
  
  pi = dacos(-1d0);
  index2 = 0;
  do i=1,nobasisfun**4
     index2 = index2+1
     read(15, *, end=200) bf1
     read(15, *, end=200) bf2
     read(15, *, end=200) bf3
     read(15, *, end=200) bf4
     read(15, *, end=200) intval

    if (basis(bf1,3) == 5) then; 
        alpha = basis(bf1,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf1,22)
    else if (basis(bf1,3) == 0) then;
        intval = 0d0;
    else if (basis(bf1,3) .ge. 6 .and. basis(bf1,3) .le. 8) then;
        intval = 0d0;
     else if (basis(bf1,3) .ge. 9 .and. basis(bf1,3) .le. 11) then;
        alpha = basis(bf1,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf1,22)
    end if

    if (basis(bf2,3) == 5) then; 
        alpha = basis(bf2,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf2,22)
    else if (basis(bf2,3) == 0) then;
        intval = 0d0;
    else if (basis(bf2,3) .ge. 6 .and. basis(bf2,3) .le. 8) then;
        intval = 0d0;
     else if (basis(bf2,3) .ge. 9 .and. basis(bf2,3) .le. 11) then;
        alpha = basis(bf2,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf2,22)
    end if

    if (basis(bf3,3) == 5) then; 
        alpha = basis(bf3,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf3,22)
    else if (basis(bf3,3) == 0) then;
        intval = 0d0;
    else if (basis(bf3,3) .ge. 6 .and. basis(bf3,3) .le. 8) then;
        intval = 0d0;
     else if (basis(bf3,3) .ge. 9 .and. basis(bf3,3) .le. 11) then;
        alpha = basis(bf3,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf3,22)
    end if

    if (basis(bf4,3) == 5) then; 
        alpha = basis(bf4,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf4,22)
    else if (basis(bf4,3) == 0) then;
        intval = 0d0;
    else if (basis(bf4,3) .ge. 6 .and. basis(bf4,3) .le. 8) then;
        intval = 0d0;
     else if (basis(bf4,3) .ge. 9 .and. basis(bf4,3) .le. 11) then;
        alpha = basis(bf4,12)
        nc = 1/(2d0*alpha**(5d0/4d0)*(2d0/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf4,22)
    end if

     normtwoelecmat2(bf1,bf2,bf3,bf4,1) = intval;
     normtwoelecmat2(bf1,bf2,bf4,bf3,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf2,bf1,bf3,bf4,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf2,bf1,bf4,bf3,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)  
     
     normtwoelecmat2(bf3,bf4,bf1,bf2,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf3,bf4,bf2,bf1,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf4,bf3,bf1,bf2,1) = normtwoelecmat2(bf3,bf4,bf2,bf1,1)  
     normtwoelecmat2(bf4,bf3,bf2,bf1,1) = normtwoelecmat2(bf3,bf4,bf1,bf2,1) 
  end do
200 print *, ""; 
  close(15)

  do m = 1,nobasisfun;  do n = 1,nobasisfun;  do l = 1,nobasisfun;   do s = 1,nobasisfun
     Galpha(m,n) = Galpha(m,n) + Ptot(l,s)* normtwoelecmat2(m,n,l,s,1) - Palpha(l,s)* normtwoelecmat2(m,l,s,n,1)
     Gbeta(m,n)  = Gbeta(m,n)  + Ptot(l,s)* normtwoelecmat2(m,n,l,s,1) -  Pbeta(l,s)* normtwoelecmat2(m,l,s,n,1)
  end do;  end do; end do;end do

  return
end subroutine readQchemtwoee


!*****************************************************
subroutine readQchemtwoee2(normtwoelecmat,Galpha,Gbeta,Palpha,Pbeta,Ptot,nobasisfun,basis,moleculename, &
   lengthtwoelec,compareprint)
  implicit none
  integer, intent(in)                     :: nobasisfun, lengthtwoelec,compareprint
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(in)                      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)                      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)                   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  real*8                                  :: normtwoelecmat(0:lengthtwoelec,0:lengthtwoelec,0:lengthtwoelec,0:lengthtwoelec,1:2)
  integer                                 :: i, index2, bf1, bf2, bf3, bf4, j, brano, ketno, braketno
  real*8                                  :: x, alpha, nc, pi, intval, var2, y
  character(len=30), intent(in)           :: moleculename
  character(len=100)                      :: filename

call exit
Galpha = 0d0; Gbeta = 0d0;
  filename = "Data/"
  filename = trim(filename)//trim(moleculename)
  filename = trim(filename) // "_twoelec_int.csv"
  open(unit=15, file = filename)
  
  pi = dacos(-1d0);
  index2 = 0;
  do i=1,nobasisfun**4
     index2 = index2+1
     read(15, *, end=200) bf1
     read(15, *, end=200) bf2
     read(15, *, end=200) bf3
     read(15, *, end=200) bf4
     read(15, *, end=200) intval

    if (basis(bf1,3) == 5) then; 
        alpha = basis(bf1,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf1,22)
    else if (basis(bf1,3) == 0 .or. (basis(bf1,3) >=6 .and. basis(bf1,3)<= 8)) then;
        intval = 0d0;
    end if

    if (basis(bf2,3) == 5) then; 
        alpha = basis(bf2,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf2,22)
    else if (basis(bf2,3) == 0 .or. (basis(bf2,3) >=6 .and. basis(bf2,3)<= 8)) then;
        intval = 0d0;
    end if

    if (basis(bf3,3) == 5) then; 
        alpha = basis(bf3,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf3,22)
    else if (basis(bf3,3) == 0 .or. (basis(bf3,3) >=6 .and. basis(bf3,3)<= 8)) then;
        intval = 0d0;
    end if

    if (basis(bf4,3) == 5) then; 
        alpha = basis(bf4,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf4,22)
    else if (basis(bf4,3) == 0 .or. (basis(bf4,3) >=6 .and. basis(bf4,3)<= 8)) then;
        intval = 0d0;
    end if

  if (bf1 == bf2) then; brano = 1; else; brano = 2; end if
  if (bf3 == bf4) then; ketno = 1; else; ketno = 2; end if
 
  if ((bf1 == bf3 .and. bf2 == bf4) .or. (bf1 == bf4 .and. bf2 == bf3)) then; braketno = 1;
  else; braketno = 2; end if;

     y =ketno*Ptot(bf3,bf4)* intval 
     Galpha(bf1,bf2) = Galpha(bf1,bf2) + y
     Gbeta(bf1,bf2)  = Gbeta(bf1,bf2)  + y
     Galpha(bf2,bf1) = Galpha(bf1,bf2)
     Gbeta(bf2,bf1)  = Gbeta(bf1,bf2)

     if (braketno > 1) then;
        x = brano*Ptot(bf1,bf2)*intval;
        Galpha(bf3,bf4) = Galpha(bf3,bf4) + x
        Gbeta(bf3,bf4)  = Gbeta(bf3,bf4)  + x
        Galpha(bf4,bf3) = Galpha(bf3,bf4)
        Gbeta(bf4,bf3)  = Gbeta(bf3,bf4)
     end if
     
     Galpha(bf1,bf4) = Galpha(bf1,bf4) - Palpha(bf2,bf3)* intval 
     Gbeta(bf1,bf4)  = Gbeta(bf1,bf4)  - Pbeta(bf2,bf3)* intval
     if (ketno > 1) then;
        Galpha(bf1,bf3) = Galpha(bf1,bf3) -Palpha(bf2,bf4)*intval
        Gbeta(bf1,bf3)  = Gbeta(bf1,bf3)  - Pbeta(bf2,bf4)*intval
     end if
     if (brano > 1) then;
        Galpha(bf2,bf4) = Galpha(bf2,bf4) -Palpha(bf1,bf3)*intval
        Gbeta(bf2,bf4)  = Gbeta(bf2,bf4)  - Pbeta(bf1,bf3)*intval
        if (ketno > 1) then;
           Galpha(bf2,bf3) = Galpha(bf2,bf3) -Palpha(bf1,bf4)*intval
           Gbeta(bf2,bf3)  = Gbeta(bf2,bf3)  - Pbeta(bf1,bf4)*intval
        end if
     end if

   if (braketno > 1) then;  
     Galpha(bf4,bf1) = Galpha(bf4,bf1) - Palpha(bf3,bf2)*intval 
     Gbeta(bf4,bf1)  = Gbeta(bf4,bf1)  - Pbeta(bf3,bf2)* intval
     if (ketno > 1) then;
        Galpha(bf3,bf1) = Galpha(bf3,bf1) -Palpha(bf4,bf2)*intval
        Gbeta(bf3,bf1)  = Gbeta(bf3,bf1)  - Pbeta(bf4,bf2)*intval
     end if
     if (brano > 1) then;
        Galpha(bf4,bf2) = Galpha(bf4,bf2) -Palpha(bf3,bf1)*intval
        Gbeta(bf4,bf2)  = Gbeta(bf4,bf2)  - Pbeta(bf3,bf1)*intval
        if (ketno > 1) then;
           Galpha(bf3,bf2) = Galpha(bf3,bf2) -Palpha(bf4,bf1)*intval
           Gbeta(bf3,bf2)  = Gbeta(bf3,bf2)  - Pbeta(bf4,bf1)*intval
        end if
     end if
   end if

  end do
200 print *, ""; 
  close(15)


do i=1,nobasisfun; do j=1,nobasisfun
  if (Palpha(i,j) .ne. Palpha(j,i)) then
   print *, i, j, "dig P differs", Palpha(i,j), Palpha(j,i)
  end if
  if (Pbeta(i,j) .ne. Pbeta(j,i)) then
   print *, i, j, "dig2", Pbeta(i,j), Pbeta(j,i)
  end if
end do; end do;

do i=1,nobasisfun; do j=1,nobasisfun
  if (Galpha(i,j) .ne. Galpha(j,i)) then
   print *, i, j, "dig", Galpha(i,j), Galpha(j,i)
  end if
  if (Gbeta(i,j) .ne. Gbeta(j,i)) then
   print *, i, j, "dig2", Gbeta(i,j), Gbeta(j,i)
  end if
end do; end do;

if (compareprint == 1) then;

  normtwoelecmat = 0d0;
  
  filename = "Data/"
  filename = trim(filename)//trim(moleculename)
  filename = trim(filename) // "_twoelec_int.csv"
  open(unit=15, file = filename)
  
  pi = dacos(-1d0);
  index2 = 0;
  do i=1,nobasisfun**4
     index2 = index2+1
     read(15, *, end=300) bf1
     read(15, *, end=300) bf2
     read(15, *, end=300) bf3
     read(15, *, end=300) bf4
     read(15, *, end=300) intval

    if (basis(bf1,3) == 5) then; 
        alpha = basis(bf1,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf1,22)
    else if (basis(bf1,3) == 0 .or. (basis(bf1,3) >=6 .and. basis(bf1,3)<= 8)) then;
        intval = 0d0;
    end if

    if (basis(bf2,3) == 5) then; 
        alpha = basis(bf2,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf2,22)
    else if (basis(bf2,3) == 0 .or. (basis(bf2,3) >=6 .and. basis(bf2,3)<= 8)) then;
        intval = 0d0;
    end if

    if (basis(bf3,3) == 5) then; 
        alpha = basis(bf3,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf3,22)
    else if (basis(bf3,3) == 0 .or. (basis(bf3,3) >=6 .and. basis(bf3,3)<= 8)) then;
        intval = 0d0;
    end if

    if (basis(bf4,3) == 5) then; 
        alpha = basis(bf4,12)
        nc = 1/((2d0*alpha/pi)**(3d0/4d0))
        intval = intval*nc*basis(bf4,22)
    else if (basis(bf4,3) == 0 .or. (basis(bf4,3) >=6 .and. basis(bf4,3)<= 8)) then;
        intval = 0d0;
    end if

     normtwoelecmat(bf1,bf2,bf3,bf4,1) = intval;
     normtwoelecmat(bf1,bf2,bf4,bf3,1) = normtwoelecmat(bf1,bf2,bf3,bf4,1)
     normtwoelecmat(bf2,bf1,bf3,bf4,1) = normtwoelecmat(bf1,bf2,bf3,bf4,1)
     normtwoelecmat(bf2,bf1,bf4,bf3,1) = normtwoelecmat(bf1,bf2,bf3,bf4,1)  
     
     normtwoelecmat(bf3,bf4,bf1,bf2,1) = normtwoelecmat(bf1,bf2,bf3,bf4,1)
     normtwoelecmat(bf3,bf4,bf2,bf1,1) = normtwoelecmat(bf1,bf2,bf3,bf4,1)
     normtwoelecmat(bf4,bf3,bf1,bf2,1) = normtwoelecmat(bf3,bf4,bf2,bf1,1)  
     normtwoelecmat(bf4,bf3,bf2,bf1,1) = normtwoelecmat(bf3,bf4,bf1,bf2,1) 
  end do
300 print *, ""; 
  close(15)
end if

  return
end subroutine readQchemtwoee2





!**********************************************************
subroutine readQchempseudona(qchemna,nobasisfun,totnaread,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(out)                     :: qchemna(1:nobasisfun**2,1:4)
  integer                                 :: i, index2, bf, atid, isinput(0:nobasisfun,0:nobasisfun)
  real*8                                  :: x, alpha, nc, pi
  integer, intent(out)                    :: totnaread
  character(len=30)                       :: moleculename
  character(len=40)                       :: filename
  
  qchemna = 0d0;
  
  isinput = 0;
  
  filename = 'Data/'//trim(moleculename)
  filename = trim(filename) // "_pseudona_int.csv"
  open(unit=15, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  
  do i=1,nobasisfun**2
     if (index2 > 0) then
        if (isinput( int( qchemna(index2,2) ), int( qchemna(index2,3) ) ) == 0) then
           isinput(int(qchemna(index2,2)),int( qchemna(index2,3))) = 1;
           isinput(int(qchemna(index2,3)),int( qchemna(index2,2))) = 1;
           index2 = index2+1
        end if
     else 
        index2 = index2 + 1
     end if
     read(15, *, end=200) qchemna(index2,2)
     read(15, *, end=200) qchemna(index2,3)
     read(15, *, end=200) qchemna(index2,1)
     qchemna(index2,4) = 1
  end do
200 totnaread = i-1
  close(15)
  
  
  return
end subroutine readQchempseudona



!**********************************************************
subroutine readQchempseudomp(qchemmp,nobasisfun,totmpread,filename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(out)                     :: qchemmp(1:nobasisfun**2,1:4)
  integer                                 :: i, index2, bf, atid, isinput(0:nobasisfun,0:nobasisfun)
  real*8                                  :: x, alpha, nc, pi
  integer, intent(out)                    :: totmpread
  character(len=40)                       :: filename
  
  qchemmp = 0d0;
  
  isinput = 0;
  
  open(unit=15, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  
  do i=1,nobasisfun**2
     if (index2 > 0) then
        if (isinput( int( qchemmp(index2,2) ), int( qchemmp(index2,3) ) ) == 0) then
           isinput(int(qchemmp(index2,2)),int( qchemmp(index2,3))) = 1;
           isinput(int(qchemmp(index2,3)),int( qchemmp(index2,2))) = 1;
           index2 = index2+1
        end if
     else 
        index2 = index2 + 1
     end if
     read(15, *, end=200) qchemmp(index2,2)
     read(15, *, end=200) qchemmp(index2,3)
     read(15, *, end=200) qchemmp(index2,1)
     qchemmp(index2,4) = 1
  end do
200 totmpread = i-1
  close(15)
  
  
  return
end subroutine readQchempseudomp





!*****************************************************
subroutine readQchemtwoeecompare(normtwoelecmat2,nobasisfun,basis,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8                                  :: normtwoelecmat2(0:nobasisfun,0: nobasisfun,0: nobasisfun,0: nobasisfun,1:1)
  integer                                 :: i, index2, bf1, bf2, bf3, bf4, m, n, l, s
  real*8                                  :: x, alpha, nc, pi, intval
  character(len=30), intent(in)           :: moleculename
  character(len=100)                      :: filename

  normtwoelecmat2 = 0d0;
  
print *, "moleculename = ", moleculename
  filename = "Data/"
  filename = trim(filename)//trim(moleculename)
  filename = trim(filename) // "_twoelec_int.csv"
  open(unit=15, file = filename)

  
  pi = dacos(-1d0);
  index2 = 0;
  do i=1,nobasisfun**4
     index2 = index2+1
     read(15, *, end=200) bf1
     read(15, *, end=200) bf2
     read(15, *, end=200) bf3
     read(15, *, end=200) bf4
     read(15, *, end=200) intval

     normtwoelecmat2(bf1,bf2,bf3,bf4,1) = intval;
     normtwoelecmat2(bf1,bf2,bf4,bf3,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf2,bf1,bf3,bf4,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf2,bf1,bf4,bf3,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)  
     
     normtwoelecmat2(bf3,bf4,bf1,bf2,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf3,bf4,bf2,bf1,1) = normtwoelecmat2(bf1,bf2,bf3,bf4,1)
     normtwoelecmat2(bf4,bf3,bf1,bf2,1) = normtwoelecmat2(bf3,bf4,bf2,bf1,1)  
     normtwoelecmat2(bf4,bf3,bf2,bf1,1) = normtwoelecmat2(bf3,bf4,bf1,bf2,1) 
  end do
200 print *, ""; 
  close(15)

  return
end subroutine readQchemtwoeecompare























!*************************



subroutine writeoverlap(normolmat,nobasisfun,basis,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(in)                     :: normolmat(0:nobasisfun,0:nobasisfun,1:2)
  real*8                                  :: x, alpha, nc, pi
  integer                                 :: i, index2, bf, isinput(1:nobasisfun,1:nobasisfun), bf1, bf2
  character(len=30)                       :: moleculename
  character(len=40)                       :: filename, filename2
  
  isinput = 0;
  filename = 'Data/'//trim(moleculename)
  filename = trim(filename) // "_ol_int.csv"
  open(unit=15, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  
  filename2 = "Integrals/"
!  filename2 = trim(filename2)//trim(moleculename)
  filename2 = trim(filename2) // "ol.txt"
  open(unit=20, file = filename2);

  do bf1=1,nobasisfun
   do bf2=1,nobasisfun
     write(20,"(i3,a3,i3,a3,f30.15)") int(bf1), ",", int(bf2),",", normolmat(bf1,bf2,1)
   end do 
  end do
200 close(15)
  
  return
end subroutine writeoverlap



!*************************



subroutine writekinetic(normkemat,nobasisfun,basis, moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(in)                     :: normkemat(0:nobasisfun,0:nobasisfun,1:2)
  real*8                                  :: x, alpha, nc, pi
  integer                                 :: i, index2, bf, isinput(1:nobasisfun,1:nobasisfun), bf1, bf2
  character(len=30)                       :: moleculename
  character(len=40)                       :: filename, filename2
  
  
  isinput = 0;
  filename = 'Data/'//trim(moleculename)
  filename = trim(filename) // "_ol_int.csv"
  open(unit=15, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  
  filename2 = "Integrals/"
!  filename2 = trim(filename2)//trim(moleculename)
  filename2 = trim(filename2) // "ke.txt"
  open(unit=20, file = filename2);

  do bf1=1,nobasisfun
   do bf2=1,nobasisfun
     write(20,"(i3,a3,i3,a3,f30.15)") int(bf1), ",", int(bf2),",", normkemat(bf1,bf2,1)
   end do 
  end do

200 close(15)
  
  return
end subroutine writekinetic


!*************************



subroutine writena(normnamat,nobasisfun,basis,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(in)                     :: normnamat(0:nobasisfun,0:nobasisfun,1:2)
  real*8                                  :: x, alpha, nc, pi
  integer                                 :: i, index2, bf, isinput(1:nobasisfun,1:nobasisfun), bf1, bf2
  character(len=30)                       :: moleculename
  character(len=40)                       :: filename, filename2
 
  
  isinput = 0;
  filename = 'Data/'//trim(moleculename)
  filename = trim(filename) // "_ol_int.csv"
  open(unit=15, file = filename)
  pi = dacos(-1d0);
  index2 = 0;
  
  filename2 = "Integrals/"
  filename2 = trim(filename2) // "na.txt"
  open(unit=20, file = filename2);

  do bf1=1,nobasisfun
   do bf2=1,nobasisfun
     write(20,"(i3,a3,i3,a3,f30.15)") int(bf1), ",", int(bf2),",", normnamat(bf1,bf2,1)
   end do 
  end do

200 close(15)
  
  return
end subroutine writena








!*****************************************************
subroutine writetwoee(newtwo,nobasisfun,basis,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(in)                      :: basis(1:nobasisfun,1:30)
  real*8, intent(in)                      :: newtwo(0:nobasisfun,0:nobasisfun,0:nobasisfun,0:nobasisfun,1:1)
  integer                                 :: i, index2, bf1, bf2, bf3, bf4, firstbf1, firstbf2, firstbf3, firstbf4
  real*8                                  :: x, alpha, nc, pi, intval
  character(len=30), intent(in)           :: moleculename
  character(len=100)                      :: filename, filename2
   Character( Len = 10 ) :: ca



do bf1=1,nobasisfun
if (int(basis(bf1,3)) == 0 .or. int(basis(bf1,3)) == 5) then
do bf2 = 1,nobasisfun
  filename2 = "Integrals/"
  filename2 = trim(filename2) // "ramptwoelec_int_"
  if (bf1 .lt. 10) then    
    write(ca,'(i1)') bf1
  else
    write(ca,'(i2)') bf1
  end if
  filename2 = trim(filename2) // trim(ca) 
  filename2 = trim(filename2) // "_"
  if (bf2 .lt. 10) then
    write(ca,'(i1)') bf2
  else
    write(ca,'(i2)') bf2
  end if
  filename2 = trim(filename2) // trim(ca)
  filename2 = trim(filename2) // ".txt"
  open(unit=20, file = filename2);
  do bf3 = 1,nobasisfun
   do bf4 = 1,nobasisfun
    if (abs(newtwo(bf1,bf2,bf3,bf4,1)) .gt. 1d-10) then 
     write(20,"(i8,i8,i8,i8,f30.15)") int(bf1), int(bf2), int(bf3), int(bf4), newtwo(bf1,bf2,bf3,bf4,1)
    end if
   end do
  end do
close(20)
end do
end if
end do



  return
end subroutine writetwoee














!*****************************************************
subroutine readQchempseudotwoee(qchemee,nobasisfun,toteeread,moleculename)
  implicit none
  integer, intent(in)                     :: nobasisfun
  real*8, intent(out)                     :: qchemee(1:6,1:nobasisfun**4)
  integer                                 :: i, index2, bf
  real*8                                  :: x, alpha, nc, pi
  integer, intent(out)                    :: toteeread
  character(len=30), intent(in)           :: moleculename
  character(len=100)                      :: filename
  
  qchemee = 0d0;        toteeread = 0;
  
  filename = "Data/";     filename = trim(filename)//trim(moleculename)
  filename = trim(filename) // "_pseudotwoelec_int.csv";
  open(unit=15, file = filename);
  
  index2 = 0;
  do i=1,nobasisfun**4;
     index2 = index2+1;
     read(15, *, end=200) qchemee(2,index2);
     read(15, *, end=200) qchemee(3,index2);
     read(15, *, end=200) qchemee(4,index2);
     read(15, *, end=200) qchemee(5,index2);
     read(15, *, end=200) qchemee(1,index2);
     qchemee(6,index2) = 1;
  end do;
200 toteeread = i-1;
  close(15);
  return;
end subroutine readQchempseudotwoee





