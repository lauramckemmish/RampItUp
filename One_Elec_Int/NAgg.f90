subroutine calcggssna(naggss,noggss,ggss,noatoms,atoms)
 implicit none
  integer, intent(in) :: noggss,noatoms
  real*8, intent(in)  :: ggss(1:30,1:noggss),atoms(1:5,1:noatoms)
  real*8              :: naggss(1:noggss,1:4,1:noatoms)
  real*8              :: pi, zeta, alpha, beta, RPC2, Cx, Cy, Cz, f0eval
  integer             :: i, j

 pi = dacos(-1d0);
 do j=1,noatoms
  Cx = atoms(3,j); Cy = atoms(4,j); Cz = atoms(5,j);
  do i=1,noggss
     RPC2 = (ggss(4,i)-Cx)**2+(ggss(5,i)-Cy)**2+(ggss(6,i)-Cz)**2;
     zeta = ggss(8,i)
     call calcf0(f0eval,zeta*RPC2)
     naggss(i,1,j) = 2*pi/zeta*f0eval; 
     naggss(i,2,j) = ggss(1,i)
     naggss(i,3,j) = ggss(2,i)
     naggss(i,4,j) = ggss(9,i)
  end do
end do

!call exit
end subroutine calcggssna

!***************************

subroutine calcggspna(naggsp,noggsp,ggsp,noatoms,atoms)
 implicit none
  integer, intent(in) :: noggsp,noatoms
  real*8, intent(in)  :: ggsp(1:30,1:noggsp),atoms(1:5,1:noatoms)
  real*8              :: naggsp(1:3*noggsp,1:4,1:noatoms)
  real*8              :: pi, zeta, alpha, beta, RPC2, Cx, Cy, Cz, f0eval, pi3r2
  integer             :: i,j
  real*8              :: ABx, ABy, ABz, RAB2, Px, Py, Pz, PCx, PCy, PCz, t1, t2
 
 pi = dacos(-1d0);
 pi3r2 = pi**(3d0/2d0);
 do j=1,noatoms
  Cx = atoms(3,j); Cy = atoms(4,j); Cz = atoms(5,j);
  do i=1,noggsp
     ABx = ggsp(11,i); ABy = ggsp(12,i); ABz = ggsp(13,i);
     RAB2 = ABx**2+ABy**2+ABz**2;
     alpha = ggsp(6,i); beta = ggsp(7,i); zeta = ggsp(8,i);
     Px = ggsp(3,i); Py = ggsp(4,i); Pz = ggsp(5,i);
     PCx = Px-Cx; PCy = Py-Cy; PCz = Pz-Cz;
     RPC2 = (Px-Cx)**2+(Py-Cy)**2+(Pz-Cz)**2;
     call calcf0(f0eval,zeta*RPC2)
     t1 = exp(-RPC2*zeta)*pi/(RPC2*zeta**2);
     t2 = pi*f0eval/(RPC2*zeta**2);
     if (RPC2 .lt. 1d-10) then; t1 = 0d0; t2 = 0d0; end if;
     naggsp(3*i-2    ,1,j) = PCx*(t1-t2) + 2d0*pi*ABx*alpha*f0eval/(zeta**2);
     naggsp(3*i-1    ,1,j) = PCy*(t1-t2) + 2d0*pi*ABy*alpha*f0eval/(zeta**2);
     naggsp(3*i-0    ,1,j) = PCz*(t1-t2) + 2d0*pi*ABz*alpha*f0eval/(zeta**2); 
     naggsp(3*i-2:3*i,2,j) = ggsp(1,i)
     naggsp(3*i-2    ,3,j) = ggsp(2,i)
     naggsp(3*i-1    ,3,j) = ggsp(2,i)+1
     naggsp(3*i      ,3,j) = ggsp(2,i)+2
     naggsp(3*i-2:3*i,4,j) = ggsp(9,i)

  end do
 end do

end subroutine calcggspna


!***************************

subroutine calcggppna(naggpp,noggpp,ggpp,noatoms,atoms)
 implicit none
  integer, intent(in) :: noggpp,noatoms
  real*8, intent(in)  :: ggpp(1:30,1:noggpp),atoms(1:5,1:noatoms)
  real*8              :: naggpp(1:9*noggpp,1:4,1:noatoms)
  real*8              :: pi, zeta, alpha, beta, RPC2, Cx, Cy, Cz, f0eval
  integer             :: i,j
  real*8              :: ABx, ABy, ABz, RAB2, Px, Py, Pz, PCx, PCy, PCz
  real*8              :: t1, t2, t3, ttot, t4, t5, t6, t7

pi = dacos(-1d0);

 do j=1,noatoms
  Cx = atoms(3,j); Cy = atoms(4,j); Cz = atoms(5,j);
  do i=1,noggpp
     ABx = ggpp(11,i); ABy = ggpp(12,i); ABz = ggpp(13,i);
     RAB2 = ABx**2+ABy**2+ABz**2;
     alpha = ggpp(6,i); beta = ggpp(7,i); zeta = ggpp(8,i);
     Px = ggpp(3,i); Py = ggpp(4,i); Pz = ggpp(5,i)
     PCx = Px-Cx; PCy = Py-Cy; PCz = Pz-Cz;
     RPC2 = PCx**2+PCy**2+PCz**2;
 if (RPC2 .lt. 1d-10) then; 
    naggpp(9*i-8:9*i,1,j) = 0d0; 
    naggpp(9*i-8,1,j) = (2*pi-6*pi*ABx**2*alpha*beta/zeta)/(3*zeta**2);
    naggpp(9*i-4,1,j) = (2*pi-6*pi*ABy**2*alpha*beta/zeta)/(3*zeta**2);
    naggpp(9*i-0,1,j) = (2*pi-6*pi*ABz**2*alpha*beta/zeta)/(3*zeta**2);
    naggpp(9*i-7,1,j) = -2*ABx*ABy*pi*alpha*beta/(zeta**3)
    naggpp(9*i-6,1,j) = -2*ABx*ABz*pi*alpha*beta/(zeta**3)
    naggpp(9*i-5,1,j) = -2*ABx*ABy*pi*alpha*beta/(zeta**3)
    naggpp(9*i-3,1,j) = -2*ABy*ABz*pi*alpha*beta/(zeta**3)
    naggpp(9*i-2,1,j) = -2*ABx*ABz*pi*alpha*beta/(zeta**3)
    naggpp(9*i-1,1,j) = -2*ABz*ABy*pi*alpha*beta/(zeta**3)
    

     else

     t1 = -3d0*exp(-RPC2*zeta)*pi/(2*RPC2**2*zeta**3)
     t2 = -exp(-RPC2*zeta)*pi/(RPC2*zeta**2)
     call calcf0(f0eval,zeta*RPC2)
     t3 = 3*pi*f0eval/(2d0*RPC2**2*zeta**3)
     ttot = t1+t2+t3
     t4 = exp(-RPC2*zeta)*pi/(RPC2*zeta**3)
     t5 = -pi*f0eval/(RPC2*zeta**3)
     t6 = pi*f0eval/(zeta**3)

     t7 = exp(-RPC2*zeta)*pi/(2d0*RPC2*zeta**3);
 
     naggpp(9*i-8      ,1,j) = PCx**2*ttot + (t4+t5)*ABx*PCx*(alpha-beta)+t5/2d0&
                                    + (zeta-2*alpha*beta*ABx**2)*t6 + t7 ;
     naggpp(9*i-7      ,1,j) = PCx*PCy*ttot + (t4+t5)*(ABy*PCx*alpha-ABx*PCy*beta)&
                                    -2*alpha*beta*ABx*ABy*t6 ;
     naggpp(9*i-6      ,1,j) = PCx*PCz*ttot + (t4+t5)*(ABz*PCx*alpha-ABx*PCz*beta)&
                                    -2*alpha*beta*ABx*ABz*t6;
     naggpp(9*i-5      ,1,j) = PCx*PCy*ttot + (t4+t5)*(ABx*PCy*alpha-ABy*PCx*beta)&
                                    -2*alpha*beta*ABx*ABy*t6;
     naggpp(9*i-4      ,1,j) = PCy**2*ttot + (t4+t5)*ABy*PCy*(alpha-beta)+t5/2d0&
                                    + (zeta-2*alpha*beta*ABy**2)*t6 + t7;
     naggpp(9*i-3      ,1,j) = PCy*PCz*ttot + (t4+t5)*(ABz*PCy*alpha-ABy*PCz*beta)&
                                    -2*alpha*beta*ABz*ABy*t6;
     naggpp(9*i-2      ,1,j) = PCx*PCz*ttot + (t4+t5)*(ABx*PCz*alpha-ABz*PCx*beta)&
                                    -2*alpha*beta*ABz*ABy*t6;
     naggpp(9*i-1      ,1,j) = PCy*PCz*ttot + (t4+t5)*(ABy*PCz*alpha-ABz*PCy*beta)&
                                    -2*alpha*beta*ABz*ABy*t6;
     naggpp(9*i-0      ,1,j) = PCz**2*ttot + (t4+t5)*ABz*PCz*(alpha-beta)+t5/2d0&
                                    + (zeta-2*alpha*beta*ABz**2)*t6 + t7;

    end if;

if (ggpp(1,i) == ggpp(2,i)) then
   naggpp(9*i-5,1,j) = 0d0;
   naggpp(9*i-2,1,j) = 0d0;
   naggpp(9*i-1,1,j) = 0d0;
end if

     naggpp(9*i-8:9*i-6,2,j) = ggpp(1,i)
     naggpp(9*i-5:9*i-3,2,j) = ggpp(1,i)+1
     naggpp(9*i-2:9*i-0,2,j) = ggpp(1,i)+2
     naggpp(9*i-8      ,3,j) = ggpp(2,i)
     naggpp(9*i-7      ,3,j) = ggpp(2,i)+1
     naggpp(9*i-6      ,3,j) = ggpp(2,i)+2
     naggpp(9*i-5      ,3,j) = ggpp(2,i)
     naggpp(9*i-4      ,3,j) = ggpp(2,i)+1
     naggpp(9*i-3      ,3,j) = ggpp(2,i)+2
     naggpp(9*i-2      ,3,j) = ggpp(2,i)
     naggpp(9*i-1      ,3,j) = ggpp(2,i)+1
     naggpp(9*i-0      ,3,j) = ggpp(2,i)+2
     naggpp(9*i-8:9*i  ,4,j) = ggpp(9,i)


  end do
 end do

end subroutine calcggppna


!*********

subroutine calcf0(f0eval,x)
 implicit none
 real*8, intent(out) :: f0eval
 real*8, intent(in)  :: x

 if (x .lt. 1d-8) then
   f0eval = 1d0-x/3d0;
 else
    f0eval = 0.5d0*sqrt(dacos(-1d0)/x)*erf(sqrt(x))
 end if

end subroutine calcf0
