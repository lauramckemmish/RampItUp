subroutine calcggmain(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggss,ggsp,ggpp,noggss,noggsp,noggpp,nobasisfun,&
               basis,moleculename,lengthtwoelec,newtwo,noatoms)
 implicit none
  real*8, intent(in)      :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
  real*8, intent(in)      :: basis(1:nobasisfun,1:30)
  character(len=30), intent(in)           :: moleculename
  real*8, intent(in)      :: Ptot(1:nobasisfun,1:nobasisfun)
  real*8, intent(inout)   :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
  integer                 :: noggss, noggsp,noggpp,nobasisfun, noatoms
  real*8, intent(in)      :: ggss(1:20,1:noggss), ggsp(1:20,1:noggsp), ggpp(1:20,1:noggpp)
  integer                 :: lengthtwoelec, i, j, k, l
  integer     :: shortnoggss, shortnoggsp, shortnoggpp
  integer     :: stss, stsp, stpp, finss, finsp, finpp
  real*8 :: newtwo(0:lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,0: lengthtwoelec,1:1)


  shortnoggss = ceiling(dble(noggss)/dble(noatoms))
  shortnoggsp = ceiling(dble(noggsp)/dble(noatoms))
  shortnoggpp = ceiling(dble(noggpp)/dble(noatoms))


  call calcggggssss(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggss,noggss,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
  call calcggggspsp(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggsp,noggsp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
  call calcggggpppp(Galpha,Gbeta,Palpha,Pbeta,Ptot,ggpp,noggpp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)


finss = 0;   finsp = 0;   finpp = 0;
 do i=1,noatoms
  stss = finss;          finss = min(stss+shortnoggss,noggss);
  finsp = 0;
do j=1,noatoms
   stsp = finsp;          finsp = min(stsp+shortnoggsp,noggsp);
  call calcggggsssp(Galpha,Gbeta,Palpha,Pbeta,Ptot, &
    ggss(1:20,stss+1:finss),ggsp(1:20,stsp+1:finsp),finss-stss,finsp-stsp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
  end do
 end do



finss = 0;   finsp = 0;   finpp = 0;
 do i=1,noatoms
  stsp = finsp;          finsp = min(stsp+shortnoggsp,noggsp);
  finpp = 0;
do j=1,noatoms
   stpp = finpp;          finpp = min(stpp+shortnoggpp,noggpp);


  call calcggggsppp(Galpha,Gbeta,Palpha,Pbeta,Ptot, &
    ggsp(1:20,stsp+1:finsp),ggpp(1:20,stpp+1:finpp),finsp-stsp,finpp-stpp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
  end do
 end do


finss = 0;   finsp = 0;   finpp = 0;
 do i=1,noatoms
  stss = finss;          finss = min(stss+shortnoggss,noggss);

  finpp = 0
  do j=1,noatoms 
   stpp = finpp;          finpp = min(stpp+shortnoggpp,noggpp);

   call calcggggsspp(Galpha,Gbeta,Palpha,Pbeta,Ptot, &
     ggss(1:20,stss+1:finss),ggpp(1:20,stpp+1:finpp),finss-stss,finpp-stpp,nobasisfun,basis,moleculename, &
           newtwo,lengthtwoelec)
  end do
 end do


  end subroutine calcggmain
