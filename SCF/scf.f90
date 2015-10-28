subroutine scf(shortol,coreH,nobasisfun,noatoms, atoms,noelec, multiplicity, &
invintatomdist, coreHprint, consXprint,moleculename,kmax, flatdigmp,Tintmat,&
totnomodel, nLpure, algchange, noggss, noggsp, noggpp,  &
ggss,ggsp, ggpp, basis,RRRRpreftable,L, maxrampdegree,Tcutoff,quadpts, compareprint, &
newsortRg,lenmodelnewsort,maxrgpa,norampsperatom, totnoramps ,rampinfo)

implicit none
integer, intent(in)           :: norampsperatom(1:noatoms), totnoramps, rampinfo(1:2,1:totnoramps)
integer, intent(in)           :: nobasisfun, noatoms, noelec, kmax, multiplicity, coreHprint, consXprint
real*8, intent(in)            :: atoms(1:5,1:noatoms), invintatomdist(1:noatoms,1:noatoms)
real*8, intent(in)            :: shortol(1:nobasisfun,1:nobasisfun), coreH(1:nobasisfun,1:nobasisfun)
integer                       :: i, run, j, k,m, stat, maxcycles, alphaoccupied, betaoccupied,n
integer                       :: compareprint,maxrgpa
real*8                        :: coreHprime(1:nobasisfun,1:nobasisfun), Stemp(1:nobasisfun,1:nobasisfun)
real*8                        :: X(1:nobasisfun,1:nobasisfun), oneelecenergy, nucrepenergy, tolerance
real*8                        :: coalphaprime(1:nobasisfun,1:nobasisfun), cobetaprime(1:nobasisfun,1:nobasisfun)
real*8                        :: coalpha(1:nobasisfun,1:nobasisfun), cobeta(1:nobasisfun,1:nobasisfun)
real*8                        :: Palpha(1:nobasisfun,1:nobasisfun),Pbeta(1:nobasisfun,1:nobasisfun)
real*8                        :: Ptot(1:nobasisfun,1:nobasisfun)
real*8                        :: Palphaold(1:nobasisfun,1:nobasisfun), Pbetaold(1:nobasisfun,1:nobasisfun)
real*8                        :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
real*8                        :: Fockalpha(1:nobasisfun,1:nobasisfun), Fockbeta(1:nobasisfun,1:nobasisfun)
real*8                        :: totenergy, totpastenergy, minerror, Hartrampenergy, Hartpseudoenergy, Hartgausenergy
real*8                        :: GausAE, RampAE
character(len=30)             :: moleculename, basisset
character(len=80)             :: filename
character(len=2)              :: kmaxstr
character(len=5)              :: betachangealgstr
integer, intent(in)           :: L, maxrampdegree, nLpure, algchange(1:20), quadpts, noggss, noggsp, noggpp, totnomodel
real*8, intent(in)            :: RRRRpreftable(0:maxrampdegree,0:maxrampdegree,0:L), Tcutoff
real*8, intent(in)            :: basis(1:nobasisfun,1:30), flatdigmp(1:nobasisfun*(nobasisfun+1)/2,1:nLpure+2)
real*8, intent(in)            :: ggss(1:20,1:noggss), ggsp(1:20,1:noggsp), ggpp(1:20,1:noggpp)
real*8                        :: oldggpp(1:20,1:noggpp)
real*8, intent(in)            :: Tintmat(1:nLpure,1:nLpure,0:noatoms,0:noatoms), newsortRg(1:4,1:totnomodel)
integer, intent(in)           :: lenmodelnewsort(1:4,1:maxrgpa,1:noatoms)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  This code takes input of the integrals (overlap, kinetic, n.a. and two electron) and performs a SCF calculation
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "DEBUG.h"
Palphaold = 0d0;
Pbetaold = 0d0;
#ifdef INTPRINT
do i=1,1
Palpha = 0d0;
Pbeta = 0d0;
Ptot = 0d0;
oldggpp = ggpp
call twoelec(Galpha,Gbeta,Palpha,Pbeta,Ptot,flatdigmp,Tintmat,&
totnomodel, nLpure, algchange, noggss, noggsp, noggpp,  &
ggss,ggsp, oldggpp, nobasisfun, basis, noatoms, atoms,&
moleculename, RRRRpreftable,L, maxrampdegree, Tcutoff,quadpts, compareprint,&
newsortRg,lenmodelnewsort,maxrgpa,norampsperatom, totnoramps ,rampinfo)
end do
#endif
#ifdef WRITETWOEE
do i=1,1
Palpha = 0d0;
Pbeta = 0d0;
Ptot = 0d0;
oldggpp = ggpp
call twoelec(Galpha,Gbeta,Palpha,Pbeta,Ptot,flatdigmp,Tintmat,&
totnomodel, nLpure, algchange, noggss, noggsp, noggpp,  &
ggss,ggsp, oldggpp, nobasisfun, basis, noatoms, atoms,&
moleculename, RRRRpreftable,L, maxrampdegree, Tcutoff,quadpts, compareprint,&
newsortRg,lenmodelnewsort,maxrgpa,norampsperatom, totnoramps ,rampinfo)
end do
#endif
#ifdef TIMINGRUN
do i=1,1
Palpha = 0d0;
Pbeta = 0d0;
Ptot = 0d0;
oldggpp = ggpp
call twoelec(Galpha,Gbeta,Palpha,Pbeta,Ptot,flatdigmp,Tintmat,&
totnomodel, nLpure, algchange, noggss, noggsp, noggpp,  &
ggss,ggsp, oldggpp, nobasisfun, basis, noatoms, atoms,&
moleculename, RRRRpreftable,L, maxrampdegree, Tcutoff,quadpts, compareprint,&
newsortRg,lenmodelnewsort,maxrgpa,norampsperatom, totnoramps ,rampinfo)
end do
#endif


#ifdef TIMINGRUN
call exit
#endif
!****** Calculates the orthogonalisation matrix *****
call calcX(shortol, X, nobasisfun, consXprint)

!*** Calculates the core Hamiltonian ****
call calccore(oneelecenergy, coalphaprime,cobetaprime,coreHprime, coreH,X, nobasisfun,coreHprint)

!**** Determines number of alpha and beta electrons in line with number of electrons and multiplicity of molecule.
alphaoccupied = ceiling(dble(noelec)/2d0); betaoccupied = floor(dble(noelec)/2d0);
if (multiplicity == 3 .or. multiplicity ==4) then;     alphaoccupied = alphaoccupied + 1;   betaoccupied = betaoccupied - 1
else if (multiplicity > 3) then;     print *, "Please change scf.f90 since this multiplicity can't be accommodated";   call exit;
end if
!**** Calculates the density matrix guess from the core Hamiltonian ***
call consdensmat(Palpha,Pbeta,Ptot,coalphaprime,cobetaprime, X, alphaoccupied, betaoccupied, nobasisfun)

totenergy = oneelecenergy;
totpastenergy = oneelecenergy + 1000d0;
run = 1;
tolerance = 10d0;
if (noelec == 1) then;
maxcycles = 0;
else if (compareprint == 1) then
maxcycles = 2;
else
maxcycles = 400;
end if
do while (abs(totenergy - totpastenergy) .gt. 10**(-tolerance) .and. run .lt. maxcycles)
run = run + 1;
if (run == 30 .or. &
run == 100 .or. run==125 .or. run == 150 .or. (run > 300 .and. run < 360)) then
Palpha = 0.5*(Palpha+Palphaold);
Pbeta = 0.5*(Pbeta+Pbetaold);
Ptot = Palpha + Pbeta;
end if
Palphaold = Palpha;
Pbetaold = Pbeta;
totpastenergy = totenergy;
oldggpp = ggpp;

call twoelec(Galpha,Gbeta,Palpha,Pbeta,Ptot,flatdigmp,Tintmat,&
totnomodel, nLpure, algchange, noggss, noggsp, noggpp,  &
ggss,ggsp, oldggpp, nobasisfun, basis, noatoms, atoms,&
moleculename, RRRRpreftable,L, maxrampdegree, Tcutoff,quadpts, compareprint,&
newsortRg,lenmodelnewsort,maxrgpa,norampsperatom, totnoramps ,rampinfo)

call Fock(Fockalpha,Fockbeta,coalphaprime,cobetaprime,Galpha,Gbeta,coreH,X,nobasisfun)
call consdensmat(Palpha,Pbeta,Ptot,coalphaprime,cobetaprime, X, alphaoccupied, betaoccupied, nobasisfun)
call calctotenergy(totenergy,Palpha,Pbeta,Ptot,Fockalpha,Fockbeta,coreH, nobasisfun)
print *, "run no = ", run, "   tot energy = ", totenergy



#ifdef GPRT
do i=1,nobasisfun
do j=1,nobasisfun
if (abs(Palpha(i,j)) .gt. 100) then
print *, i, j, Palpha(i,j), Pbeta(i,j), Ptot(i,j)
end if
end do
end do
#endif

#ifdef GPRT
do i=1,nobasisfun
do j=1,nobasisfun
if (abs(Galpha(i,j)) .gt. 100) then
print *, i, j, Galpha(i,j), Gbeta(i,j)
end if
end do
end do
#endif
end do

if (run == maxcycles) then;
print *, totpastenergy, totenergy;
print *, "**************** ERROR *****************";  print *, "Not converged in ", maxcycles , "cycles";
end if;

print "(i4,a15)", run, " scf cycles"

!*** calculates the nuclear-nuclear energy contribution
call calcnucrepenergy(nucrepenergy,noatoms, atoms, invintatomdist)


Hartrampenergy = totenergy + nucrepenergy;
print *, "Hartrampenergy = ", Hartrampenergy
filename = 'Data/'//trim(moleculename); filename = trim(filename) // 'qchemenergies.out'
open(unit=18,file = filename);   read(18,*) Hartgausenergy;
close(18)
print "(a20,f25.15)", "R-31+G energy is :", Hartrampenergy


print "(a20,f20.15,f20.6)", "6-31+G energy is : ", Hartgausenergy ,Hartgausenergy;
print "(a40,f20.15,f10.6,a15,f10.5,a10)", "Ramp energy - All-Gaussian energy= ",&
(-Hartgausenergy + Hartrampenergy),  (-Hartgausenergy + Hartrampenergy), " mHartree = ", &
2625.49962d0*(-Hartgausenergy + Hartrampenergy), "kJ/mol";print *, "";

print "(a20,f20.15,f20.5)", "6-31+G energy is : ", Hartgausenergy ,Hartgausenergy;
print "(a40,f20.15,f10.5,a15,f10.5,a10)", "Ramp energy - All-Gaussian energy= ",&
(-Hartgausenergy + Hartrampenergy),  (-Hartgausenergy + Hartrampenergy), " mHartree = ", &
2625.49962d0*(-Hartgausenergy + Hartrampenergy), "kJ/mol";print *, "";


open(unit = 16, file = "Data/basisset.csv");      read(16,*) basisset;     close(16)

if (basisset == 'R31pG') then;
GausAE = Hartgausenergy;
RampAE = Hartrampenergy;
do i=1,noatoms
select case (int(atoms(2,i)))
case(1)
GausAE = GausAE + 0.498232910700000d0;
RampAE = RampAE + 0.498232910700000d0;
case(3)
GausAE = GausAE + 7.431472743300000d0;
RampAE = RampAE + 7.371522197264073d0;
case(4)
GausAE = GausAE + 14.569596292000000d0;
RampAE = RampAE + 14.535052830936646d0
case(5)
GausAE = GausAE + 24.523727778800001d0;
RampAE = RampAE + 24.500210243204869d0;
case(6)
GausAE = GausAE + 37.680924592200000d0;
RampAE = RampAE + 37.661453054861845d0;
case(7)
GausAE = GausAE + 54.386294103300003d0;
RampAE = RampAE + 54.368334644700511d0;
case(8)
GausAE = GausAE + 74.783457490399996d0;
RampAE = RampAE + 74.765877434045066d0;
case(9)
GausAE = GausAE + 99.367521440800004d0;
RampAE = RampAE + 99.350460770374454d0;
case(10)
GausAE = GausAE + 128.483549728700012d0;
RampAE = RampAE + 128.466281575255323d0;
case default
print *, "this element does not have default energy in scf.f90 6-31+G basis";
end select
end do

GausAE = - GausAE;
RampAE = - RampAE;

print "(a20,f20.15,f20.5,a10,f20.5,a10)", "6-31+G atomisation energy is : ", GausAE , GausAE, "Hartree", &
2625.49962d0*GausAE, "kJ/mol"
print "(a20,f20.15,f20.5,a10,f20.5,a10)", "R-31+G atomisation energy is : ", RampAE , RampAE, "Hartree", &
2625.49962d0*GausAE, "kJ/mol"

print "(a40,f20.15,f10.5,a15,f10.5,a10)", "Ramp - All-Gaussian atomisation energy= ",&
(-GausAE + RampAE),  (-GausAE + RampAE), " mHartree = ", &
2625.49962d0*(-GausAE + RampAE), "kJ/mol";print *, "";

end if

return
end subroutine scf

!***************************************************
subroutine consGmat(Galpha, Gbeta, Palpha, Pbeta, Ptot, shorttwoelec,nobasisfun)
implicit none
real*8, intent(in)                        :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
real*8, intent(in)                        :: Ptot(1:nobasisfun,1:nobasisfun)
real*8, intent(in)                        :: shorttwoelec(1:nobasisfun,1:nobasisfun,1:nobasisfun,1:nobasisfun)
real*8, intent(out)                       :: Galpha(1:nobasisfun,1:nobasisfun), Gbeta(1:nobasisfun,1:nobasisfun)
integer                                   :: m,l,s,n, i, j, k
integer, intent(in)                       :: nobasisfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the G matrix (see Szabo and Ostland)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Galpha = 0d0;        Gbeta = 0d0;

do m = 1,nobasisfun;  do n = 1,nobasisfun;  do l = 1,nobasisfun;   do s = 1,nobasisfun
Galpha(m,n) = Galpha(m,n) + Ptot(l,s)* shorttwoelec(m,n,l,s) - Palpha(l,s)* shorttwoelec(m,l,s,n)
Gbeta(m,n)  = Gbeta(m,n)  + Ptot(l,s)* shorttwoelec(m,n,l,s) -  Pbeta(l,s)* shorttwoelec(m,l,s,n)
end do;  end do; end do;end do

return
end subroutine consGmat



!********************************************************
subroutine calccore(oneelecenergy, coalphaprime,cobetaprime,coreHprime, coreH,X, nobasisfun,coreHprint)
implicit none
real*8, intent(out)           :: coreHprime(1:nobasisfun,1:nobasisfun)
real*8, intent(in)            :: X(1:nobasisfun,1:nobasisfun), coreH(1:nobasisfun,1:nobasisfun)
integer, intent(in)           :: coreHprint, nobasisfun
real*8, intent(out)           :: oneelecenergy
real*8, intent(out)           :: coalphaprime(1:nobasisfun,1:nobasisfun), cobetaprime(1:nobasisfun,1:nobasisfun)
real*8                        :: coreHprimetemp(1:nobasisfun,1:nobasisfun), coreHint(1:nobasisfun,1:nobasisfun)
real*8                        :: oneelecenergies(1:nobasisfun), temp(1:200*nobasisfun)
integer                       :: info, i

call dgemm('N','N',nobasisfun,nobasisfun,nobasisfun,1d0, coreH, nobasisfun,X, nobasisfun,0d0,coreHint,nobasisfun)
call dgemm('T','N',nobasisfun,nobasisfun,nobasisfun,1d0, X, nobasisfun,coreHint, nobasisfun,0d0,coreHprime,nobasisfun)

oneelecenergies = 0d0;    coreHprimetemp = coreHprime;
call dsyev('V','U', nobasisfun, coreHprimetemp, nobasisfun, oneelecenergies, temp, 100*nobasisfun, info)

coalphaprime = coreHprimetemp;     cobetaprime = coreHprimetemp

if (coreHprint == 1) then
print *, "Core Hamiltonian"
do i=1,nobasisfun;      print "(13f10.5)", coreH(i,1:13);
end do
print *, "Core Hamiltonian prime";
do i=1,nobasisfun;      print "(13f10.5)", coreHprime(i,1:13);
end do
print *, "Eigenvalues of the Core Hamiltonian prime";    print "(13f10.5)", oneelecenergies;
end if

oneelecenergy = oneelecenergies(1);
return
end subroutine calccore

!********************************************************
subroutine calcX(S, X, nobasisfun,consXprint)
implicit none
integer, intent(in)           :: nobasisfun, consXprint
real*8, intent(in)            :: S(1:nobasisfun,1:nobasisfun)
real*8, intent(out)           :: X(1:nobasisfun,1:nobasisfun)
real*8                        :: Stemp(1:nobasisfun,1:nobasisfun)
integer                       :: i

Stemp = S;   X = 0d0;
call MatPow(X, Stemp, nobasisfun,-1d0/2d0)

if (consXprint ==1) then
print *, "Overlap matrix";    do i=1,nobasisfun;      print "(13f10.5)", S(i,1:13);  end do
print *, "Orthonormalisation matrix"; do i=1,nobasisfun; print "(13f10.5)", X(i,1:13); end do ;
end if
return
end subroutine calcX



!************************************************************

subroutine calcnucrepenergy(nucrepenergy,noatoms, atoms, invintatomdist)
implicit none
integer, intent(in)           :: noatoms
real*8, intent(in)            :: atoms(1:5,1:noatoms),invintatomdist(1:noatoms,1:noatoms)
real*8, intent(out)           :: nucrepenergy
integer                       :: i,j

nucrepenergy = 0d0;
do i=1,noatoms; do j=i+1,noatoms
nucrepenergy =nucrepenergy + atoms(2,i)*atoms(2,j)*invintatomdist(i,j)
end do; end do

return
end subroutine calcnucrepenergy



!*************************************************************************
subroutine consdensmat(Palpha,Pbeta,Ptot,coalphaprime,cobetaprime, X, alphaoccupied, betaoccupied, nobasisfun)
implicit none
integer, intent(in)                       :: nobasisfun
real*8, intent(in)                        :: coalphaprime(1:nobasisfun,1:nobasisfun), cobetaprime(1:nobasisfun,1:nobasisfun)
real*8, intent(inout)                     :: Palpha(1:nobasisfun,1:nobasisfun),Pbeta(1:nobasisfun,1:nobasisfun)
real*8, intent(inout)                     :: Ptot(1:nobasisfun,1:nobasisfun)
integer, intent(in)                       :: alphaoccupied, betaoccupied
real*8                                    :: coalpha(1:nobasisfun,1:nobasisfun), cobeta(1:nobasisfun,1:nobasisfun)
real*8                                    :: X(1:nobasisfun,1:nobasisfun)
integer                                   :: i, k, m

call dgemm('N','N',nobasisfun,nobasisfun,nobasisfun,1d0, X, nobasisfun, coalphaprime, nobasisfun,0d0,coalpha,nobasisfun)
call dgemm('N','N',nobasisfun,nobasisfun,nobasisfun,1d0, X, nobasisfun, cobetaprime, nobasisfun,0d0,cobeta,nobasisfun)

Palpha = 0d0;
Pbeta = 0d0;
do k = 1,nobasisfun; do m=1,nobasisfun
do i=1,alphaoccupied;   Palpha(k,m) = Palpha(k,m) + coalpha(k,i)*coalpha(m,i);  end do
do i=1,betaoccupied;   Pbeta(k,m)  = Pbeta(k,m)  + cobeta(k,i) *cobeta(m,i);
end do
end do;
end do
Ptot = Palpha + Pbeta;

return
end subroutine consdensmat


!*********************************************************
subroutine Fock(Fockalpha,Fockbeta,coalphaprime,cobetaprime,Galpha,Gbeta,coreH,X,nobasisfun)
implicit none
integer, intent(in)                       :: nobasisfun
real*8, intent(in)                        :: Galpha(1:nobasisfun,1:nobasisfun),Gbeta(1:nobasisfun,1:nobasisfun)
real*8, intent(in)                        :: coreH(1:nobasisfun,1:nobasisfun),  X(1:nobasisfun,1:nobasisfun)
real*8, intent(out)                       :: coalphaprime(1:nobasisfun,1:nobasisfun),cobetaprime(1:nobasisfun,1:nobasisfun)
real*8, intent(out)                       :: Fockalpha(1:nobasisfun,1:nobasisfun), Fockbeta(1:nobasisfun,1:nobasisfun)
real*8                                    :: Fockalphaprime(1:nobasisfun,1:nobasisfun), Fockbetaprime(1:nobasisfun,1:nobasisfun)
real*8                                    :: Fockalphaprimetemp(1:nobasisfun,1:nobasisfun)
real*8                                    :: Fockbetaprimetemp(1:nobasisfun,1:nobasisfun)
real*8                                    :: Fockalphaint(1:nobasisfun,1:nobasisfun), Fockbetaint(1:nobasisfun,1:nobasisfun)
real*8                                    :: betaenergies(1:nobasisfun), alphaenergies(1:nobasisfun),temp(1:200*nobasisfun)
integer                                   :: info


Fockalpha = coreH + Galpha
call dgemm('N','N',nobasisfun,nobasisfun,nobasisfun,1d0, Fockalpha, nobasisfun,X, nobasisfun,0d0,Fockalphaint,nobasisfun)
call dgemm('T','N',nobasisfun,nobasisfun,nobasisfun,1d0, X, nobasisfun,Fockalphaint, nobasisfun,0d0,Fockalphaprime,nobasisfun)
Fockalphaprimetemp = Fockalphaprime
call dsyev('V','U', nobasisfun, Fockalphaprimetemp, nobasisfun, alphaenergies, temp, 100*nobasisfun, info)
coalphaprime = Fockalphaprimetemp

Fockbeta = coreH + Gbeta
call dgemm('N','N',nobasisfun,nobasisfun,nobasisfun,1d0, Fockbeta, nobasisfun,X, nobasisfun,0d0,Fockbetaint,nobasisfun)
call dgemm('T','N',nobasisfun,nobasisfun,nobasisfun,1d0, X, nobasisfun,Fockbetaint, nobasisfun,0d0,Fockbetaprime,nobasisfun)
Fockbetaprimetemp = Fockbetaprime
call dsyev('V','U', nobasisfun, Fockbetaprimetemp, nobasisfun, betaenergies, temp, 100*nobasisfun, info)
cobetaprime = Fockbetaprimetemp
return
end subroutine Fock


!********************************
subroutine calctotenergy(totenergy,Palpha,Pbeta,Ptot,Fockalpha,Fockbeta,coreH,nobasisfun)
implicit none
real*8, intent(out)                       :: totenergy
integer, intent(in)                       :: nobasisfun
real*8, intent(in)                        :: Palpha(1:nobasisfun,1:nobasisfun), Pbeta(1:nobasisfun,1:nobasisfun)
real*8, intent(in)                        :: Ptot(1:nobasisfun,1:nobasisfun)
real*8, intent(in)                        :: Fockalpha(1:nobasisfun,1:nobasisfun), Fockbeta(1:nobasisfun,1:nobasisfun)
real*8, intent(in)                        :: coreH(1:nobasisfun,1:nobasisfun)
integer                                   :: i,j

totenergy = 0d0;
do i=1,nobasisfun;  do j=1,nobasisfun
totenergy = totenergy + (1d0/2d0)*(Ptot(i,j)*coreH(i,j)+Palpha(i,j)*Fockalpha(i,j) + Pbeta(i,j)*Fockbeta(i,j))
end do; end do
return
end subroutine calctotenergy



!************************************

SUBROUTINE MatPow(B,A,N,Y)
!     ******************************************************************
!     *                                                                *
!     *         MatPow finds B = A**Y and E = eigenvalues of A         *
!     *         A is a real symmetric positive-definite matrix         *
!     *                                                                *
!     *         PMWG (12/94)                                           *
!     *                                                                *
!     ******************************************************************
implicit none
integer, intent(in)   :: N
real*8, intent(out)   :: B(1:N,1:N)
integer               :: i , info
real*8                :: A(1:N,1:N), E(1:N), Y, temp(100*N)

! input = A;

call dsyev('V','U', N, A, N, E, temp, 100*N, info)
! During this routine, A is initially the input matrix, then becomes its eigvec

! eigvec = A


if (E(1) < 0) then
print *, "Some of the overlap matrix eigenvalues are less than zero --> something wrong has happened"
call exit
end if
do i = 1,N
call dscal(N, E(i)**(Y/2d0),A(1,i), 1)
end do

call dgemm('N','T',N, N, N, 1d0, A, N, A, N, 0d0, B, N)

return

end subroutine MatPow



