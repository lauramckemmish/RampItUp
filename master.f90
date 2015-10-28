program master
  implicit none
  character(len=2)                        :: atomtype
  real*8, dimension(:,:), allocatable     :: atoms,basis, basright,tempRR, tempRgSs,tempRgSp, tempggss, tempggsp, tempggpp
  real*8, dimension(:,:), allocatable     :: RR, RgSs, RgSp, ggss, ggsp, ggpp, tempnewsortRg, newsortRg
  real*8, dimension(:,:), allocatable     :: shortol, shortcore,  invintatomdist, flatdigmp, tempmpprim, tempmpcont, mpprim, mpcont
  real*8,dimension(:,:,:), allocatable    :: normolmat, normkemat, combnormnamat, RRRRpreftable
  real*8, dimension(:,:,:,:), allocatable :: Tintmat
  integer                                 :: algchange(1:30), totmpprim, totmpcont, quadpts
  integer                                 :: noRR, noRgSs, noRgSp, noggss, noggsp, noggpp,L, nLpure, nLCart 
  integer                                 :: suppressall,totnomodel, concmodeloption, coreHprint, consXprint, distprint
  integer                                 :: sigspprint,naprint, olprint, keprint, mpprint, Tintprint, modelRgprint
  integer                                 :: twoelecprint, moleculeprint,thresh
  integer                                 :: compareprint, mdlchk, i, j, atid, cnt, a, kmax, maxrampdegree, maxmodellen,lenmdlconc
  integer                                 :: timecnt, clock_start, clock_end, clock_rate, maxrgpa, k, bf3, bf4
  integer                                 :: nobasisfun, noatoms, noelec , multiplicity, toadd, lengthtwoelec, totnoramps
  character(len=150)                      :: filename  
  character(len=30)                       :: moleculename, basisset
  real*8                                  :: t(1:100), Tcutoff
  integer, dimension(:), allocatable      :: norampsperatom
  integer, dimension(:,:),allocatable     :: temprampinfo, rampinfo
  integer, dimension(:,:,:),allocatable   :: lenmodelnewsort, templenmodelnewsort, calcmp

#include "DEBUG.h"

#ifdef DEBUG_PRINT 
 print *, "Debug version running"
#endif




  timecnt = 1;
  call cpu_time(t(timecnt))  
  call SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the 
  call SYSTEM_CLOCK(COUNT=clock_start) ! Start timing
    
  open(unit=20,file = "Data/ramprem.csv");
  read(20,*) kmax;     read(20,*) maxmodellen;    read(20,*) lenmdlconc;  read(20,*) Tcutoff;  read(20,*) quadpts;                
  close(20)
  L=4;   maxmodellen =(lenmdlconc*(1+L)**2);  nLpure = (1+L)**2; nLCart = (L+1)*(L+2)*(L+3)/6;
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !     PRELIMINARIES
  !******** Preliminaries that actually can be pre-computed and put in table or memory but precomputed here for convenience
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !******** Start molecule-independent and basis-function-indepenent preliminaries
  maxrampdegree = 3*lenmdlconc+4;  !Length of model + degree of F ramp

  allocate(RRRRpreftable(0:maxrampdegree,0:maxrampdegree,0:L));  call calcconstants(L, maxrampdegree,RRRRpreftable)
  !******** End molecule-independent and basis-function-indepenent preliminaries
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !    READ MOLECULE INFO FROM FILES
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  open(unit = 16, file = "Data/moleculename.csv");  read(16,*) moleculename; close(16)
  open(unit = 16, file = "Data/basisset.csv");      read(16,*) basisset;     close(16)
  print *, "Calculation for ", moleculename  
  filename = 'Data/geom.csv'
  open(unit=18, file = filename);  read(18,*) noatoms; read(18,*) noelec; read(18,*) multiplicity
  nobasisfun = 0; toadd = 0
  do i=1,noatoms
     read(18,*)          ; read(18,*) atomtype
     if (trim(atomtype) == '3' .or. trim(atomtype) == '4') then
        print *, "************************ ALERT *****************************"
        print *, "Probably can't trust pseudoramp results"
        print *, "************************ ALERT *****************************"
     end if
     read(18,*)  ;  read(18,*)    ;  read(18,*)                ! don't need this line here
     filename = 'Data/'//trim(basisset);  filename = trim(filename)//'/Basis_fortran/'
     filename = trim(filename)//atomtype;  filename = trim(filename)//'.bas'
     open(unit = 25, file = filename);  read(25,*) toadd   ;   close(25)  
     nobasisfun = nobasisfun + toadd
  end do
  close(18)
  
  print "(a30,i4)", "Number of Atoms is ", noatoms;  print "(a30,i4)", "Number of Electrons is ", noelec
  print "(a30,i4)", "Number of Basis functions is", nobasisfun; print "(a30,i4)", "Multiplicity is ", multiplicity
  
  allocate(basis(1:nobasisfun,1:30),basright(1:30,1:nobasisfun), atoms(1:5,1:noatoms))
  allocate(invintatomdist(1:noatoms,1:noatoms))
  
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !     SET PRINT OPTIONS
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************

suppressall = 0
#ifdef DEBUG_PRINT
 suppressall =0
#endif
  
  if (suppressall .ne. 1) then
     moleculeprint     = 0;
     sigspprint        = 0;
#ifdef BASISPRT
     moleculeprint     = 1;                     ! prints out info about the molecule and basis set
     sigspprint        = 1;                     ! prints the number of each kind of significant shell pair
#endif
     distprint         = 0;                     ! prints interatomic distances
     olprint           = 0;                     ! prints the overlap integrals
     keprint           = 0;                     ! prints out kinetic eneRgy integrals
     naprint           = 0;                     ! prints out nuclear attraction integrals
     mpprint           = 0;                     ! If you want to print up to (L)-poles, put in L. 0 prints nothing
     Tintprint         = 0;                     ! prints out the unit-multipole-unit-multipole interaction matrix
     modelRgprint      = 0;
#ifdef MDLPRT
     modelRgprint      = 1;                     ! prints out information about Rg --> sum R modelling
#endif
     coreHprint        = 0;                     ! prints out core Hamiltonian
     consXprint        = 0;                     ! prints out orthogonalisation matrix (i.e. X = S^(-1/2), where S is overlap matrix)
     compareprint      = 0;                     ! compares the pseudoramp with the ramp two-electron integrals 
     mdlchk            = 0;                     ! checks accuracy of Rg models point-by-point.
  end if

  timecnt = timecnt+1;  call cpu_time(t(timecnt))
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !    GET MOLECULE AND BASIS INFO
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  call getmolecule(nobasisfun,noatoms,atoms,basis,moleculeprint,moleculename,kmax,10)
  do i=1,nobasisfun; basright(1:30,i) = basis(i,1:30);  end do
  call getinvinteratomdist(invintatomdist,noatoms,atoms,distprint); print *, "Got molecule information"

  allocate(norampsperatom(1:noatoms))
  allocate(temprampinfo(1:2,1:4*noatoms))
  call getramps(norampsperatom, totnoramps ,temprampinfo,basright,nobasisfun,noatoms)  
  allocate(rampinfo(1:2,1:totnoramps))
  rampinfo(1:2,1:totnoramps) = temprampinfo(1:2,1:totnoramps)
  deallocate(temprampinfo)

  timecnt = timecnt+1;  call cpu_time(t(timecnt));
  print "(a120,f20.12,a5)", "Time for getmolecule is", t(timecnt) -t(timecnt-1), "secs"
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !     GET SHELL PAIR DATA
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************

thresh = 10
  allocate(tempRR(1:30,1:4*nobasisfun),tempRgSs(1:30,1:20*nobasisfun),tempRgSp(1:30,1:20*nobasisfun))

call getsigspRx(algchange,noRR, noRgSs, noRgSp, tempRR, tempRgSs, tempRgSp,nobasisfun,noatoms,atoms,basis,sigspprint, &
thresh)
  print *, "Constructed list of significant shell pairs"
  
  timecnt = timecnt+1;  call cpu_time(t(timecnt))
  print "(a120,f20.12,a5)", "Time for getsigspRx is", t(timecnt) -t(timecnt-1), "secs"
  
  allocate(RR(1:30,1:noRR),RgSs(1:30,1:noRgSs),RgSp(1:30,1:noRgSp))
  RR(1:30,1:noRR) = tempRR(1:30,1:noRR);            RgSs(1:30,1:noRgSs) = tempRgSs(1:30,1:noRgSs)
  RgSp(1:30,1:noRgSp) = tempRgSp(1:30,1:noRgSp);
  deallocate(tempRgSs,tempRgSp, tempRR)

  allocate(tempggss(1:30,1:20*nobasisfun**2))
  call getsigspss(algchange, noggss,tempggss,  nobasisfun,noatoms,atoms,basis,sigspprint,thresh)
  allocate(ggss(1:30,1:noggss));  ggss(1:30,1:noggss) = tempggss(1:30,1:noggss);  deallocate(tempggss)

  allocate(tempggsp(1:30,1:20*nobasisfun**2))
  call getsigspsp(algchange, noggsp,tempggsp,  nobasisfun,noatoms,atoms,basis,sigspprint,thresh)
  allocate(ggsp(1:30,1:noggsp));  ggsp(1:30,1:noggsp) = tempggsp(1:30,1:noggsp);  deallocate(tempggsp)

  allocate(tempggpp(1:30,1:20*nobasisfun**2))
  call getsigsppp(algchange, noggpp,tempggpp,  nobasisfun,noatoms,atoms,basis,sigspprint,thresh)
  allocate(ggpp(1:30,1:noggpp));  ggpp(1:30,1:noggpp) = tempggpp(1:30,1:noggpp);  deallocate(tempggpp)

  timecnt = timecnt+1;  call cpu_time(t(timecnt))
  print "(a120,f20.12,a5)", "Time for getsigspgg is", t(timecnt) -t(timecnt-1), "secs"

  call calcusefulSPinter(algchange,noRR, noRgSs, noRgSp,noggss,noggsp,noggpp,RR, RgSs, RgSp,ggss,ggsp,ggpp)
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !    ONE ELECTRON INTEGRALS
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  allocate(calcmp(1:nobasisfun,1:nobasisfun,1:noatoms), normolmat(0:nobasisfun,0:nobasisfun,1:2))
  allocate(normkemat(0:nobasisfun,0:nobasisfun,1:2), combnormnamat(0:nobasisfun,0:nobasisfun,1:2))

  call overlap(normolmat,algchange,noRR, noRgSs, noRgSp, RR, RgSs, RgSp, nobasisfun, basis, olprint,moleculename,kmax,&
        2*quadpts, ggss, ggsp, ggpp, noggss, noggsp, noggpp)
  timecnt = timecnt+1;  call cpu_time(t(timecnt))
  print "(a120,f20.12,a5)", "Time for overlap is", t(timecnt) -t(timecnt-1), "secs"

  call kinetic(normkemat,algchange,noRR, noRgSs, noRgSp, RR, RgSs, RgSp, nobasisfun,basis, keprint,moleculename,kmax,&
      2*quadpts, ggss, ggsp, ggpp, noggss, noggsp, noggpp)

  print *, "Calculated overlap and kinetic eneRgy matrices"
  timecnt = timecnt+1;  call cpu_time(t(timecnt)); 
 print "(a120,f20.12,a5)", "Time for kinetic is", t(timecnt) -t(timecnt-1), "secs";
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !     MULTIPOLES
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  
  allocate(flatdigmp(1:nobasisfun*(nobasisfun+1)/2,1:nLpure+2))
  call multipoles(flatdigmp,calcmp, L,nLpure,noRR, noRgSs, noRgSp, RR, RgSs,RgSp, nobasisfun,noatoms,&
          atoms, mpprint,kmax, algchange, moleculename,2*quadpts,2)
  print *, "Calculated multipoles"
  
  allocate(Tintmat(1:nLpure,1:nLpure,0:noatoms,0:noatoms))
  Tintmat = 0d0; call calcTintmat(Tintmat,L,nLpure,noatoms,atoms,invintatomdist,Tintprint);
  timecnt = timecnt+1;  call cpu_time(t(timecnt)); 
 print "(a120,f20.12,a5)", "Time for multipoles is", t(timecnt) -t(timecnt-1), "secs"
  
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !     NUCLEAR ATTRACTION
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
 
  call na(combnormnamat,calcmp,algchange,noRR, noRgSs, noRgSp, RR, RgSs, RgSp,nobasisfun,basis,noatoms,atoms, &
       invintatomdist, Tintmat, flatdigmp,L,nLpure, naprint,moleculename,kmax,2*quadpts,&
       ggss, ggsp, ggpp, noggss, noggsp, noggpp)

  print *, "Calculated nuclear attraction matrixes";    deallocate(calcmp)
  
  timecnt = timecnt+1;  call cpu_time(t(timecnt))
  print "(a120,f20.12,a5)", "Time for na ints is", t(timecnt) -t(timecnt-1), "secs"
  
  allocate(shortol(1:nobasisfun,1:nobasisfun),shortcore(1:nobasisfun,1:nobasisfun));          
  shortol   = normolmat(1:nobasisfun,1:nobasisfun,1)
  shortcore = normkemat(1:nobasisfun,1:nobasisfun,1) + combnormnamat(1:nobasisfun,1:nobasisfun,1)
  deallocate(normkemat, normolmat)

  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !   MODEL SHELL PAIRS
  !********* Model Shell Pairs by Sums of Ramps **********************
  !********!!!!!!!!!!!!************!!!!!!!!!!!!************** 
  allocate(tempnewsortRg(1:4,1:nobasisfun*noatoms*maxmodellen),templenmodelnewsort(1:4,1:10*nobasisfun,1:noatoms))
  totnomodel = 0;
  call modelRg(totnomodel, noRR,noRgSs,noRgSp, RR,RgSs,RgSp,nobasisfun,noatoms,maxmodellen,modelRgprint, algchange, &
          maxrampdegree,L,kmax,lenmdlconc,1, rrrrpreftable, L, mdlchk,&
          tempnewsortRg, templenmodelnewsort,maxrgpa)

  allocate(newsortRg(1:4,1:totnomodel));  newsortRg = tempnewsortRg(1:4,1:totnomodel);
  deallocate(tempnewsortRg); print *, "Finished modelRg"

  allocate(lenmodelnewsort(1:4,1:maxrgpa,1:noatoms)); 
  lenmodelnewsort = templenmodelnewsort(1:4,1:maxrgpa,1:noatoms);
  deallocate(templenmodelnewsort); print *, "Finished modelRg"


  !*****************************************
  timecnt = timecnt+1;  call cpu_time(t(timecnt));  
print "(a120,f20.12,a5)", "Time for modelsp is", t(timecnt) -t(timecnt-1), "secs"
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !   TWO ELECTRON INTEGRALS
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  deallocate(RgSs,RgSp)

  timecnt = timecnt+1;  call cpu_time(t(timecnt));  CALL SYSTEM_CLOCK(COUNT=clock_end) ! End timing
  
  print "(a50,f20.10,a10)", "Total cpu_time up to start scf is ", t(timecnt) -t(1), "secs"
  print "(a50,f20.10,a10)", "Total system_clock time up to start scf", real((clock_end-clock_start)/clock_rate), "secs"
  
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  !   SCF CODE
  !********!!!!!!!!!!!!************!!!!!!!!!!!!**************
  
  j = 0;  do bf3 = 1,nobasisfun; do bf4= bf3,nobasisfun;     j = j+1;
     if (bf3 == bf4) then;   flatdigmp(j,1:nLpure) = 0.5d0*flatdigmp(j,1:nLpure);     end if
  end do; end do

  do i=1,noggss;       if (int(ggss(1,i)) == int(ggss(2,i))) then;        ggss(9,i) = 0.5d0*ggss(9,i);         end if ;  end do
  do i=1,noggpp;       if (int(ggpp(1,i)) == int(ggpp(2,i))) then;        ggpp(9,i) = 0.5d0*ggpp(9,i);         end if ;  end do

! The scaling of newsortRg is done in the original routine!!
  call scf(shortol, shortcore,nobasisfun,noatoms, atoms,noelec, multiplicity, invintatomdist, coreHprint, consXprint,moleculename,&
          kmax, flatdigmp,Tintmat, totnomodel, nLpure, algchange, noggss, noggsp, noggpp, &
          ggss(1:20,1:noggss),ggsp(1:20,1:noggsp), ggpp(1:20,1:noggpp), basis,RRRRpreftable,L,&
          maxrampdegree,Tcutoff,quadpts, compareprint, newsortRg,lenmodelnewsort,maxrgpa,norampsperatom, totnoramps ,rampinfo)

  timecnt = timecnt+1;  call cpu_time(t(timecnt));print "(a120,f20.12,a5)", "Time for scf cal is", t(timecnt) -t(timecnt-1), "secs";
  deallocate(shortol, shortcore,basis,basright,atoms,invintatomdist,RRRRpreftable,ggss,ggsp,ggpp, Tintmat,flatdigmp);

  print *, "kmax = ", kmax; print *, "kmax is constant throughout the program"
  print *, "lenmdlconc = ", lenmdlconc;
  print *, "lenmdlconc is the length of the |Rg) --> Sum(R ) concentric that is used"
 
  print *, ""; print "(a44,a4,a20,a10)", "It is the end of the program calculating ", moleculename, " using basis set ", basisset
  print *, ""; print *, "****************** SUCCESSFUL FINISH************************************"; print *, ""

  print *, "totnomodel = ", totnomodel  

end program master
!************************************************************************************************
