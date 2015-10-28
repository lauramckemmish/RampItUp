subroutine secondsortmodelRG(newsortRg, lenmodelnewsort, totnomodel,normmodelRGmat, lennormmodel, &
     nobasisfun, noatoms, maxmodellen, maxL,maxlenbf, maxrgpa)
  implicit none
  integer, intent(in)                    :: nobasisfun, noatoms,maxmodellen, maxL, lennormmodel(1:nobasisfun,1:nobasisfun,1:2)
  real*8, intent(in)                     :: normmodelRGmat(1:nobasisfun,1:nobasisfun,1:3,1: maxlenbf)
  integer, intent(in)                    :: maxlenbf
  real*8, intent(out)                    :: newsortRg(1:4,1:nobasisfun*noatoms*maxmodellen)
  integer, intent(out)                   :: totnomodel, lenmodelnewsort(1:4,1:10*nobasisfun,1:noatoms)
  integer                                :: i, j, k, Kmodel, st, fin, a, rampcenter, lensortedmodel(1:(1+maxL)**2,1:noatoms)
  integer   :: cnt(1:noatoms), stindex, finindex, bf1, bf2, maxrgpa, w, rc2
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This function sorts the model so thrampcenter it the list is in terms of the type of ramps in K-form pure ramps
!! then whrampcenter rampcenterom centre the ramp is on. 

!! The indexes of newsortRg are (coefficient, ramp degree, bf1, bf2, ramp centre)
!! The indexes of inclensortedmodel are (start index newsortRg, finish index newsortRg, number of this kind of ramp 
!!  (with same angular momentum and rampcenteromic center) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *, "in secondsortmodelRg"
  lenmodelnewsort = 0; totnomodel = 0;
  
  cnt = 0;
  do bf1=1,nobasisfun;  do bf2=bf1,nobasisfun
     if (lennormmodel(bf1,bf2,1) > 0) then;
         rampcenter = lennormmodel(bf1,bf2,2)
         if (rampcenter < 0) then
            rampcenter = -rampcenter;
               rc2 = lennormmodel(bf2,bf1,2)
               cnt(rc2) = cnt(rc2) + 1;
               lenmodelnewsort(1,cnt(rc2), rc2) = bf2
               lenmodelnewsort(2,cnt(rc2), rc2) = bf1
               print *, "IMPORTANT - Working in secondsortmodelRG", bf1,bf2,rampcenter,rc2
         end if        
         cnt(rampcenter) = cnt(rampcenter) + 1;
         lenmodelnewsort(1,cnt(rampcenter), rampcenter) = bf1
         lenmodelnewsort(2,cnt(rampcenter), rampcenter) = bf2
     end if
  end do; end do;

print *, "here 2 in secondsortRg"
  stindex = 0; 
  finindex = 0
  do rampcenter = 1,noatoms
    do j=1,cnt(rampcenter)
      bf1 = lenmodelnewsort(1,j,rampcenter)
      bf2 = lenmodelnewsort(2,j,rampcenter)
 if (bf1 == 0) then; call exit; end if
 if (bf2 == 0) then; call exit; end if
      stindex = finindex
      finindex = stindex + lennormmodel(bf1,bf2,1)
      lenmodelnewsort(3,j, rampcenter) = stindex
      lenmodelnewsort(4,j, rampcenter) = finindex

      do i=1,lennormmodel(bf1,bf2,1)
         newsortRg(1:2,stindex+i) = normmodelRGmat(bf1,bf2,1:2,i)
         newsortRg(3,stindex+i) = int(normmodelRGmat(bf1,bf2,3,i))    ! k angular momentum
         k = newsortRg(3,stindex+i)
         if (k==1) then;
             newsortRg(4,stindex+i) = 0            ! l angular momentum no
         else if (k > 1 .and. k <= 4) then;
             newsortRg(4,stindex+i) = 1    ! l angular momentum no
         else if (k > 4 .and. k <= 9) then;
             newsortRg(4,stindex+i) = 2    ! l angular momentum no
         else if (k > 9 .and. k <= 16) then;
             newsortRg(4,stindex+i) = 3    ! l angular momentum no
         else if (k > 16 .and. k <= 25) then;
             newsortRg(4,stindex+i) = 4    ! l angular momentum no
         else if (k > 25 .and. k <= 36) then;
             newsortRg(4,stindex+i) = 5    ! l angular momentum no
         else if (k > 36 .and. k <= 49) then;
             newsortRg(4,stindex+i) = 6    ! l angular momentum no
         else if (k > 49 .and. k <= 64) then;
             newsortRg(4,stindex+i) = 7    ! l angular momentum no
         end if
      end do  
      if (bf1==bf2) then;
         newsortRg(1,stindex+1:finindex) = 0.5*newsortRg(1,stindex+1:finindex)
      end if
    end do

   do j=cnt(rampcenter)+1,nobasisfun;
      lenmodelnewsort(1,j, rampcenter) = 0
      lenmodelnewsort(2,j, rampcenter) = 0
      lenmodelnewsort(3,j, rampcenter) = finindex
      lenmodelnewsort(4,j, rampcenter) = finindex
   end do
  end do
 
print *, "here near end of secondsortRg"
  totnomodel = finindex;
  maxrgpa = 0;
  do i=1,noatoms
     if (cnt(i) > maxrgpa) then; 
        maxrgpa = cnt(i)
     end if
  end do

print *, "end of secondsortRg"

  return
end subroutine secondsortmodelRG
