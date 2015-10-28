subroutine modelRGdigestion(normmodelRGmat, lennormmodel, basisfunmodelRGmat, lenbfmodel, &
     bfmodelids, maxbflen,nobasisfun,no) 
  implicit none
  real*8, intent(out)                     :: normmodelRGmat(1:nobasisfun,1:nobasisfun,1:3,1:maxbflen)
  integer, intent(out)                    :: lennormmodel(1:nobasisfun,1:nobasisfun,1:2)
  real*8, intent(in)                      :: basisfunmodelRGmat(1:no,1:3,1:maxbflen), bfmodelids(1:no,1:4)
  integer, intent(in)                     :: lenbfmodel(1:no), maxbflen, nobasisfun, no
  real*8                                  :: nc, StoIconversion, temparray(1:3,1:maxbflen)
  integer                                 :: nomodel, extracnt, i, j, k, m, bf1, bf2, swit, minbf, maxbf
  integer :: actk, actn, cnti, cntx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine digests all the various components of a basis function pair model and 
!!! makes one set of model components for each basis function pair (that requires them)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 lennormmodel = 0;
  StoIconversion = 1d0/(2d0*sqrt(dacos(-1d0)));
  
  do i=1,no
     bf1 = bfmodelids(i,1);     bf2 = bfmodelids(i,2)
     minbf = min(bf1, bf2);     maxbf = max(bf1, bf2)
     nc  = bfmodelids(i,4)*StoIconversion
     if (bf1 > 0 .and. bf2 > 0) then   
        ! this is the atomic centre
        if (bfmodelids(i,3) .ne. 0 .and. lennormmodel(minbf,maxbf,2) .gt. 0 .and.&
                           bfmodelids(i,3) .ne. lennormmodel(minbf,maxbf,2)) then;
          lennormmodel(minbf,maxbf,2) = - lennormmodel(minbf,maxbf,2)
          minbf = max(bf1,bf2); maxbf = min(bf1,bf2);
        end if
        lennormmodel(minbf,maxbf,2) = bfmodelids(i,3)
        nomodel = lennormmodel(minbf,maxbf,1);
        
        if (nomodel ==0) then
           lennormmodel(minbf,maxbf,1) = lenbfmodel(i)
           normmodelRGmat(minbf,maxbf,1  ,1:lenbfmodel(i)) = basisfunmodelRGmat(i,1  ,1:lenbfmodel(i))*nc
           normmodelRGmat(minbf,maxbf,2:3,1:lenbfmodel(i)) = basisfunmodelRGmat(i,2:3,1:lenbfmodel(i))
           do j=1,lenbfmodel(i)
              if (int(basisfunmodelRGmat(i,3,j)) == 0) then; print *, "K = 0 in modelRGdigestion";              end if
              if (int(normmodelRGmat(minbf,maxbf,3,j)) == 0) then; print *, "K = 0 in modelRGdigestion";               end if; 
           end do
        else
           extracnt = 0;
           do k=1,lenbfmodel(i)
              swit = 0;
              do m = 1,nomodel 
                 if (normmodelRGmat(minbf,maxbf,2,m) == basisfunmodelRGmat(i,2,k) .and. &
                      normmodelRGmat(minbf,maxbf,3,m) == basisfunmodelRGmat(i,3,k)) then
                    if (swit == 0) then
                       swit = 1; normmodelRGmat(minbf,maxbf,1,m) = normmodelRGmat(minbf,maxbf,1,m) + basisfunmodelRGmat(i,1,k)*nc
                    else;
                       print *, "there is a problem in modelRGdigestion";     call exit;
                    end if;
                 end if;
              end do;
              if (swit == 0) then;
                 extracnt = extracnt + 1;
                 normmodelRGmat(minbf,maxbf,1,nomodel+extracnt) = basisfunmodelRGmat(i,1,k)*nc
                 normmodelRGmat(minbf,maxbf,2:3,nomodel+extracnt) = basisfunmodelRGmat(i,2:3,k)
              end if;

        end do;
! This part of the code ensures that all models are ordered correctly
            temparray(1:3,1:nomodel+extracnt) = normmodelRgmat(minbf,maxbf,1:3,1:nomodel+extracnt);
            cntx  = 0;
            do actk=1,25
             do actn=0,100
              do cnti=1,nomodel+extracnt
                if (temparray(2,cnti) == actn .and. temparray(3,cnti) == actk) then;
                   cntx = cntx + 1;
                   normmodelRgmat(minbf,maxbf,1:3,cntx) = temparray(1:3,cnti)
                end if
              end do
             end do
           end do;
          if (cntx .ne. lennormmodel(minbf,maxbf,1) + extracnt) then;
              print *, "NOT EQUAL In MODELRGDIGESTION"
          end if;

           lennormmodel(minbf,maxbf,1) = lennormmodel(minbf,maxbf,1) + extracnt

           if (lennormmodel(minbf,maxbf,1) > maxbflen) then; print *, "Model length out of control";
              print *, "Calling exit"; call exit; 
           end if
        end if
     end if
  end do
  
  do i=1,nobasisfun; do j=i+1,nobasisfun;
     if (lennormmodel(j,i,1) > 0) then
        print *, i, j; print *, "modelling the wrong shell pair"
        print *, "IMPORTANT - working for ", i, j
     end if;
  end do; end do;
  
  return
  
end subroutine modelRGdigestion
