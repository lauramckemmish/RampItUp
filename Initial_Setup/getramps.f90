subroutine getramps(norampsperatom, totnoramps ,temprampinfo,basright,nobasisfun,noatoms)
 implicit none
 integer, intent(in)   :: noatoms, nobasisfun
 real*8, intent(in)    :: basright(1:30,1:nobasisfun)
 integer, intent(out)   :: temprampinfo(1:2,1:4*noatoms), norampsperatom(1:noatoms), totnoramps
 integer               :: bfid, typebf, atid

 totnoramps = 0; 
 norampsperatom = 0;
 do bfid=1,nobasisfun
     typebf = int(basright(3,bfid))
     atid   = int(basright(1,bfid))
     if (typebf == 0 .or. (typebf .ge. 5 .and. typebf .le. 11)) then
         totnoramps = totnoramps + 1;
         norampsperatom(atid) = norampsperatom(atid) + 1
         temprampinfo(1:2, totnoramps) = (/bfid,atid/)
     end if
 end do

end subroutine getramps