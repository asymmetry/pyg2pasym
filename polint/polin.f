      program polin
      implicit none
      character(len=20) arg,addfn
      character(len=50) infile,outfile
      integer   pol,ios,i,Ebeam
      real*8    Ei,Ef,theta,xs_born,xs_rad,W2_START
      COMMON/XGOOF_2/W2_START 
      
      addfn="000"
      pol=1
ccc    pol=0: parallel unpolarized
ccc        1: parallel polarized
ccc        2: perpendicular unpolarized
ccc        3: perpendicular polarized

      do i=1,iargc()
        call getarg(i,arg)
        if(i.eq.1) then
            read(arg,'(I4)') Ebeam
        elseif(i.eq.2) then
            read(arg,'(I2)') pol
        elseif(i.eq.3) then
            addfn=arg
        endif
      end do

      write(infile,"(A13A3A4)") "polinEiEflist",addfn,".dat"
      write(outfile,"(A5A3A4)") "polin",addfn,".dat"

      open(unit=1,file=infile)
      open(unit=2,file=outfile)

      call OPENDB(Ebeam)
      i=0
      do 
        read(1,"(F7.2,X,F7.2,X,F7.2,X,F7.5,X,F8.5)",iostat=ios)
     .      Ei,Ef,theta,xs_born
        if (ios/=0) exit
c        W2_START=1.2483011590040012
        W2_START=1.1664
        if(abs(theta).lt.1e-5.or.Ef.le.0) then
            write (2,*) 0
        else
            call polsig_2(Ei,Ef,theta,0,pol,xs_born,xs_rad)
            write (2, *) xs_rad
            if(mod(i,1000).eq.0)then
c               print*,Ei,Ef,theta,pol,xs_born
c               print*,xs_rad
                print*,i
            endif
        endif
        i=i+1
      end do

      close(1)
      close(2)


      end

