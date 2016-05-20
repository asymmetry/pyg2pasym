      program polin
      implicit none
      character(len=20) arg,addfn
      character(len=50) infile,outfile
      integer   pol,ios,i
      real*8    Ei,Ef,theta,xs_born,xs_rad
      
      addfn="000"
      pol=1
ccc    pol=0: parallel unpolarized
ccc        1: parallel polarized
ccc        2: perpendicular unpolarized
ccc        3: perpendicular polarized

      do i=1,iargc()
        call getarg(i,arg)
        if(i.eq.1) then
            read(arg,'(I2)') pol
        elseif(i.eq.2) then
            addfn=arg
        endif
      end do

      write(infile,"(A13A3A4)") "polelEiEflist",addfn,".dat"
      write(outfile,"(A5A3A4)") "polel",addfn,".dat"

      open(unit=1,file=infile)
      open(unit=2,file=outfile)

      i=0
      do 
        read (1,"(F7.2,X,F7.2,X,F7.5)",iostat=ios) 
     .      Ei,Ef,theta
        if (ios/=0) exit
        if(abs(theta).lt.1e-5.or.Ef.le.0) then
            write (2,*) 0,0
        else
            call polsig_el(Ei,Ef,theta,pol,xs_born,xs_rad)
            write (2,*) xs_born,xs_rad
            if(mod(i,1000).eq.0) then
c               print*,Ei,Ef,theta,pol
c               print*,xs_born,xs_rad
                print*,i
            endif
        endif
        i=i+1
      end do

      close(1)
      close(2)

      end

