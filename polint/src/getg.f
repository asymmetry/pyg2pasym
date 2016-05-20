***********************get g1 g2 from sqlite********************
c for efficiency consideration, read txt file to memory first
c instead of reading the sqlite database
c please run outg12.py to generate data/wQg12.txt first before
c using it
c author: Pengjia Zhu

      subroutine GET_G(X,Q2,G1TMP,G2TMP)
      implicit none
      integer DB_N,DB_NW
      parameter(DB_N=126170)
      real*8 db_W(DB_N),db_Q2GeV(DB_N)
      real*8 db_g1(DB_N),db_g2(DB_N)
      real*8 X,Q2,W,nu,M,G1TMP,G2TMP
      real*8 Q2fc(2),Wfc(2),kv(4,4),kv2(2,4),kv3(4)
      real*8 x2x1,y2y1,x2x,xx1,y2y,yy1
      integer iQ2fc(2),iWfc(2),ifc(4)
      integer i,j
      common/DBVARS/ DB_NW,db_W,db_Q2GeV,db_g1,db_g2
      
      G1TMP=0
      G2TMP=0
      M=0.938272013
      nu=Q2/(2*M*X)
      W=sqrt(M*M+2*M*nu-Q2)*1000

      if(W.lt.1080.or.W.gt.2000.or.Q2.lt.0.or.Q2.gt.5) then
          return
      endif
     
      Wfc(1)=floor(W*2./10.)*5
      Wfc(2)=ceiling(W*2./10.)*5
      iWfc(1)=nint((Wfc(1)-1080)/5)
      iWfc(2)=nint((Wfc(2)-1080)/5)

      if(Q2.lt.0.0001) then
        Q2fc(1)=0
        Q2fc(2)=0.0001
        iQ2fc(1)=0
        iQ2fc(2)=1
      else if(Q2.lt.0.01) then
        Q2fc(1)=floor(Q2*1.e4)/1.e4
        Q2fc(2)=ceiling(Q2*1.e4)/1.e4
        iQ2fc(1)=nint((Q2fc(1)-0.0001)/0.0001+1)
        iQ2fc(2)=nint((Q2fc(2)-0.0001)/0.0001+1)
      else if(Q2.lt.0.1) then
        Q2fc(1)=floor(Q2*1.e3)/1.e3
        Q2fc(2)=ceiling(Q2*1.e3)/1.e3
        iQ2fc(1)=nint((Q2fc(1)-0.01)/0.001+100)
        iQ2fc(2)=nint((Q2fc(2)-0.01)/0.001+100)
      else
        Q2fc(1)=floor(Q2*1.e2)/1.e2
        Q2fc(2)=ceiling(Q2*1.e2)/1.e2
        iQ2fc(1)=nint((Q2fc(1)-0.1)/0.01+191)
        iQ2fc(2)=nint((Q2fc(2)-0.1)/0.01+191)
      endif
        
      ifc(1)=iQ2fc(1)*DB_NW+iWfc(1)+1 
      ifc(2)=iQ2fc(2)*DB_NW+iWfc(1)+1
      ifc(3)=iQ2fc(1)*DB_NW+iWfc(2)+1
      ifc(4)=iQ2fc(2)*DB_NW+iWfc(2)+1
   
      do i=1,4
c        kv(i,1)=db_W(ifc(i))
c        kv(i,2)=db_Q2GeV(ifc(i))
        kv(i,3)=db_g1(ifc(i))
        kv(i,4)=db_g2(ifc(i))
      enddo
    
      x2x1=Wfc(2)-Wfc(1)
      y2y1=Q2fc(2)-Q2fc(1)
      x2x=Wfc(2)-W
      xx1=W-Wfc(1)
      y2y=Q2fc(2)-Q2
      yy1=Q2-Q2fc(1)

      do i=1,2
        do j=3,4
            if(x2x1.eq.0) then
                kv2(i,j)=kv(i,j)
            else
                kv2(i,j)=x2x/x2x1*kv(i,j)+xx1/x2x1*kv(i+2,j)
            endif
        enddo
      enddo

      do j=3,4
        if(y2y1.eq.0)then
            kv3(j)=kv2(1,j)
        else
            kv3(j)=y2y/y2y1*kv2(1,j)+yy1/y2y1*kv2(2,j)
        endif
      enddo

      G1TMP=-kv3(3)
      G2TMP=-kv3(4)
      
      end

      subroutine OPENDB(Ebeam)
        implicit none
        integer DB_N,DB_NW,Ebeam
        parameter(DB_N=126170)
        real*8 db_W(DB_N),db_Q2GeV(DB_N)
        real*8 db_g1(DB_N),db_g2(DB_N)
        character line*100
        character g12file*100
        integer i,ios
        common/DBVARS/ DB_NW,db_W,db_Q2GeV,db_g1,db_g2

        DB_NW=185
        i=1
        write(g12file,"(A11I4A4)") "data/WQg12_",Ebeam,".txt"
        print*,"open g12 file: ",g12file
        open(unit=15,FILE=g12file)

        do
            read(15,'(A)',iostat=ios) line
            if (ios/=0) exit
            read(line,"(F7.2,X,F7.5,X,F12.9,X,F12.9)") 
     .          db_W(i),db_Q2GeV(i),db_g1(i),db_g2(i)
            i=i+1 
        enddo
        
      end
      
