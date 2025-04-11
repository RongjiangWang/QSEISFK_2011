      program qseisfk
      implicit none
c
      include 'qsglobal.h'
c
c     work space
c
      integer i,istp,l,n,lup,llw,lf,lv,runtime
      double precision pi2,rad2deg,f,k,v,am,ph,beta
      double complex ypsv(4,3),ysh(2)
      integer time
c
      double precision eps
      data eps/1.0d-50/
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#    QQQ    SSSS  EEEEE  III   SSSS   FFFFF  K  K    #'
      print *,'#   Q   Q  S      E       I   S       F      K K     #'
      print *,'#   Q Q Q   SSS   EEEE    I    SSS    FFFF   KK      #'
      print *,'#   Q  QQ      S  E       I       S   F      K K     #'
      print *,'#    QQQQ  SSSS   EEEEE  III  SSSS    F      K  K    #'
      print *,'#                                                    #'
      print *,'#                  (Version 2011)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           Deutsches GeoForschungsZentrum           #'
      print *,'#         Last modified: Potsdam, April, 2011        #'
      print *,'######################################################'
      print *,'                          '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      pi2=8.d0*datan(1.d0)
      rad2deg=360.d0/pi2
c
      open(10,file=inputfile,status='old')
      call qsgetinp(10)
      close(10)
c
      open(20,file=outputfile,status='unknown')
      if(iout.eq.1)then
        write(20,'(a)')'           f[Hz]         v[km/s]'
     &               //'         AmUz_Ex         PhUz_Ex'
     &               //'         AmUr_Ex         PhUr_Ex'
     &               //'         AmUz_Fz         PhUz_Fz'
     &               //'         AmUr_Fz         PhUr_Fz'
     &               //'         AmUz_Fh         PhUz_Fh'
     &               //'         AmUr_Fh         PhUr_Fh'
     &               //'         AmUt_Fh         PhUt_Fh'
      else
        write(20,'(a)')'           f[Hz]         v[km/s]'
     &               //'         ReUz_Ex         ImUz_Ex'
     &               //'         ReUr_Ex         ImUr_Ex'
     &               //'         ReUz_Fz         ImUz_Fz'
     &               //'         ReUr_Fz         ImUr_Fz'
     &               //'         ReUz_Fh         ImUz_Fh'
     &               //'         ReUr_Fh         ImUr_Fh'
     &               //'         ReUt_Fh         ImUt_Fh'
      endif
      do lf=0,nf
        f=f0+dble(lf)*df
        call qsqmodel(f)
        do lv=0,nv
          v=v0+dble(lv)*dv
          k=pi2*f/(v*km2m)
          call qswaveno(f,k)
c
          do istp=1,3
            do i=1,4
              ypsv(i,istp)=(0.d0,0.d0)
            enddo
          enddo
          do i=1,2
            ysh(i)=(0.d0,0.d0)
          enddo
c
          lup=1
          beta=0.d0
          do l=ls-1,lzr,-1
            n=nno(l)
            if(vs(n).gt.(1.d0+eps)*vspmin*vp(n))then
              beta=beta-dreal(ks(n))*hp(l)
            else
              beta=beta-dreal(kp(n))*hp(l)
            endif
            if(beta.le.dlog(eps))then
              goto 200
            endif
          enddo
c
          llw=lp
          beta=0.d0
          do l=ls+1,lp-1
            n=nno(l)
            if(vs(n).gt.(1.d0+eps)*vspmin*vp(n))then
              beta=beta-dreal(ks(n))*hp(l)
            else
              beta=beta-dreal(kp(n))*hp(l)
            endif
            if(beta.le.dlog(eps))then
              llw=l
              goto 100
            endif
          enddo
100       continue
c
          call qspsv(ypsv,k,lup,llw)
          call qssh(ysh,k,lup,llw)
c
          do istp=1,3
            do i=1,4
              ypsv(i,istp)=ypsv(i,istp)/dcmplx(k,0.d0)
            enddo
          enddo
          do i=1,2
            ysh(i)=ysh(1)/dcmplx(k,0.d0)
          enddo
c
200       continue
c
          if(iout.eq.1)then
            do istp=1,3
              do i=1,3,2
                am=cdabs(ypsv(i,istp))
                if(am.le.eps)then
                  ypsv(i,istp)=(0.d0,0.d0)
                else
                  ph=datan2(dimag(ypsv(i,istp)),dreal(ypsv(i,istp)))
     &              *rad2deg
                  ypsv(i,istp)=dcmplx(am,ph)
                endif
              enddo
            enddo
            am=cdabs(ysh(1))
            if(am.le.eps)then
              ysh(1)=(0.d0,0.d0)
            else
              ph=datan2(dimag(ysh(1)),dreal(ysh(1)))
     &          *rad2deg
              ysh(1)=dcmplx(am,ph)
            endif
            write(20,'(2E16.8,7(E16.8,f16.8))')f,v,
     &               ((ypsv(i,istp),i=1,3,2),istp=1,3),ysh(1)
          else
            do istp=1,3
              do i=1,3,2
                am=cdabs(ypsv(i,istp))
                if(am.le.eps)then
                  ypsv(i,istp)=(0.d0,0.d0)
                endif
              enddo
            enddo
            am=cdabs(ysh(1))
            if(am.le.eps)then
              ysh(1)=(0.d0,0.d0)
            endif
            write(20,'(16E16.8)')f,v,
     &               ((ypsv(i,istp),i=1,3,2),istp=1,3),ysh(1)
          endif
        enddo
      enddo
      close(20)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with qseisfk     #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
1001  format(2i7,E12.4,a)
1002  format(i4,a,E12.4,a,$)
1003  format(E12.5,$)
1004  format(2E12.4,$)
1005  format(2E12.4)
 500  stop
      end
