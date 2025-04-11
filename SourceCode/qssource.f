      subroutine qssource(ros,vps,vss)
      implicit none
c
      double precision ros,vps,vss
c
      include 'qsglobal.h'
c
      integer i,istp
      double precision pi,pi2
c
      do istp=1,3
        do i=1,6
          sfct0(i,istp)=0.d0
          sfct1(i,istp)=0.d0
        enddo
      enddo
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
c     istp = 1
c     explosion source (m11=m22=m33=1)
c
      sfct0(1,1)=-1.d0/(pi2*ros*vps*vps)
      sfct1(4,1)=-(vss/vps)**2/pi
c
c     istp = 2
c     vertical-single-force (fz=1)
c
      sfct0(2,2)=1.d0/pi2
c
c     istp = 3
c     horizontal-single-force (fx=1)
c
      sfct0(4,3)=1.d0/pi2
      sfct0(6,3)=1.d0/pi2
c
      return
      end
