      subroutine qswaveno(f,k)
      implicit none
c
      double precision f,k
c
      include 'qsglobal.h'
c
      integer n
      double complex ck,ck2
c
      ck=dcmplx(k,0.d0)
      ck2=dcmplx(k*k,0.d0)
c
      if(roair.gt.0.d0)then
        kpair=cdsqrt((ck+comega/cvpair)*(ck-comega/cvpair))
      else
        kpair=(0.d0,0.d0)
      endif
c
      do n=1,n0
        kp(n)=cdsqrt((ck+comega/cvp(n))*(ck-comega/cvp(n)))
        ks(n)=cdsqrt((ck+comega/cvs(n))*(ck-comega/cvs(n)))
        wb(n)=(2.d0,0.d0)*cmu(n)*ck
        wa(n)=cmu(n)*(ck2+ks(n)*ks(n))
      enddo
c
      return
      end
