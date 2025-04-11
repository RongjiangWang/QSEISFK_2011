      subroutine qsgetinp(unit)
      implicit none
      integer unit
c
      include 'qsglobal.h'
c
c     work space
c
      integer i,istp,j,l,l1,n,ierr
      double precision rr,z
      double precision depth,vsliquid
      double precision f1,v1
      double precision ros,vps,vss
      double precision resolut(3)
      character*180 comments
c
c     source parameters
c     =================
c
      call getdata(unit,comments)
      read(comments,*)zs
      zs=dmax1(0.d0,km2m*zs)
c
c     receiver parameters
c     ===================
c
      call getdata(unit,comments)
      read(comments,*)zr
      zr=km2m*zr
      call getdata(unit,comments)
      read(comments,*)f0,f1,df
      if(f0.le.0.d0)then
        stop ' Error: frequency <= 0!'
      else if(f0.ge.f1.or.df.le.0.d0)then
        nf=0
      else
        nf=idnint((f1-f0)/df)
        if(nf.gt.0)then
          df=(f1-f0)/dble(nf)
        endif
      endif
      call getdata(unit,comments)
      read(comments,*)v0,v1,dv
      if(v0.le.0.d0)then
        stop ' Error: slowness <= 0!'
      else if(v0.ge.v1.or.dv.le.0.d0)then
        nv=0
      else
        nv=idnint((v1-v0)/dv)
        if(nv.gt.0)then
          dv=(v1-v0)/dble(nv)
        endif
      endif
c
c     partial solution parameters
c     ===========================
c
      call getdata(unit,comments)
      read(comments,*)isurf
      if(isurf.lt.0.or.isurf.gt.2)then
        stop 'Error: wrong switch for filtering surface reflection!'
      else if(isurf.eq.2.and.zr.gt.0.d0)then
        stop 'Error: filtering surface multiples for zr > 0!'
      endif
      call getdata(unit,comments)
      read(comments,*)ipath,pathdepth
      pathdepth=pathdepth*km2m
      if(ipath.eq.1.and.(pathdepth.lt.zs.or.pathdepth.lt.zr))then
        print *,'warning: condition for path filter is not satisfied,'
        print *,'==> path filter will not be selected!'
        ipath=0
      endif
      call getdata(unit,comments)
      read(comments,*)npar
      if(npar.ge.1)then
        ipartial=1
        do i=1,npar
          call getdata(unit,comments)
          read(comments,*)zup(i),zlow(i),ipsv(i)
          if(ipsv(i).le.0.or.ipsv(i).ge.5)then
            stop ' Error in qsmain: wrong partial solution selection!'
          endif
          zup(i)=zup(i)*km2m
          zlow(i)=zlow(i)*km2m
        enddo
      endif
c
c     output file
c     ===========
c
      call getdata(unit,comments)
      read(comments,*)outputfile
      call getdata(unit,comments)
      read(comments,*)iout
c
c     global model parameters
c     =======================
c
      call getdata(unit,comments)
      read(comments,*)iflat
      call getdata(unit,comments)
      read(comments,*)(resolut(i),i=1,3)
      do i=1,3
        if(resolut(i).le.0.d0)resolut(i)=0.1d0
        resolut(i)=1.d-02*resolut(i)
      enddo
      call getdata(unit,comments)
      read(comments,*)l
      if(l.gt.lmax)then
        stop ' Error in input: to large number of layers!'
      endif
c
c     multilayered model parameters
c     =============================
c
      do i=1,l
        call getdata(unit,comments)
        read(comments,*)j,h(i),vp(i),vs(i),ro(i),qp(i),qs(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        h(i)=km2m*h(i)
        vp(i)=km2m*vp(i)
        vs(i)=km2m*vs(i)
        ro(i)=km2m*ro(i)
        if(vs(i).le.vspmin*vp(i))vs(i)=0.9d0*vspmin*vp(i)
      enddo
c
c     end of inputs
c     =============
c
      if(iflat.eq.1)then
c
c       flat earth transformation (Mueller, 1985)
c
        zs=rr0*dlog(rr0/(rr0-zs))
        zr=rr0*dlog(rr0/(rr0-zr))
        if(ipartial.eq.1)then
          do i=1,npar
            zup(i)=rr0*dlog(rr0/(rr0-zup(i)))
            zlow(i)=rr0*dlog(rr0/(rr0-zlow(i)))
          enddo
        endif
        if(ipath.eq.1)then
          pathdepth=rr0*dlog(rr0/(rr0-pathdepth))
        endif
c
        do i=1,l
          rr=rr0-h(i)
          h(i)=rr0*dlog(rr0/rr)
          vp(i)=vp(i)*rr0/rr
          vs(i)=vs(i)*rr0/rr
          ro(i)=ro(i)*(rr/rr0)**ndens
        enddo
      endif
c
c     determine upper und lower parameter values of each layer
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h(i).gt.h(i-1))then
          z1(l0)=h(i-1)
          vp1(l0)=vp(i-1)
          vs1(l0)=vs(i-1)
          ro1(l0)=ro(i-1)
          qp1(l0)=qp(i-1)
          qs1(l0)=qs(i-1)
c
          z2(l0)=h(i)
          vp2(l0)=vp(i)
          vs2(l0)=vs(i)
          ro2(l0)=ro(i)
          qp2(l0)=qp(i)
          qs2(l0)=qs(i)
          l0=l0+1
        else
          z1(l0)=h(i)
          vp1(l0)=vp(i)
          vs1(l0)=vs(i)
          ro1(l0)=ro(i)
          qp1(l0)=qp(i)
          qs1(l0)=qs(i)
        endif
      enddo
      z1(l0)=h(l)
      vp1(l0)=vp(l)
      vs1(l0)=vs(l)
      ro1(l0)=ro(l)
      qp1(l0)=qp(l)
      qs1(l0)=qs(l)
c
c     construction of sublayers at the cutoff frequency
c
      call qssublay(resolut,f1)
      write(*,*)' The layered model of source site:'
      write(*,'(7a)')'    no ','  z(km)  ',
     &               '  vp(km/s) ','  vs(km/s) ',' ro(g/cm^3)',
     &               '    qp   ','    qs'
      depth=0.d0
      do i=1,n0
        if(vs(i).le.vspmin*vp(i))then
          vsliquid=0.d0
        else
          vsliquid=vs(i)
        endif
        write(*,1000)i,depth/km2m,vp(i)/km2m,
     &               vsliquid/km2m,ro(i)/km2m,qp(i),qs(i)
        depth=depth+h(i)
        if(i.lt.n0)then
          write(*,1000)i,depth/km2m,vp(i)/km2m,
     &               vsliquid/km2m,ro(i)/km2m,qp(i),qs(i)
        endif
      enddo
c
      call qslayer(ierr)
      n=nno(ls)
      ros=ro(n)
      vps=vp(n)
      vss=vs(n)
      if(iflat.eq.1)then
        rr=rr0*dexp(-zs/rr0)
        ros=ros*(rr0/rr)**ndens
        vps=vps*rr/rr0
        vss=vss*rr/rr0
      endif
c
      call qssource(ros,vps,vss)
c
c     for partial solution only
c
      do i=1,lp
        n=nno(i)
        pup(i)=.true.
        pdw(i)=.true.
        if(vs(n).gt.vspmin*vp(n))then
          svup(i)=.true.
          svdw(i)=.true.
          sh(i)=.true.
        else
          svup(i)=.false.
          svdw(i)=.false.
          sh(i)=.false.
        endif
      enddo
      if(ipartial.eq.1)then
        z=zr
        do i=1,lp-1
          z=z+0.5d0*hp(i)
          do j=1,npar
            if(z.ge.zup(j).and.z.le.zlow(j))then
              if(ipsv(j).eq.1)pup(i)=.false.
              if(ipsv(j).eq.2)pdw(i)=.false.
              if(ipsv(j).eq.3)svup(i)=.false.
              if(ipsv(j).eq.4)svdw(i)=.false.
            endif
          enddo
          z=z+0.5d0*hp(i)
        enddo
      endif
      if(ipath.eq.1)then
        z=dmax1(zs,zr)
        lpath=max0(ls,lzr)
        do i=max0(ls,lzr)+1,lp
          z=z+hp(i-1)
          if(pathdepth.ge.z)lpath=i
        enddo
        if(lpath.eq.lp)then
          print *,'the depth limit for path filter is too large!'
          print *,'=> no signals in the given time window!'
          stop
        endif
      else
        lpath=0
      endif
1000  format(i5,f12.2,3f11.4,2f8.1)
c
      return
      end
