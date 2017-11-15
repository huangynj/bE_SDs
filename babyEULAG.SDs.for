      program bubble_sd
cc
      include 'param.grid'
      include 'param.sds'

cc  nember of time steps:
      parameter(ntime0=10*60)
cc  major dimensions
cc HAVE TO BE CHANGED IN ALL ROUTINES:
cc  grid (for regular dzeta isosurfaces)
cc small time step required for SDs:
      data dx,dz,dt /20.,20.,1.0/
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      dimension xx(nx),zz(nz)

      dimension xcntrm(ntime0),zcntrm(ntime0),cntti(ntime0)
      dimension qcav(ntime0),qcst(ntime0)
      dimension concm(ntime0),radcm(ntime0),stdcm(ntime0)

c SD atributes: two-time level
c SD atributes: position
      dimension xp(np),zp(np)
      dimension xpn(np),zpn(np)
c               number of SDs, cloud water, multiplicity, radius
      dimension qcp(np),plic(np),radp(np)
      dimension qcpn(np),plicn(np),radpn(np)

CC  MODEL VARIABLES
cc thermodynamics
      dimension theta(nx,nz),qv(nx,nz),qc(nx,nz),nc(nx,nz)
      dimension t(nx,nz)
cc limiting through S; routine adjust
      dimension sup(nx,nz),fsup(nx,nz),dpdz(nz)
      common /p_deriv/ dpdz
cc dynamics
      dimension ux(nx,nz),uz(nx,nz)          ! current time level
      dimension uxp(nx,nz),uzp(nx,nz)       ! previous time level
cc  forces for model variables
      dimension ft(nx,nz),fx(nx,nz),fy(nx,nz),fz(nx,nz),fqv(nx,nz)
cc  advective velocities:
      dimension uxa(nx+1,nz),uza(nx,nz+1)

cc profiles  
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      common /prof_a/ tau(nz)

cc arrays for w diagnosis; adjust scheme...
      real wdum(nx,nz)
      real thwd(nx,nz),twd(nx,nz)  ! for w doagnosis...
      real thwd0(nx,nz),twd0(nx,nz)  ! for w doagnosis...

cc required variables:
      dimension den(nx,nz),p(nx,nz),scr1(nx,nz),scr2(nx,nz)

cc constants
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

c grid in supers to get the same mulitplicity
      dimension sup_sdb(nppg),rad_sdb(nppg)
c parameters for activation:
      common /activ/ sup_sdb,rad_sdb,dcon

      call opngks
      call gsclip(0)

cc adjustment as in Grabowski and Morrison (MWR 2008):
cc  1 - use adjust
cc  0 - no adjust
       iadjust=1
c       iadjust=0
cc w for adjustment: uz (0) or diagnosed through T and theta (1):
c       iwdiag=1
       iwdiag=0

cc grid:
      time=0.
      dxi=1./dx
      dzi=1./dz
      dti=1./dt
      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
cc zz is regular (ustretched) dzeta grid
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo

      call zstrtch(zz,nz,dz)

cc initialize moisture parameters:
      call moist_init

cc initialize model profiles:
      call prof_init

           if(iadjust.eq.1) then
       e=-cp/rg
       do k=2,nz-1
       thetmep=th_e(k+1)/tm_e(k+1)
       prekp=1.e5*thetmep**e
       thetmem=th_e(k-1)/tm_e(k-1)
       prekm=1.e5*thetmem**e
       dpdz(k)=(prekp-prekm)/(height(k+1)-height(k-1))
       enddo
       dpdz(1)=dpdz(2)
       dpdz(nz)=dpdz(nz-1)
           endif
cc read in and initialize CCN info: 
      call activat(0.,0.,0)
      sup00=ssmax*1.e-2
      call activat(sup00,conmax,1)
      print*,'CCN: max conc (kg**-1): ',conmax
      dcon=conmax/nppg
cc grid in S for SDs:
       do ip=1,nppg
       con=float(ip)*dcon
       call activat(supsat,con,2)
       sup_sdb(ip)=supsat
       rad_sdb(ip)=8.e-10/supsat
       print*,'--ip,con,sup(%),r(mic): ',
     1   ip,dcon,con,sup_sdb(ip)*1.e2,rad_sdb(ip)*1.e6
       enddo

       do k=1,nz
       do i=1,nx
       den(i,k)=rho0(k)*gac(k)
       enddo
       enddo

cc initial fields
      a=rg/rv
      c=hlatv/cp
      b=hlatv/(rv*tt0)
      d=hlatv/rv
      e=-cp/rg

       zcen=800.
       xcen=(nx-1)*dx/2.
c       rad1=200.
c       rad2=300.
       rad1=250.
       rad2=350.
       pi=4.*atan(1.)

       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       do i=1,nx

        radd=(xx(i)-xcen)**2+(zz(k)-zcen)**2
        radd=sqrt(radd)
         if(radd.gt.rad2) rh=0.2
         if(radd.le.rad1) rh=1.
         if(radd.gt.rad1 .and. radd.le.rad2)
     1   rh=.2 + .8*(cos(pi/2. * (radd-rad1)/(rad2-rad1)))**2.

        theta(i,k)=th_e(k)

         thi=1./theta(i,k)
         y=b*thetme*tt0*thi
         ees=ee0*exp(b-y)
         qvs=a*ees/(prek-ees)
         del=rh*qvs - qv_e(k)

        qv(i,k)=qv_e(k)    + del

        sup(i,k)=qv(i,k)/qvs - 1.

        qc(i,k)=0.

        ux(i,k)=ux_e(k)
        uz(i,k)=0.
        uxp(i,k)=ux_e(k)
        uzp(i,k)=0.

        ft(i,k)=0.
        fx(i,k)=0.
        fy(i,k)=0.
        fz(i,k)=0.

        fqv(i,k)=0.
        fsup(i,k)=0.

        p(i,k)=0.

       enddo
       enddo

         npx=0
cc initial SDs:
         do ip=1,np
         xp(ip)=0.
         zp(ip)=0.
         qcp(ip)=0.
         plic(ip)=0.
         radp(ip)=0.
         qcpn(ip)=0.
         plicn(ip)=0.
         radpn(ip)=0.
         enddo

cc total water and mse:
       sum=0.
       sum1=0.
       sum2=0.
       do i=1,nx-1
       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       coe=1.
       if(k.eq.1 .or. k.eq.nz) coe=0.5
       sum=sum+coe*den(i,k)*(qv(i,k)+qc(i,k))
       ttt=theta(i,k)/thetme
       amse=cp*ttt + g*height(k) + hlatv*qv(i,k)
       thee=theta(i,k) + hlatv/cp*thetme*qv(i,k)
       sum1=sum1+amse
       sum2=sum2+thee
       enddo
       enddo
       print*,'------ total water,mse,the:',sum,sum1,sum2

cc plot initial fields:
       call diagno_1(ux,uz,theta,nx,nz,scr1,scr2,den)
       call diagno_2(qv,qc,nx,nz,scr1,scr2,den)
       call diagno_3(theta,qv,scr1)
       call plot_1(ux,uz,theta)
       call plot_2(theta,qv,qc)

ccc save initial data:
       write(17) time,nx,nz,nppg,np
       write(17) ux,uz,theta,qv,qc
       
CCCC MARCH FORWARD IN TIME:
ccc different random number realization:
ccc commented: run1
c            do ir=1,1171    ! run2
cc            do ir=1,133171    ! run3
c            r1=rand()
c            enddo

              ntime=ntime0
              do itime=1,ntime   ! TIME LOOP
               print*,'*** itime, time: ',itime,time

cc extrapolate in time to get advective momentums:
       call velprd_1(ux,uxp,uxa,uz,uzp,uza,nx,nz,den)

cc transport SDs:
      call sd_adv(ux,uxp,uz,uzp,xp,zp,npx)
          
cc save previous velocities:
       do i=1,nxz
        uxp(i,1)=ux(i,1)
        uzp(i,1)=uz(i,1)
       enddo

           if(iadjust.eq.1) then
cc derive new absolute S
       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       do i=1,nx
         thi=1./theta(i,k)
         y=b*thetme*tt0*thi
         ees=ee0*exp(b-y)
         qvs=a*ees/(prek-ees)
         sup(i,k)=qv(i,k)-qvs  ! absolute S
c
           if(sup(i,k)/qvs.gt.ssmax*1.e-2) then
           print*,'sup(i,k).gt.ssmax'
           print*,'i,k,s,smax: ',i,k,sup(i,k),ssmax
           call clsgks
           stop 'sup(i,k).gt.ssmax'
           endif
c
       sup(i,k)=sup(i,k)+dt*fsup(i,k)
       enddo
       enddo
             endif

cc add half of the force for velocities:
       do i=1,nxz
        ux(i,1)   =   ux(i,1)+.5*dt*fx(i,1)
        uz(i,1)   =   uz(i,1)+.5*dt*fz(i,1)
cc theta and qv: 1st order
        theta(i,1)=theta(i,1)+dt*ft(i,1)
        qv(i,1)   =   qv(i,1)+dt*fqv(i,1)
       enddo

CC ADVECTION:
c liner: 1-iord=1, 0-iord prescribed inside mpdata
        liner=0
c        if(itime/10*10.eq.itime) liner=1

cc advect velocities:
       call mpdat_2d(uxa,uza,   ux,den,1,liner)
       call mpdat_2d(uxa,uza,   uz,den,2,liner)
cc get new thermodynamic variables:
       call mpdat_2d(uxa,uza,theta,den,3,liner)
       call mpdat_2d(uxa,uza,   qv,den,4,liner)

       if(iadjust.eq.1) then
       call mpdat_2d(uxa,uza,  sup,den,4,liner)

       if(iwdiag.eq.1) then
cc w diagnosis for adjustment...
      do k=1,nz
       thetme=th_e(k)/tm_e(k)
      do i=1,nx
      thwd(i,k)=theta(i,k)
      twd(i,k)=theta(i,k)/thetme
      thwd0(i,k)=theta(i,k)
      twd0(i,k)=twd(i,k)
      enddo
      enddo

      call mpdat_2d(uxa,uza,thwd,den,3,liner)
      call mpdat_2d(uxa,uza,twd,den,3,liner)

      do i=1,nxz
      thwd(i,1)=(thwd(i,1)-thwd0(i,1))/dt
      twd(i,1)=(twd(i,1)-twd0(i,1))/dt
      enddo

c wdum = grid-scale vertical velocity
      do k=1,nz
      thetme=th_e(k)/tm_e(k)
      do i=1,nx
      wdum(i,k)=cp/gg * (-thwd(i,k)/thetme + twd(i,k))
      enddo
      enddo
       else    ! if(iwdiag.eq.1) then
      do k=1,nz
      do i=1,nx
      wdum(i,k)=uz(i,k)
      enddo
      enddo
       endif    ! if(iwdiag.eq.1) then

       endif ! iadjust

ccc 1st, finish dynamics...
cc save velocities after advection into advective velocities:
       do i=1,nx
       do k=1,nz
       uxa(i,k)=ux(i,k)
       uza(i,k)=uz(i,k)
       enddo
       enddo
cc derive qc for buoyancy:
       do i=1,nxz
       qc(i,1)=0.
       enddo
           do ip=1,npx
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
          qc(i,k)=qc(i,k)+qcp(ip)
           enddo

cc add buoyancy
       epsb=rv/rg-1.
       do k=1,nz
       do i=1,nx
       scr1(i,k)=gg*( (theta(i,k)-th_e(k))/th0(k)
     *   + epsb*(qv(i,k)-qv_e(k))-qc(i,k) )
       enddo
       enddo

cc filter in vertical
cc       call integz(scr1,scr2,nx,nz)
       call integxz(scr1,scr2,nx,nz)

cc apply
       do i=1,nxz
       uz(i,1) = uz(i,1)+.5*dt*scr1(i,1)
       enddo

cc calculate pressure gradient force:
c      epp=1.e-6
      epp=1.e-7
      itp=100
      call gcrk_1(p,scr1,scr2,ux,uz,nx,nz,itp,epp)
      call prforc_1(p,scr1,scr2,ux,uz,nx,nz)
      do i=1,nxz
      ux(i,1)=scr1(i,1)
      uz(i,1)=scr2(i,1)
      enddo

cc calculate velocity forces (using saved velocities after advection):
       do k=1,nz
       do i=1,nx
       fx(i,k)=(ux(i,k)-uxa(i,k))  *2./dt
       fz(i,k)=(uz(i,k)-uza(i,k))  *2./dt
       enddo
       enddo

ccc 2nd, derive new forces for thermodynamics

cc adjust temp and qv to match sss:
       if(iadjust.eq.1)
     1 call adjust(theta,qv,qc,sup,npx,xp,zp,qcp,plic,radp)

cc derive new S:
        if(iadjust.eq.0) then
       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       do i=1,nx
         thi=1./theta(i,k)
         y=b*thetme*tt0*thi
         ees=ee0*exp(b-y)
         qvs=a*ees/(prek-ees)
         sup(i,k)=qv(i,k)/qvs-1.
           if(sup(i,k).gt.ssmax*1.e-2) then
           print*,'sup(i,k).gt.ssmax'
           print*,'i,k,s,smax: ',i,k,sup(i,k),ssmax
           call clsgks
           stop 'sup(i,k).gt.ssmax'
           endif
       enddo 
       enddo 
         endif   ! iadjust

        if(iadjust.eq.1) then
       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       do i=1,nx
         thi=1./theta(i,k)
         y=b*thetme*tt0*thi
         ees=ee0*exp(b-y)
         qvs=a*ees/(prek-ees)
         sup(i,k)=qv(i,k)-qvs
           if(sup(i,k)/qvs.gt.ssmax*1.e-2) then
           print*,'sup(i,k).gt.ssmax'
           print*,'i,k,s,smax: ',i,k,sup(i,k),ssmax
           call clsgks
           stop 'sup(i,k).gt.ssmax'
           endif
       enddo 
       enddo 
         endif   ! iadjust

cc set forces to zero:
       do i=1,nxz
       ft(i,1)=0.
       fqv(i,1)=0.
       fsup(i,1)=0.
       enddo
cc SDs thermodynamics; NOTE: theta and qv tendencies are brought back,
cc     SD properties are imediately updated...
       call sd_phy(theta,qv,sup,fsup,wdum,ft,fqv,
     1             npx,xp,zp,qcp,plic,radp,iadjust)

cc randomly re-position droplets
cc or every drep sec
c       drep=dt
c       drep=60.
       drep=600000000. ! never...
       irep=nint(drep/dt)
       if(itime/irep*irep.eq.itime) call sd_rep(xp,zp,npx)

cc update clock (in minutes...)
       time = float(itime)*dt/60. 

ccc diagno every time step
c       call diagno_1(ux,uz,theta,nx,nz,scr1,scr2,den)
c       call diagno_2(qv,qc,nx,nz,scr1,scr2,den)
c       call diagno_3(theta,qv,scr1)
c       call diagno_4(npx,xp,zp,qcp,plic,radp)

cc output and plot:
       dtout=1. ! in min
       dtape=1.  ! in min
       if(amod(time+.1*dt/60.,dtout).lt.0.5*dt/60.) then
ccc plot selected fields:
       call plot_1(ux,uz,theta)
       call plot_2(theta,qv,qc,npx)
       call plot_3(qc,xp,zp,qcp,plic,radp,npx)
cc analysis of output:
       print*,'   '
       call diagno_1(ux,uz,theta,nx,nz,scr1,scr2,den)
       call diagno_2(qv,qc,nx,nz,scr1,scr2,den)
       call diagno_3(theta,qv,scr1)
       call diagno_4(npx,xp,zp,qcp,plic,radp)
cc 
cc total water and mse:
       sum=0.
       sum1=0.
       sum2=0.
       do i=1,nx-1
       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       coe=1.
       if(k.eq.1 .or. k.eq.nz) coe=0.5
       sum=sum+coe*den(i,k)*(qv(i,k)+qc(i,k))
       ttt=theta(i,k)/thetme
       amse=cp*ttt + g*height(k) + hlatv*qv(i,k)
       thee=theta(i,k) + hlatv/cp*thetme*qv(i,k)
       sum1=sum1+amse
       sum2=sum2+thee
       enddo
       enddo
       print*,'------ total water,mse,the:',sum,sum1,sum2

cccc mean qc and S in cloudy volumes:
       sum1=0.
       sum2=0.
       sum3=0.
       sum4=0.
       do i=1,nx-1
       do k=1,nz
       if(qc(i,k).gt.1.e-5) then
       sum1=sum1+qc(i,k)
       sum2=sum2+1.
       endif
       if(qc(i,k).gt.1.e-4) then
       sum3=sum3+sup(i,k)*100.
       sum4=sum4+1.
       endif
       enddo
       enddo
       if(sum2.gt.0.) then
       qcmean=sum1/sum2*1.e3
       else
       qcmean=0.
       endif
       if(sum4.gt.0.) then
       smean=sum3/sum4
       else
       smean=0.
       endif
       print*,'------ mean qc in cloudy(10**-5) volumes:',qcmean
       print*,'------ mean S(%) in cloudy(10**-4) volumes:',smean

       endif
       if(amod(time+.1*dt/60.,dtape).lt.0.5*dt/60.) then
       write(17) time,nx,nz,nppg,np
       write(17) ux,uz,theta,qv,qc
       write(17) npx
        do ip=1,npx
        write(17) xp(ip),zp(ip)
        write(17) qcp(ip),plic(ip),radp(ip)
        enddo
       print*,'wrote tape for t = ',time
       endif

      cntti(itime)=time
cc center of mass of qc:
      sum1=0.
      sum2=0.
      sum3=0.
      do i=1,nx
      xxx=float(i-1)*dx
      do k=1,nz
      sum1=sum1 + qc(i,k)
      sum2=sum2 + height(k)*qc(i,k)
      sum3=sum3 +       xxx*qc(i,k)
      enddo
      enddo
      if(sum1.gt.1.e-3) then
      zcntrm(itime)=sum2/sum1  *1.e-3
      xcntrm(itime)=sum3/sum1  *1.e-3
      else
      zcntrm(itime)=0.
      xcntrm(itime)=xx(nx)/2. * 1.e-3
      endif

cc conditionally-averaged qc:
      sum1=0.
      sum2=0.
      do i=1,nx
      do k=1,nz
      if(qc(i,k).ge.1.e-5) then
      sum1=sum1 + 1.
      sum2=sum2 + qc(i,k)
      endif
      enddo
      enddo
      if(sum1.gt.1.) then
      qcav(itime)=sum2/sum1  *1.e3
      else
      qcav(itime)=0.
      endif
cc std dev:
      sum1=0.
      sum2=0.
      do i=1,nx
      do k=1,nz
      if(qc(i,k).ge.1.e-5) then
      sum1=sum1 + 1.
      sum2=sum2 + (qc(i,k)*1.e3-qcav(itime))**2
      endif
      enddo
      enddo
      if(sum1.gt.1.) then
      qcst(itime)=sqrt(sum2/sum1) /qcav(itime) 
      else
      qcst(itime)=0.
      endif

       sum1=0.
       sum2=0.
          ic=(xcntrm(itime)*1.e3 + dx/2.)/dz + 1
          kc=(zcntrm(itime)*1.e3 + dz/2.)/dz + 1
       do ip=1,npx
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
        if(i.eq.ic .and. k.eq.kc) then
          con=plic(ip)
          sum1=sum1 + con
          sum2=sum2 + radp(ip)*con
        endif
       enddo
      if(sum1.gt.1.) then
      concm(itime)=sum1
      radcm(itime)=sum2/sum1
      else
      concm(itime)=0.
      radcm(itime)=0.
      endif

       sum1=0.
       sum2=0.
          ic=(xcntrm(itime)*1.e3 + dx/2.)/dz + 1
          kc=(zcntrm(itime)*1.e3 + dz/2.)/dz + 1
       do ip=1,npx
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
        if(i.eq.ic .and. k.eq.kc) then
          sum1=sum1 + 1.
          sum2=sum2 + (radp(ip)-radcm(itime))**2
        endif
       enddo
      if(sum1.gt.1.) then
      stdcm(itime)=sqrt(sum2/sum1)
      else
      stdcm(itime)=0.
      endif

      concm(itime)=concm(itime)*1.e-6 
      radcm(itime)=radcm(itime)*1.e6
      stdcm(itime)=stdcm(itime)*1.e6
      print*,'conc,rad,std: ',concm(itime),radcm(itime),stdcm(itime)

           enddo      ! TIME LOOP

cc    finished...

cc plot cntrm:
      call setusv('LW',2000)
      call set(.12,.52,.1,.5,0.,10.,0.,2.,1)           
      call labmod('(f3.0)','(f3.0)',4,4,2,2,20,20,0)     
      call periml(2,5,2,5)                           
      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)      
      call curved(cntti,zcntrm,ntime)                  
      call set(.58,.98,.1,.5,0.,10.,1.,2.,1)           
      call labmod('(f3.0)','(f3.1)',4,4,2,2,20,20,0)     
      call periml(2,5,2,5)                           
      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)      
      call curved(cntti,xcntrm,ntime)                  
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.32,0.03, 'time (min)', 0.014,0.,0)
      CALL plchhq(.78,0.03, 'time (min)', 0.014,0.,0)
      CALL plchhq(.03,0.30, 'x,z qc center of mass (km)', 0.014,90.,0)

cc plot qcav:
      call setusv('LW',2000)
      call set(.12,.52,.55,.95,0.,10.,0.,1.,1)           
      call labmod('(f3.0)','(f3.0)',4,4,2,2,20,20,0)     
      call periml(2,5,2,5)                           
      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)      
      call curved(cntti,qcav,ntime)                  
      call set(.58,.98,.55,.95,0.,10.,0.,1.,1)           
      call labmod('(f3.0)','(f3.1)',4,4,2,2,20,20,0)     
      call periml(2,5,2,5)                           
      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)      
      call curved(cntti,qcst,ntime)                  
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.03,0.75, 'av and st dev qc (g/kg)', 0.014,90.,0)
      call frame                                   

        call clsgks
        stop
        end

      subroutine noise(ff,nx,nz,nzn)
      dimension ff(nx,nz)
c      double precision rand

      do i=1,nx
      do k=1,nz
      ff(i,k)=0.
      enddo
      enddo

      do i=1,nx-1
      do k=2,nzn
      ff(i,k)=2.*(rand()-.5)
      enddo
      enddo
      do k=1,nzn
      ff(nx,k)=ff(1,k)
      enddo

      return
      end

      subroutine velprd_1(ux,uxp,uxa,uz,uzp,uza,nx,nz,rho)
      dimension ux(nx,nz),uz(nx,nz)     
      dimension uxp(nx,nz),uzp(nx,nz)   
      dimension uxa(nx+1,nz),uza(nx,nz+1),rho(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

       do k=1,nz
       do i=2,nx
       uxa(i,k) =(0.75*(ux(i-1,k)*rho(i-1,k)+ux(i,k)*rho(i,k)) 
     .          - 0.25*(uxp(i-1,k)*rho(i-1,k)+uxp(i,k)*rho(i,k)) )*dt/dx 
       enddo
cc cyclic in horizontal
       uxa(1,k) = uxa(nx,k)
       uxa(nx+1,k) = uxa(2,k)
       enddo
          
       do i=1,nx
       do k=2,nz
       uza(i,k) =(0.75*(uz(i,k-1)*rho(i,k-1)+uz(i,k)*rho(i,k)) 
     .          - 0.25*(uzp(i,k-1)*rho(i,k-1)+uzp(i,k)*rho(i,k)) )*dt/dz 
       enddo
cc zero flux in vertical
       uza(i,1) = - uza(i,2)
       uza(i,nz+1) = - uza(i,nz)
       enddo

       return
       end

      subroutine mpdat_2d(u1,u2,x,h,iflg,liner)
      include 'param.grid'
      parameter(n1=nx+1,n2=nz+1)
      parameter(n1m=n1-1,n2m=n2-1)
      dimension u1(n1,n2m),u2(n1m,n2),x(n1m,n2m),h(n1m,n2m)
      common// v1(n1,n2m),v2(n1m,n2),f1(n1,n2m),f2(n1m,n2),
     *         cp(n1m,n2m),cn(n1m,n2m),
     *         mx(n1m,n2m),mn(n1m,n2m)
      real mx,mn
      parameter(iord0=2,isor=1,nonos=1,idiv=0)
      data ep/1.e-12/
c
      donor(y1,y2,a)=cvmgm(y2,y1,a)*a
      vdyf(x1,x2,a,r)=(abs(a)-a**2/r)*(abs(x2)-abs(x1))
     1                               /(abs(x2)+abs(x1)+ep)
      vcorr(a,b,y1,y2,r)=-0.125*a*b*y1/(y2*r)
      vcor31(a,x0,x1,x2,x3,r)= -(a -3.*abs(a)*a/r+2.*a**3/r**2)/3.
     1                         *(abs(x0)+abs(x3)-abs(x1)-abs(x2))
     2                         /(abs(x0)+abs(x3)+abs(x1)+abs(x2)+ep)
      vcor32(a,b,y1,y2,r)=0.25*b/r*(abs(a)-2.*a**2/r)*y1/y2
      vdiv1(a1,a2,a3,r)=0.25*a2*(a3-a1)/r
      vdiv2(a,b1,b2,b3,b4,r)=0.25*a*(b1+b2-b3-b4)/r
      pp(y)= amax1(0.,y)
      pn(y)=-amin1(0.,y)
 
      iord=iord0
      if(isor.eq.3) iord=max0(iord,3)
      if(liner.eq.1) iord=1

      do j=1,n2-1
        do i=1,n1
          v1(i,j) = u1(i,j)
        end do 
      end do 
      do i=1,n1-1
        do j=1,n2
          v2(i,j) = u2(i,j)
        end do 
      enddo

      if(nonos.eq.1) then
      do j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do i=1,n1m
      im=(i-1+(n1-i)/n1m*(n1-2))
      ip=(i+1    -i /n1m*(n1-2))
      mx(i,j)=amax1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
      mn(i,j)=amin1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp))
      end do
      end do
      endif
 
                         do 3 k=1,iord
 
      do 331 j=1,n2-1
      do 331 i=2,n1-1
  331 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      do j=1,n2-1
      f1(1 ,j)=f1(n1-1,j)
      f1(n1,j)=f1(2,j)
      enddo
      do 332 j=2,n2-1
      do 332 i=1,n1-1
  332 f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
      if (iflg.eq.6) then
        do i=1,n1m
          f2(i, 1)=donor(x(i,  1),x(i,  1),v2(i, 1))
          f2(i,n2)=donor(x(i,n2m),x(i,n2m),v2(i,n2))
        end do
      else
        do i=1,n1m
          f2(i, 1)=-f2(i,  2)
          f2(i,n2)=-f2(i,n2m)
        end do
      end if
  
      do 333 j=1,n2-1
      do 333 i=1,n1-1
  333 x(i,j)=x(i,j)-(f1(i+1,j)-f1(i,j)+f2(i,j+1)-f2(i,j))/h(i,j)
 
      if(k.eq.iord) go to 6

      do 49 j=1,n2-1
      do 49 i=1,n1
      f1(i,j)=v1(i,j)
   49 v1(i,j)=0.
      do 50 j=1,n2
      do 50 i=1,n1-1
      f2(i,j)=v2(i,j)
   50 v2(i,j)=0.
      do 51 j=2,n2-2
      do 51 i=2,n1-1
   51 v1(i,j)=vdyf(x(i-1,j),x(i,j),f1(i,j),.5*(h(i-1,j)+h(i,j)))
     *       +vcorr(f1(i,j), f2(i-1,j)+f2(i-1,j+1)+f2(i,j+1)+f2(i,j),
     *   abs(x(i-1,j+1))+abs(x(i,j+1))-abs(x(i-1,j-1))-abs(x(i,j-1)),
     *   abs(x(i-1,j+1))+abs(x(i,j+1))+abs(x(i-1,j-1))+abs(x(i,j-1))+ep,
     *                 .5*(h(i-1,j)+h(i,j)))
      if(idiv.eq.1) then
      do 511 j=2,n2-2
      do 511 i=2,n1-1
  511 v1(i,j)=v1(i,j)
     *    -vdiv1(f1(i-1,j),f1(i,j),f1(i+1,j),.5*(h(i-1,j)+h(i,j)))
     *    -vdiv2(f1(i,j),f2(i-1,j+1),f2(i,j+1),f2(i-1,j),f2(i,j),
     *                 .5*(h(i-1,j)+h(i,j)))
      endif
      do 52 j=2,n2-1
      do 52 i=2,n1-2
   52 v2(i,j)=vdyf(x(i,j-1),x(i,j),f2(i,j),.5*(h(i,j-1)+h(i,j)))
     *       +vcorr(f2(i,j), f1(i,j-1)+f1(i,j)+f1(i+1,j)+f1(i+1,j-1),
     *   abs(x(i+1,j-1))+abs(x(i+1,j))-abs(x(i-1,j-1))-abs(x(i-1,j)),
     *   abs(x(i+1,j-1))+abs(x(i+1,j))+abs(x(i-1,j-1))+abs(x(i-1,j))+ep,
     *                 .5*(h(i,j-1)+h(i,j)))
      i0=n1-2
      do j=2,n2-1
      v2(1,j)=vdyf(x(1,j-1),x(1,j),f2(1,j),.5*(h(1,j-1)+h(1,j)))
     *       +vcorr(f2(1,j), f1(1,j-1)+f1(1,j)+f1(2,j)+f1(2,j-1),
     *   abs(x(2,j-1))+abs(x(2,j))-abs(x(i0,j-1))-abs(x(i0,j)),
     *   abs(x(2,j-1))+abs(x(2,j))+abs(x(i0,j-1))+abs(x(i0,j))+ep,
     *                 .5*(h(1,j-1)+h(1,j)))
      v2(n1-1,j)=v2(1,j)
      enddo

      if(idiv.eq.1) then
      do 521 j=2,n2-1
      do 521 i=1,n1-1
  521 v2(i,j)=v2(i,j)
     *    -vdiv1(f2(i,j-1),f2(i,j),f2(i,j+1),.5*(h(i,j-1)+h(i,j)))
     *    -vdiv2(f2(i,j),f1(i+1,j),f1(i+1,j-1),f1(i,j-1),f1(i,j),
     *                 .5*(h(i,j-1)+h(i,j)))
      endif
      if(isor.eq.3) then
      do 61 j=2,n2-2
      do 61 i=3,n1-2
   61 v1(i,j)=v1(i,j)     +vcor31(f1(i,j),
     1        x(i-2,j),x(i-1,j),x(i,j),x(i+1,j),.5*(h(i-1,j)+h(i,j)))
      do j=2,n2-2
      v1(2,j)=v1(2,j)     +vcor31(f1(2,j),
     1        x(n1-2,j),x(1,j),x(2,j),x(3,j),.5*(h(1,j)+h(2,j)))
      v1(n1-1,j)=v1(n1-1,j) +vcor31(f1(n1-1,j),x(n1-3,j),x(n1-2,j),
     1                  x(n1-1,j),x(2,j),.5*(h(n1-2,j)+h(n1-1,j)))
      enddo
      do 62 j=2,n2-2
      do 62 i=2,n1-1
   62 v1(i,j)=v1(i,j)
     1 +vcor32(f1(i,j),f2(i-1,j)+f2(i-1,j+1)+f2(i,j+1)+f2(i,j),
     *   abs(x(i,j+1))-abs(x(i,j-1))-abs(x(i-1,j+1))+abs(x(i-1,j-1)),
     *   abs(x(i,j+1))+abs(x(i,j-1))+abs(x(i-1,j+1))+abs(x(i-1,j-1))+ep,
     *                   .5*(h(i-1,j)+h(i,j)))
      do 63 j=3,n2-2
      do 63 i=1,n1-1
   63 v2(i,j)=v2(i,j)     +vcor31(f2(i,j),
     1        x(i,j-2),x(i,j-1),x(i,j),x(i,j+1),.5*(h(i,j-1)+h(i,j)))
      do 64 j=3,n2-2
      do 64 i=2,n1-2
   64 v2(i,j)=v2(i,j)
     1 +vcor32(f2(i,j),f1(i,j-1)+f1(i+1,j-1)+f1(i+1,j)+f1(i,j),
     *   abs(x(i+1,j))-abs(x(i-1,j))-abs(x(i+1,j-1))+abs(x(i-1,j-1)),
     *   abs(x(i+1,j))+abs(x(i-1,j))+abs(x(i+1,j-1))+abs(x(i-1,j-1))+ep,
     *                   .5*(h(i,j-1)+h(i,j)))
      do 641 j=3,n2-2
      v2(1,j)=v2(1,j)
     1 +vcor32(f2(1,j),f1(1,j-1)+f1(2,j-1)+f1(2,j)+f1(1,j),
     *   abs(x(2,j))-abs(x(n1-2,j))-abs(x(2,j-1))+abs(x(n1-2,j-1)),
     *   abs(x(2,j))+abs(x(n1-2,j))+abs(x(2,j-1))+abs(x(n1-2,j-1))+ep,
     *                   .5*(h(1,j-1)+h(1,j)))
  641 v2(n1-1,j)=v2(1,j)
      endif
 
        do j=1,n2m
          v1( 1,j)=v1(n1m,j)
          v1(n1,j)=v1(  2,j)
        end do

      if (iflg.ne.6) then
        do i=1,n1m
          v2(i, 1)=-v2(i,  2)
          v2(i,n2)=-v2(i,n2m)
        end do
      end if

                  if(nonos.eq.1) then
c                 non-osscilatory option
      do 401 j=1,n2m
      jm=max0(j-1,1  )
      jp=min0(j+1,n2m)
      do 401 i=1,n1m
      im=(i-1+(n1-i)/n1m*(n1-2))
      ip=(i+1    -i /n1m*(n1-2))
      mx(i,j)=amax1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mx(i,j))
  401 mn(i,j)=amin1(x(im,j),x(i,j),x(ip,j),x(i,jm),x(i,jp),mn(i,j))

      do 402 j=1,n2m 
      do 4021 i=2,n1-1
 4021 f1(i,j)=donor(x(i-1,j),x(i,j),v1(i,j))
      f1(1 ,j)=f1(n1m,j)
      f1(n1,j)=f1(2  ,j)
  402 continue
     
      do 403 i=1,n1m
      do 4031 j=2,n2m
 4031 f2(i,j)=donor(x(i,j-1),x(i,j),v2(i,j))
      if(iflg.ne.6) then
      f2(i, 1)=-f2(i,  2)
      f2(i,n2)=-f2(i,n2m)
      else
      f2(i, 1)=0.
      f2(i,n2)=0.
      endif
  403 continue

      do 404 j=1,n2m   
      do 404 i=1,n1m
      cp(i,j)=(mx(i,j)-x(i,j))*h(i,j)/
     1(pn(f1(i+1,j))+pp(f1(i,j))+pn(f2(i,j+1))+pp(f2(i,j))+ep)
      cn(i,j)=(x(i,j)-mn(i,j))*h(i,j)/
     1(pp(f1(i+1,j))+pn(f1(i,j))+pp(f2(i,j+1))+pn(f2(i,j))+ep)
  404 continue
      do 405 j=1,n2m 
      do 4051 i=2,n1m 
 4051 v1(i,j)=pp(v1(i,j))*
     1  ( amin1(1.,cp(i,j),cn(i-1,j))*pp(sign(1., x(i-1,j)))
     1   +amin1(1.,cp(i-1,j),cn(i,j))*pp(sign(1.,-x(i-1,j))) )
     2       -pn(v1(i,j))*
     2  ( amin1(1.,cp(i-1,j),cn(i,j))*pp(sign(1., x(i ,j )))
     2   +amin1(1.,cp(i,j),cn(i-1,j))*pp(sign(1.,-x(i ,j ))) )
      v1( 1,j)=v1(n1m,j)
      v1(n1,j)=v1( 2 ,j)
  405 continue

      do 406 i=1,n1m 
      do 406 j=2,n2m 
  406 v2(i,j)=pp(v2(i,j))*
     1  ( amin1(1.,cp(i,j),cn(i,j-1))*pp(sign(1., x(i,j-1)))
     1   +amin1(1.,cp(i,j-1),cn(i,j))*pp(sign(1.,-x(i,j-1))) )
     1       -pn(v2(i,j))*
     2  ( amin1(1.,cp(i,j-1),cn(i,j))*pp(sign(1., x(i ,j )))
     2   +amin1(1.,cp(i,j),cn(i,j-1))*pp(sign(1.,-x(i ,j ))) )
                  endif
    3                      continue
    6 continue
      return
      end   


      subroutine gcrk_1(p,pfx,pfz,u,w,n1,n3,itr,eps0)
      real p(*),pfx(*),pfz(*),u(*),w(*)
      include 'param.grid'
      parameter(n=nx,l=nz)
      parameter(nn=n*l,nl=n*l)
      common// r(nn),qr(nn),ar(nn)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /strtch/ height(nz),gac(nz)
      common/itero/ niter,nitsm,icount,eer,eem
      parameter (lord=3)
      dimension x(nn,lord),ax(nn,lord),ax2(lord),axar(lord),del(lord)
      dimension rho2d(nx,nz)
convergence test modes **************************************************
      logical ctest                                                     *
      data ctest/.false./                                               *
c     data ctest/.true./                                                *
      parameter (nplt=100)                                              *
      dimension err(0:nplt),xitr(0:nplt)                                *
      if(ctest) then                                                    *
      itr=6000/lord                                                     *
      ner=60                                                            *
      snorm=1./float(n*l)                                               *
      eps0=1.e-15                                                       *
      endif                                                             *
convergence test modes **************************************************
 
      do k=1,l
      do i=1,n
      rho2d(i,k)=rho0(k)*gac(k)
      enddo
      enddo

      eps=eps0*dti
      epa=1.e-30
      nlc=0
      itmn=5
cc iprc 0-no preconditioning, 1-with precon
      iprc=1

      call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,0)

      do k=1,nl
        r(k)=0.
       ar(k)=0.
       qr(k)=0.
      enddo
      do i=1,lord
       do k=1,nl
         x(k,i)=0.
        ax(k,i)=0.
       enddo
      enddo
      call prforc_1(p,pfx,pfz,u,w,n1,n3)
       call rhsdiv_1(pfx,pfz,rho2d,r,n1,n3,-1)
        call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1)
      eer0=0.
      eem0=-1.e15
      rl20=0.
      do k=1,nl
      eer0=eer0+qr(k)**2
      eem0=amax1(eem0,abs(qr(k)))
      rl20=rl20+r(k)**2
      enddo
      eer0=amax1(eer0,epa)
      eem0=amax1(eem0,epa)
      rl20=amax1(rl20,epa)
convergence test modes **************************************************
      if(ctest) then                                                    *
      do ier=0,nplt                                                     *
      err(ier)=eps                                                      *
      enddo                                                             *
      eer=-1.e15                                                        *
      do 3 k=1,nl                                                       *
    3 eer=amax1(eer,abs(r(k)))                                          *
      err(0)=eer                                                        *
      print 300,  err(0)                                                *
  300 format(4x,e11.4,' residual error at it=1')                        *
      endif                                                             *
convergence test modes **************************************************
       do k=1,nl
        x(k,1)=qr(k)
       enddo
      call laplc_1(x(1,1),ax(1,1),pfx,pfz,n1,n3)

      do 100 it=1,itr
       do i=1,lord
        ax2(i)=0.
        rax=0.
         do k=1,nl
          rax=rax+r(k)*ax(k,i)
          ax2(i)=ax2(i)+ax(k,i)*ax(k,i)
         enddo
        ax2(i)=amax1(epa,ax2(i))
        beta=-rax/ax2(i)
        dvmx=-1.e15
        rl2=0.
         do k=1,nl
          p(k)=p(k)+beta* x(k,i)
          r(k)=r(k)+beta*ax(k,i)
          dvmx=amax1(dvmx,abs(r(k)))
          rl2=rl2+r(k)*r(k)
         enddo
       if(dvmx.le.eps.and.it.ge.itmn) go to 200
       if(rl2.ge.rl20.and..not.ctest) go to 200
          rl20=amax1(rl2,epa)
       call precon_1(r,qr,ar,pfx,pfz,rho2d,iprc,1)
       call laplc_1(qr,ar,pfx,pfz,n1,n3)
        nlc=nlc+1
         do ii=1,i
          axar(ii)=0.
           do k=1,nl
            axar(ii)=axar(ii)+ax(k,ii)*ar(k)
           enddo
          del(ii)=-axar(ii)/ax2(ii)
         enddo
        if(i.lt.lord) then
          do k=1,nl
            x(k,i+1)=qr(k)
           ax(k,i+1)=ar(k)
          enddo
           do ii=1,i
            do k=1,nl
              x(k,i+1)= x(k,i+1)+del(ii)* x(k,ii)
             ax(k,i+1)=ax(k,i+1)+del(ii)*ax(k,ii)
            enddo
           enddo
        else
          do k=1,nl
            x(k,1)=qr(k)+del(1)* x(k,1)
           ax(k,1)=ar(k)+del(1)*ax(k,1)
          enddo
           do ii=2,i
            do k=1,nl
              x(k,1 )= x(k,1)+del(ii)* x(k,ii)
              x(k,ii)=0.
             ax(k,1 )=ax(k,1)+del(ii)*ax(k,ii)
             ax(k,ii)=0.
            enddo
           enddo
        endif
convergence test modes **************************************************
      if(ctest) then                                                    *
      if(nlc/ner*ner.eq.nlc) then                                       *
      ier=nlc/ner                                                       *
      eer=-1.e15                                                        *
      do 50 k=1,nl                                                      *
   50 eer=amax1(eer,abs(r(k)))                                          *
      err(ier)=eer                                                      *
      endif                                                             *
      endif                                                             *
convergence test modes **************************************************
       enddo
  100 continue
  200 continue
      eer=0.
      eem=-1.e15
      do k=1,nl
      eer=eer+qr(k)**2
      eem=amax1(eem,abs(qr(k)))
      enddo
      eer=eer/eer0
      eem=eem/eem0
      niter=nlc
      nitsm=nitsm+niter
      icount=icount+1

convergence test modes **************************************************
      if(ctest) then                                                    *
      print 301, (err(ier),ier=1,nplt,1)                                *
  301 format(4x,5e11.4)                                                 *
      do 400 ier=0,nplt                                                 *
      xitr(ier)=ier*ner                                                 *
  400 err(ier)=alog10(err(ier)*dt )                                     *
      plmx=float(itr*lord)                                              *
      call set(.1,.9,.1,.9,0.,plmx,-10.,0.,1)                           *
      call labmod('(i4)','(f5.0)',4,4,2,2,20,20,0)                      *
      call periml(4,10,5,2)                                             *
      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)                         *
      call curved(xitr,err,nplt+1)                                      *
      i1=int(102.4+409.6)                                               *
      call wtstr(cpux(i1),cpuy(50),'niter',3,0,0)                       *
      call wtstr(cpux(17),cpuy(i1),'log e',3,90,0)                      *
      call frame                                                        *
      endif                                                             *
convergence test modes **************************************************
      return
      end

      subroutine precon_1(rhs,p,r,c11,c33,d,iflg,jfl)
      include 'param.grid'
      parameter(n=nx,l=nz,nl=nx*nz)
      dimension rhs(n,l),p(n,l),r(n,l),
     .          c11(n,l),c33(n,l),d(n,l)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
compute available storage 
      dimension e(n,0:l-1),f(n,0:l-1)
     .       ,px(n,l),dgh(n,l),po(n,l)

      data beta/-1.e15/
      data itr,line/2,1/

      if(iflg.eq.0) then
       do i=1,nl
        p(i,1)=rhs(i,1)
       enddo
      return
      endif
    
      omg=.7
      oms=1.-omg
      dxi2=0.25*dxi*dxi
      dzi2=0.25*dzi*dzi
      do i=1,nl
       c33(i,1)=d(i,1)*dzi2
       c11(i,1)=d(i,1)*dxi2
       dgh(i,1)=0.
        po(i,1)=0.
         p(i,1)=0.
         r(i,1)=0.
      enddo
      if(line.eq.1) then
       do k=1,l
         do i=2,n-1
          dgh(i,k)=c11(i+1,k)+c11(i-1,k)
         enddo
          dgh(1,k)=c11(2,k)+c11(n-1,k)
          dgh(n,k)=c11(2,k)+c11(n-1,k)
       enddo
      endif

      if(jfl.eq.0) then
      if(line.eq.0) then
      beta=-1.e15
      do i=1,nl
      beta=amax1(beta,abs(c11(i,1))/d(i,1))
      enddo
      beta=0.5/beta
      else
      beta=1.
      endif
      return
      endif
      beti=1./beta*(1-line)

      do 100 it=1,itr
      do i=1,nl
       r(i,1)=r(i,1)+d(i,1)*(beti*p(i,1)-rhs(i,1))
     .                  +dgh(i,1)*p(i,1)
      enddo
      do i=1,n
       e(i,0)=1.
       f(i,0)=0.
       dn=d(i,1)*beti+2.*c33(i,2)+dgh(i,1)
       e(i,1)=2.*c33(i,2)/dn
       f(i,1)=     r(i,1)/dn
      enddo
      do k=2,l-1
       do i=1,n
        dn=c33(i,k+1)+c33(i,k-1)*(1.-e(i,k-2))+d(i,k)*beti
     .                                + dgh(i,k)
        e(i,k)=                      c33(i,k+1)/dn
        f(i,k)=(c33(i,k-1)*f(i,k-2)+r(i,k))/dn
       enddo
      enddo
       do i=1,n
        dn=d(i,l)*beti+2.*(1.-e(i,l-2))*c33(i,l-1)
     .                                + dgh(i,l)
        p(i,l)=(r(i,l)+2.*f(i,l-2)*c33(i,l-1))/dn
        p(i,l-1)=f(i,l-1)/(1.-e(i,l-1))
       enddo
      do k=l-2,1,-1
       do i=1,n
        p(i,k)=e(i,k)*p(i,k+2)+f(i,k)
       enddo
      enddo


      if(line.eq.1) then
       do i=1,nl
        p(i,1)=oms*po(i,1)+omg*p(i,1)
       po(i,1)=     p(i,1)
       enddo
      endif

      if(it.eq.itr) go to 101
      do k=1,l
      do i=2,n-1
      px(i,k)=c11(i,k)*(p(i+1,k)-p(i-1,k))
      enddo
      px(1,k)=c11(1,k)*(p(2,k)-p(n-1,k))
      px(n,k)=c11(n,k)*(p(2,k)-p(n-1,k))
      enddo

      do k=1,l
      do 91 i=2,n-1
   91 r(i,k)=px(i+1,k)-px(i-1,k)
      r(1,k)=(px(2,k)-px(n-1,k))
      r(n,k)=(px(2,k)-px(n-1,k))
      enddo

  100 continue
  101 continue

      return
      end


      subroutine laplc_1(p,r,u,w,n1,l1)
      dimension p(n1,l1),r(n1,l1),u(n1,l1),w(n1,l1)
      include 'param.grid'
      parameter(n=nx,l=nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
compute available storage in //
      parameter (nl=n*l)
      dimension px(n,l),pz(n,l)

      dxil=.5*dxi
      dzil=.5*dzi

compute pressure derivatives everywhere
      do 18 k=1,l
      do 1 i=2,n-1
    1 px(i,k)=     dxil*(p(i+1,k)-p(i-1,k))
      px(1,k)=dxil*(p(2,k)-p(n-1,k))
      px(n,k)=dxil*(p(2,k)-p(n-1,k))
   18 continue
      do 38 i=1,n
      do 3 k=2,l-1
    3 pz(i,k)=dzil*(p(i,k+1)-p(i,k-1))
      pz(i,1)= dzi*(p(i,2)-p(i,1))
   38 pz(i,l)= dzi*(p(i,l)-p(i,l-1))

compute interior pressure forces
      do 21 i=1,n
      do 10 k=2,l-1
      u(i,k)=-px(i,k)
   10 w(i,k)=-pz(i,k)
      w(i,1)=0.
      w(i,l)=0.
      u(i,1)=-px(i,1)
   21 u(i,l)=-px(i,l)

      do 99 k=1,l
      coef=rho0(k)*gac(k)
      do 99 i=1,n
      u(i,k)=coef*u(i,k)
   99 w(i,k)=coef*w(i,k)

compute laplacian
      do 911 k=1,l
      do 91  i=2,n-1
   91 r(i,k)=dxil*(u(i+1,k)-u(i-1,k))
      r(1,k)=dxil*(u(2,k)-u(n-1,k))
      r(n,k)=dxil*(u(2,k)-u(n-1,k))
  911 continue
      do 931 i=1,n
      do 93  k=2,l-1
   93 r(i,k)=r(i,k)+dzil*(w(i,k+1)-w(i,k-1))
      r(i,1)=r(i,1)+dzi *(w(i,2)-w(i,1)) 
      r(i,l)=r(i,l)+dzi *(w(i,l)-w(i,l-1))
  931 continue
      do 94 i=1,n
      do 94 k=1,l
   94 r(i,k)=-r(i,k)/(rho0(k)*gac(k))

      return
      end

      subroutine prforc_1(p,pfx,pfz,u,w,n1,n3)
      dimension p(n1,n3),u(n1,n3),w(n1,n3),pfx(n1,n3),pfz(n1,n3)
      include 'param.grid'
      parameter(n=nx,l=nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      parameter (nl=n*l)
      dimension px(n,l),pz(n,l)

      dxil=.5*dxi
      dzil=.5*dzi

compute pressure derivatives everywhere
      do 18 k=1,l
      do 1 i=2,n-1
    1 px(i,k)=     dxil*(p(i+1,k)-p(i-1,k))
      px(1,k)=dxil*(p(2,k)-p(n-1,k))
      px(n,k)=dxil*(p(2,k)-p(n-1,k))
   18 continue
      do 38 i=1,n
      do 3 k=2,l-1
    3 pz(i,k)=dzil*(p(i,k+1)-p(i,k-1))
      pz(i,1)= dzi*(p(i,2)-p(i,1))
   38 pz(i,l)= dzi*(p(i,l)-p(i,l-1))

compute interior pressure forces
      do 21 i=1,n
      do 10 k=2,l-1
      pfx(i,k)=u(i,k)-px(i,k)
   10 pfz(i,k)=w(i,k)-pz(i,k)
      pfz(i,1)=0.
      pfz(i,l)=0.
      pfx(i,1)=u(i,1)-px(i,1)
   21 pfx(i,l)=u(i,l)-px(i,l)

      return
      end

      subroutine rhsdiv_1(u,w,d,r,n1,l1,iflg)
      dimension u(n1,l1),w(n1,l1),d(n1,l1),r(n1,l1)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

      n=n1
      l=l1
      nl=n*l

      do 200 i=1,nl
  200 r(i,1)=0.

      dxil=.5*dxi
      dzil=.5*dzi

      do 1 i=2,n-1
      do 1 j=1,l
    1 r(i,j)=dxil*(u(i+1,j)*d(i+1,j)-u(i-1,j)*d(i-1,j))
      do 11 j=1,l
      r(1,j)=dxil*(u(2,j)*d(2,j)-u(n-1,j)*d(n-1,j))
      r(n,j)=dxil*(u(2,j)*d(2,j)-u(n-1,j)*d(n-1,j))
   11 continue
      do 3 k=2,l-1
      do 3 i=1,n
    3 r(i,k)=r(i,k)
     3        +dzil*(w(i,k+1)*d(i,k+1)-w(i,k-1)*d(i,k-1))
      do 13 i=1,n
      r(i,1)=r(i,1)+dzi*(w(i,2)*d(i,2)-w(i,1)*d(i,1))
   13 r(i,l)=r(i,l)+dzi*(w(i,l)*d(i,l)-w(i,l-1)*d(i,l-1))

      if(iflg.ne.0) then
      do 4 i=1,nl
    4 r(i,1)=iflg*r(i,1)/d(i,1)
      endif

      return
      end

       subroutine minmax(a,n,an,ax)
       dimension a(n)
       an= 1.e15
       ax=-1.e15
       do i=1,n
       an=amin1(a(i),an)
       ax=amax1(a(i),ax)
       enddo
       return
       end

       subroutine minmax11(a,nx,nz,an,ax,ii,kk)
       dimension a(nx,nz)
       ii=-100
       kk=-100
       an= 1.e15
       ax=-1.e15
       do i=1,nx
       do k=1,nz
       an=amin1(a(i,k),an)
       ax=amax1(a(i,k),ax)
       if(ax.eq.a(i,k)) then
       ii=i
       kk=k
       endif
       enddo
       enddo
       return
       end

      subroutine diagno_1(ux,uz,th,nx,nz,scr1,scr2,rho)
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),scr1(nx,nz),scr2(nx,nz)
      dimension rho(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
cc pressure solver diagnostics
      common/itero/ niter,nitsm,icount,eer,eem

      nxz=nx*nz
      do i=1,nx
      do k=1,nz
      scr2(i,k)=rho(i,k)
      enddo
      enddo

      print 200, time
 200  format(1x,' ****** analysis for time (min): ',f8.2)

      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max ux: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max uz: ',2e12.4)

      cour=0.
      do i=1,nxz
      cour=amax1(cour,abs(ux(i,1))*dt/dx+abs(uz(i,1))*dt/dz)
      enddo
      print 302,cour
 302  format(1x,' --> max courant: ',e12.4)

      call minmax(th,nxz,amn,amx)
      print 203,amn,amx
 203  format(1x,' --> min, max th: ',2e12.4)
      call rhsdiv_1(ux,uz,scr2,scr1,nx,nz,1)

      call minmax(scr1,nxz,amn,amx)
      print 204,amn,amx
 204  format(1x,' --> min, max div: ',2e12.4)

      nitav=nitsm/max0(icount,1)
      print 205, eer,eem,niter,nitav
  205 format(1x,'            eer,eem:',2e11.4/
     .       1x,'       niter, nitav:',2i4)

       if(cour.gt.1.) then
       call clsgks
       stop 'courant'
       endif

       return
       end

      subroutine diagno_2(ux,uz,nx,nz,scr1,scr2,rho)
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),scr1(nx,nz),scr2(nx,nz)
      dimension rho(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
cc pressure solver diagnostics
      common/itero/ niter,nitsm,icount,eer,eem

      nxz=nx*nz

      call minmax(ux,nxz,amn,amx)
      print 201,amn,amx
 201  format(1x,' --> min, max qv: ',2e12.4)

      call minmax(uz,nxz,amn,amx)
      print 202,amn,amx
 202  format(1x,' --> min, max qc: ',2e12.4)

       return
       end

      subroutine diagno_3(th,qv,rh)
      include 'param.grid'
      dimension qv(nx,nz),rh(nx,nz),th(nx,nz)
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)

cc rh:
      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

      do k=1,nz
      thetme=th_e(k)/tm_e(k)
      pre=1.e5*thetme**e
      do i=1,nx
      tt=th(i,k)/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      qvsw=a * esw /(pre-esw)
      rh(i,k)=(qv(i,k)/qvsw  - 1.)*100.
      enddo
      enddo

      call minmax11(rh,nx,nz,amn,amx,i,k)
      print 202,amn,amx
 202  format(1x,' --> min, max S(%): ',2e12.4)
      print*,' for (i,k): ',i,k

      return
      end

      subroutine diagno_4(npx,xp,zp,qcp,plic,radp)
      include 'param.grid'
      include 'param.sds'
      dimension qc(nx,nz),con(nx,nz),rad(nx,nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      dimension xp(np),zp(np)
      dimension qcp(np),plic(np),radp(np)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi

       do i=1,nxz
       rad(i,1)=0.
       con(i,1)=0.
       enddo

       do ip=1,npx
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
          con1=plic(ip)
          con(i,k)=con(i,k) + con1
          rad(i,k)=rad(i,k) + radp(ip)*con1
       enddo

      do i=1,nx
      do k=1,nz
      if(con(i,k).gt.1.) then
      rad(i,k)=rad(i,k)/con(i,k)
      else
      rad(i,k)=0.
      endif
      con(i,k)=con(i,k)*1.e-6  ! to cm**-3
      rad(i,k)=rad(i,k)*1.e6   ! to micron
      enddo
      enddo

      call minmax11(con,nx,nz,amn,amx,i,k)
      print 202,amn,amx
 202  format(1x,' --> min, max conc (cm**3): ',2e12.4)
      print*,' for (i,k): ',i,k

      call minmax11(rad,nx,nz,amn,amx,i,k)
      print 203,amn,amx
 203  format(1x,' --> min, max rad (micron): ',2e12.4)
      print*,' for (i,k): ',i,k

      return
      end

      subroutine plot_1(ux,uz,th)
      include 'param.grid'
      dimension ux(nx,nz),uz(nx,nz),th(nx,nz),xx(nx),zz(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)
      dimension pl(nx,nz),w1(nx,nz),w2(nx,nz)

      character*50 lhead

      xl=float(nx-1)*dx/1.e3
      zl=float(nz-1)*dz/1.e3

      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo
      rat=zl/xl

cc theta
      do i=1,nx
      do k=1,nz
      pl(i,k)=th(i,k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=288.
      amx=352.
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)                           
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,290.,590.,2.,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,951) time
  951 format('   theta (K)      time (min): ',f6.1)
c      CALL plchmq(.525,0.93, lhead(1:50), 0.016,0.,0)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame

cc ux
      do i=1,nx
      do k=1,nz
      pl(i,k)=ux(i,k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=-5.
      amx=5.
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)                           
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,-38.,38.,0.5,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,952) time
  952 format('     ux (m/s)      time (min): ',f6.1)
c      CALL plchmq(.525,0.93, lhead(1:50), 0.016,0.,0)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame
cc uz
      do i=1,nx
      do k=1,nz
      pl(i,k)=uz(i,k)*gac(k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=-1.
      amx=1.
      sn=1.e4
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call prof(uz,height,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)                           
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,-28.,28.,0.5,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,953) time
  953 format('     uz (m/s)      time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame

      return
      end

      subroutine plot_2(th,qv,qc)
      include 'param.grid'
      dimension qv(nx,nz),qc(nx,nz),xx(nx),zz(nz)
      dimension sc(nx,nz),th(nx,nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      character*50 lhead
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      dimension pl(nx,nz),w1(nx,nz),w2(nx,nz)


      xl=float(nx-1)*dx/1.e3
      zl=float(nz-1)*dz/1.e3

      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo
      rat=zl/xl

cc qv
      do i=1,nx
      do k=1,nz
      pl(i,k)=qv(i,k)*1.e3
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=10.
      amx=18.
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)                           
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,0.5,15.,0.5,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,951) time
  951 format('   qv (g/kg)      time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame
cc qc
      do i=1,nx
      do k=1,nz
      pl(i,k)=qc(i,k)*1.e3
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=0.
      amx=.1
      sn=1.
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)                           
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,.1,30.,.1,-1,-1,-682)
      call cpcnrc(pl,nx,nx,nz,.01,.011,.01,-1,-1,682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,952) time
  952 format('     qc (g/kg)      time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame

cc rh:
      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

      do k=1,nz
      thetme=th_e(k)/tm_e(k)
      pre=1.e5*thetme**e
      do i=1,nx
      tt=th(i,k)/thetme
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      qvsw=a * esw /(pre-esw)
      sc(i,k)=(qv(i,k)/qvsw-1.) * 100.
      enddo
      enddo
      do i=1,nx
      do k=1,nz
      pl(i,k)=sc(i,k)
      enddo
      enddo
      call gridint(pl,nx,nz,xx,height,nx,nz,xx,zz,w1,w2)
      amn=0.0
      amx=100.0
      sn=1
      call setusv('LW',2000)
c      call prof(pl,zz,nx,nz,sn,amn,amx,zl)
c      call set(.15,.9,.15,.9, 0.,xl,0.,zl,1)
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)                           
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
c      call cpcnrc(pl,nx,nx,nz,.0,120.,10.,-1,-1,-682)
      call cpcnrc(pl,nx,nx,nz,.0,20.,.05,-1,-1,-682)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,993) time
c  993 format('      rh (%)        time (min): ',f6.1)
  993 format(' positive S (%)        time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.08, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.425, 'z (km)', 0.016,90.,0)
      call frame

      return
      end

      subroutine plot_3(qc,xp,zp,qcp,plic,radp,npx)
      include 'param.grid'
      include 'param.sds'
      dimension xp(np),zp(np)
      dimension qcp(np),plic(np),radp(np)
      dimension qc(nx,nz),con(nx,nz),rad(nx,nz),xx(nx),zz(nz)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      character*50 lhead
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats
      common /strtch/ height(nz),gac(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz),qc_e(nz)
      dimension pl(nx,nz)


      xl=float(nx-1)*dx/1.e3
      zl=float(nz-1)*dz/1.e3

      do i=1,nx
      xx(i)=float(i-1)*dx
      enddo
      do k=1,nz
      zz(k)=float(k-1)*dz
      enddo
      rat=zl/xl

cc qc
      do i=1,nx
      do k=1,nz
      pl(i,k)=qc(i,k)*1.e3
      enddo
      enddo
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,.005,4.,.2,-1,-1,-682)
c      call cpcnrc(pl,nx,nx,nz,.01,.011,.01,-1,-1,682)
cc SDs:
       do ip=1,npx
       call points(xp(ip)*1.e-3,zp(ip)*1.e-3,1,-1,0)
c       call points(xp(ip)*1.e-3,zp(ip)*1.e-3,1,-2,0)
       enddo

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,952) time
  952 format(' qc (g/kg) + SDs time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.02, 'z (km)', 0.016,0.,0)
      CALL plchhq(.02,0.525, 'z (km)', 0.016,90.,0)
      call frame

cc concentration and mean radius:
       do i=1,nxz
       rad(i,1)=0.
       con(i,1)=0.
       enddo

       do ip=1,npx
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
          con1=plic(ip)
          con(i,k)=con(i,k) + con1
          rad(i,k)=rad(i,k) + radp(ip)*con1
       enddo

      do i=1,nx
      do k=1,nz
      if(con(i,k).gt.1.) then
      rad(i,k)=rad(i,k)/con(i,k)
      else
      rad(i,k)=0.
      endif
      con(i,k)=con(i,k)*1.e-6  ! to cm**-3
      rad(i,k)=rad(i,k)*1.e6   ! to micron
      enddo
      enddo

      do i=1,nx
      do k=1,nz
      pl(i,k)=con(i,k)
      enddo
      enddo
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,.01,10000.,10.,-1,-1,-682)
      call cpcnrc(pl,nx,nx,nz,.1,.11,.1,-1,-1,682)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,972) time
  972 format('  nc (1/cm**3)      time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.02, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.525, 'z (km)', 0.016,90.,0)
      call frame

cc conc over 100
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,100.01,2000.,5.,-1,-1,-682)
      call cpcnrc(pl,nx,nx,nz,.1,.11,.1,-1,-1,682)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,976) time
  976 format(' nc100 (1/cm**3)   time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.02, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.525, 'z (km)', 0.016,90.,0)
      call frame

      do i=1,nx
      do k=1,nz
      pl(i,k)=rad(i,k)
      enddo
      enddo
      call set(.15,.9,.15,.15+rat*(.9-.15), 0.,xl,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(2,9,2,6)
      call cpsetc('ILT',' ')
      call cpseti('LLP',0)
      call cpcnrc(pl,nx,nx,nz,.1,100.,2.,-1,-1,-682)
      call cpcnrc(pl,nx,nx,nz,.01,.011,.01,-1,-1,682)

      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      write(lhead,975) time
  975 format(' radius (micron)    time (min): ',f6.1)
      CALL plchmq(.525,0.70, lhead(1:50), 0.016,0.,0)
      CALL plchhq(.525,0.02, 'x (km)', 0.016,0.,0)
      CALL plchhq(.02,0.525, 'z (km)', 0.016,90.,0)
      call frame

      return
      end


      subroutine prof(aa,zz,nx,nz,sn,an,ax,zl)
      dimension zz(nz),aa(nx,nz)
      dimension z(400),prf(400)

      if(nz.gt.400) stop 'dim prof'

      do k=1,nz
      z(k)=zz(k)/1000.
      prf(k)=0.
      do i=1,nx-1
      prf(k)=prf(k)+aa(i,k)*sn/float(nx-1)
      enddo
      enddo
      call set(.30,.70,.15,.9, an,ax,0.,zl,1)
      call labmod('(f5.1)','(f5.1)',5,5,2,2,20,20,0)
      call periml(1,10,4,5)
      call curved(prf,z,nz)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.5,0.93, 'profile of ____', 0.016,0.,0)
      call frame

      return
      end

      subroutine moist_init

create parametrs for the model:

cc latent heats:
      data hlatv,hlats /2.53e6,2.84e6/
cc limiting temperatures
c      data tup,tdn /268.,253./   
      data tup,tdn /168.,153./   
cc reference temperature and saturated vapor pressure:
      data tt0,ee0 /273.16,611./

      common/temp_p/ tup,tdn
      common/latent/hlatv,hlats
      common/reference/ tt0,ee0

      common /const/ gg,cp,rg,rv
      data gg,cp,rg,rv   /9.81,1005.,287.,461./

       return
       end

      subroutine integz(a,b,n1,n2)
      dimension a(n1,n2),b(n1,n2)

cc z smoothing:
      do k=2,n2-1
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do i=1,n1
      b(i,1 )=0.5*(a(i, 2)+a(i,   1))
      b(i,n2)=0.5*(a(i,n2)+a(i,n2-1))
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
      end

      subroutine integxz(a,b,n1,n2)
      dimension a(n1,n2),b(n1,n2)

cc z smoothing:
      do k=2,n2-1
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do i=1,n1
      b(i,1 )=0.5*(a(i, 2)+a(i,   1))
      b(i,n2)=0.5*(a(i,n2)+a(i,n2-1))
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

cc x smoothing:
      do i=1,n1
      im=i-1
      if(im.eq.0) im=n1-1
      ip=i+1
      if(ip.eq.n1+1) ip=2
      do k=1,n2
      b(i,k)=0.25*(a(im,k)+2.*a(i,k)+a(ip,k))
      enddo
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
      end

      subroutine prof_init  
      include 'param.grid'
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)

      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)

      common/temp_p/ tup,tdn
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

c stability for the virtual potential temperature stab and relative humidity relhum
      parameter(stab=1.3e-5,relhum=0.2)

      dimension pres(nz)

      a=rg/rv
      c=hlatv/cp
      d=hlatv/rv
      e=-cp/rg

      cap=rg/cp
      capi=1./cap

cc surface data:
      zz=height(1)
      tm_e(1)=283.
      pres(1)=850.e2
      tt=tm_e(1)
      delt=(tt-tt0)/(tt*tt0)
      esw=ee0*exp(d * delt)
      qvs=a * esw /(pres(1)-esw)
      qv_e(1)=relhum*qvs

cc surface values for dry profiles:
      rho0(1)=pres(1)/rg/tm_e(1)
      th0(1)=tm_e(1)*(1.e5/pres(1))**(cap)

      th_e(1)=th0(1)*(1.+a*qv_e(1))

         k=1
         print*,'z,the,tme,qve: ',zz,th_e(k),tm_e(k),qv_e(k)

cc move upward:
       do k=2,nz
       zz=height(k)

       th_e(k)=th_e(1)*exp(stab*zz)

c predictor:
       rhob=pres(k-1)/rg/(tm_e(k-1)*(1.+a*qv_e(k-1)))
       pres(k)=pres(k-1) - gg*rhob*dz
cc iteration for T and qv:
       qv_e(k)=qv_e(k-1)
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               do iter=1,4
          tt=tm_e(k)
          delt=(tt-tt0)/(tt*tt0)
          esw=ee0*exp(d * delt)
          qvs=a * esw /(pres(k)-esw)
          qv_e(k)=relhum*qvs
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               enddo
         
c corrector:
       rhon=pres(k)/rg/(tm_e(k)*(1.+a*qv_e(k)))
       pres(k)=pres(k-1) - gg*.5*(rhob+rhon)*dz
cc iteration for T and qv:
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               do iter=1,4
          tt=tm_e(k)
          delt=(tt-tt0)/(tt*tt0)
          esw=ee0*exp(d * delt)
          qvs=a * esw /(pres(k)-esw)
          qv_e(k)=relhum*qvs
       tm_e(k)=th_e(k)*(pres(k)/1.e5)**(cap)
       tm_e(k)=tm_e(k)/(1.+a*qv_e(k))
               enddo
         
         print*,'z,the,tme,qve: ',zz,th_e(k),tm_e(k),qv_e(k)

              enddo

cc set constant-stability dry profiles:
        sum=0.
        do k=2,nz-1
          sum = sum + (th_e(k+1)-th_e(k-1))/th_e(k)
        enddo
      st=sum/(float(nz-2)*2.*dz)
       print*,'checking stability: ',stab,st
cc compute reference state vertical profiles
      cap=rg/cp
      capi=1./cap
      cs=gg/(cp*tm_e(1)*st)
      zz=0.
      print*,'z,th0,rho0: ',zz,th0(1),rho0(1)
      do k=2,nz
      zz=height(k)
      exs=exp(-st*zz)
      th0(k)=th0(1)/exs
      rho0(k)=rho0(1)*exs*(1.-cs*(1.-exs))**(capi-1.)
      print*,'z,th0,rho0: ',zz,th0(k),rho0(k)
      enddo

      return
      end

      subroutine prof_init_1
      include 'param.grid'
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /strtch/ height(nz),gac(nz)

      dimension xx(nx),zz(nz)
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)

      common/temp_p/ tup,tdn
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

      parameter(npin=5)
      dimension press(npin),temp(npin),zin(npin),relhum(npin),
     1          vap(npin),uu(npin),vv(npin),zinkm(npin)

       DATA ZIN   /  0. ,1.4  ,1.8  ,3.0,  5.0 /
       DATA TEMP  /295.5, 282.40, 280.00, 273.35, 261.54 /
       DATA VAP   /  8.5 ,  8.0 ,    3. ,    1.3  ,  0.6 /
       DATA UU /NPIN*0./
       DATA VV /NPIN*0./

cc statement functions:
      alim01(fi)=amax1(0.,amin1(1.,fi))
      comb(tm,td,tu,ad,au)=
     1  alim01((tm-td)/(tu-td))*au + alim01((tu-tm)/(tu-td))*ad

ccc pressure levels:
        do k=1,npin
        zin(k)=zin(k)*1.e3
        enddo
        press(1)=1000.
        do k=2,npin
          km=k-1
          tempk =temp(k ) * (1. + .61e-3* vap(k ))
          tempkm=temp(km) * (1. + .61e-3* vap(km))
          tav=.5*(tempk+tempkm)
          press(k)=press(km)*exp(-gg/rg/tav * (zin(k)-zin(km)))
            print*,' press: ',press(k)
         enddo

ccc height levels:
c        zin(1)=0.
c        do k=2,npin
c          km=k-1
c          tempk =temp(k )* (1.+.6e-3*vap(k ))
c          tempkm=temp(km)* (1.+.6e-3*vap(km))
c          delt=tempk-tempkm
c          if (delt.gt.1.e-4) then
c            tavi=alog(tempk/tempkm)/delt
c          else
c            tavi=1./tempk
c          endif
c          deltz=-rg/(tavi*gg) * alog(press(k)/press(km))
c          zin(k)=zin(km)+deltz
c        end do

           do k=1,npin
           print*,'z,p: ',zin(k),press(k)
           enddo

      a=rg/rv
      c=hlatv/cp
      b=hlats/rv
      d=hlatv/rv
      e=-cp/rg

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

compute environmental profiles from sounding assuming no topography:
cc surface data:
      iisn=1
      tm_e(1)=temp(iisn)
      th_e(1)=tm_e(1) * (1000./press(iisn))**(rg/cp)
      qv_e(1)=vap(iisn)*1.e-3
      ux_e(1)=uu(1)
      uy_e(1)=vv(1)
c      print*,'DETERMINED SURFACE DATA'
cc higher levels - interpolate:
c     print*,'INTERPOLATION TO HIGHER LEVELS'
      l=nz
      do k=2,l
      zz(k)=height(k)
c       print*,'k,z= ',k,zz(k)
        do kk=2,npin
          iisn=kk-1
          if(zin(kk).ge.zz(k)) go to 665
        enddo
c       print*,'INPUT SOUNDING DOES NOT GO HIGH ENOUGH. STOP.'
        stop 'SOUNDING'
 665    continue
c       print*,'iisn=',iisn
        coe2=(zz(k)-zin(iisn))/(zin(iisn+1)-zin(iisn))
        tm_e(k)=coe2*temp(iisn+1) + (1.-coe2)*temp(iisn)
        presnl=coe2*press(iisn+1) + (1.-coe2)*press(iisn)
        th_e(k)=tm_e(k) * (1000./presnl)**(rg/cp)
        qv_e(k)=(coe2*vap(iisn+1) + (1.-coe2)*vap(iisn))*1.e-3
        ux_e(k)=coe2*uu(iisn+1) + (1.-coe2)*uu(iisn)
        uy_e(k)=coe2*vv(iisn+1) + (1.-coe2)*vv(iisn)
      end do

compute th00,tt00,pr00,rh00 and average stability for base state profiles
      th00=th_e(1)
      tt00=tm_e(1)
      tvirt=tm_e(1)*(1.+.6*qv_e(1))
      rh00=press(1)*100./(rg*tvirt)
      pr00=press(1)*100.
      sum=0.
        do k=2,l-1
          sum = sum + (th_e(k+1)-th_e(k-1))/th_e(k)
        enddo
      st=sum/(float(l-2)*2.*dz*gac(k))
c      print*,'th00,tt00,pr00,rh00,st: ',th00,tt00,pr00,rh00,st

compute reference state vertical profiles 
      cap=rg/cp
      capi=1./cap
      cs=gg/(cp*tt00*st)
      do k=1,l
      exs=exp(-st*zz(k))
      th0(k)=th00/exs
      rho0(k)=rh00*exs*(1.-cs*(1.-exs))**(capi-1.)
      enddo

c      print*,'PROFILES'
c      do k=1,l
c       print 200,zz(k)/1.e3,th0(k),rho0(k),th_e(k),
c     .             tm_e(k),qv_e(k)*1.e3,ux_e(k)
c 200    format(1x,'z,th0,rho0,the,tme,qve,ue:',
c     .          2f7.1,f5.2,2f7.1,e10.3,f6.1)
c      enddo

      return
      end

      subroutine zstrtch(zz,nz1,dz)
cc define vertical grid for stretched coordinate, derive Jacobian
      dimension zz(nz1)
      include 'param.grid'
      common /strtch/ height(nz),gac(nz)
      dimension zb(400)
      if(nz.gt.400) stop 'grid formulation'

      top=float(nz-1)*dz

cc exponent of stretching function:
c      ex1=4./3.
c      ex1=4.1/2.9
      ex1=1.
cc
       if(ex1.ne.1.) stop 'fix surface flux for strtching'
cc

      aa=top / top**ex1

cc vertical coordinate:
      do k=1,nz
      zb(k)=aa*zz(k)**ex1
      height(k)=zb(k)
      enddo

cc jacobian:
      do k=1,nz
      if(k.eq.1 ) gac(k)=(zb(2)-zb(1))/dz
      if(k.eq.nz) gac(k)=(zb(nz)-zb(nz-1))/dz
      if(k.ne.1 .and. k.ne.nz) gac(k)=(zb(k+1)-zb(k-1))/(2.*dz)
      print*,k,zb(k),gac(k)
      enddo

cc check consistency (int gac ddzeta = H)
      sum1=.5*(gac(1)*dz + gac(nz)*dz)
      do i=2,nz-1
      sum1=sum1+gac(i)*dz
      enddo
      print*,'int Jacobian before adj: ',sum1

cc adjust numerical jacobian:
      coe=float(nz-1)*dz/sum1
      do i=1,nz
      gac(i)=gac(i)*coe
      enddo

cc check:
      sum1=.5*(gac(1)*dz + gac(nz)*dz)
      do i=2,nz-1
      sum1=sum1+gac(i)*dz
      enddo
      print*,'int Jacobian after adj: ',sum1

      return
      end


      function cvmgm(a,ab,ac)
       if(ac.lt.0) then
        cvmgm=a
       else
        cvmgm=ab
       endif
       return
       end

      subroutine gridint
     1 (ain,nx,nz,xx,zz,nx1,nz1,xx1,zz1,x0,z0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc this subroutine performs interpolation of the data given in the
cc input array ain(nx,nz) on the grid defined by xx(nx) and zz(nz)
cc into the grid given by xx1(nx1) and zz1(nz1). NIETHER OF THE
cc GRIDS HAS TO BE REGULAR. Data is returned in the ain(nx1,nz1)
cc part of the input array.
cc 
cc    levels in the input array are given in zz(nz), 
cc    levels in the output array are given in zz1(nz1)
cc      x-coordinate in the input array are in xx(nx)
cc      x-coordinate in the output array are in xx1(nx1)
cc        x0(nx,nz) and z0(nx,nz) are working arrays
cc
cc NOTE that nx1 (or nz1) must be smaller than nx (or nz) and xx1 (zz1)
cc  must be a subdomain of xx (zz)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dimension ain(nx,nz),zz(nz),xx(nx)
      dimension zz1(nz1),xx1(nx1) 
      dimension x0(nx,nz),z0(nx,nz)

c      check consistency of the input data
      ier=0
cc nx1,nz1 not larger than nx,nz
      if(nx1.gt.nx) ier=1
      if(nz1.gt.nz) ier=2
cc limits (zz1(nz1).le.zz(nz) ?) 
      if(zz1(1).lt.zz(1)) ier=3
      if(zz1(nz1).gt.zz(nz)) ier=4
cc limits (xx1(nx1).le.xx(nx) ?) 
      if(xx1(1).lt.xx(1)) ier=5
      if(xx1(nx1).gt.xx(nx)) ier=6
      if(ier.ne.0) then
      print 999,ier
 999  format(2x,' ** problems with input data. will stop.'/
     1 ' ier = ',i3,'. see code why stoped.')
      stop
      endif
cc
      nxz=nx*nz
      do 99 i=1,nxz
      z0(i,1)=1.
  99  x0(i,1)=1.
cc  map vertical grid positions:
      do 1 k1=1,nz1
      zzh=zz1(k1)
      do 2 k=1,nz
      kk=k
      if(zz(k).ge.zzh) go to 6
  2   continue
  6   kkm=max0(1,kk-1)
      z0(1,k1)=float(kkm)+(zzh-zz(kkm))/(zz(kk)-zz(kkm)+1.e-6)
  1   continue
      do 3 i1=2,nx1
      do 3 k1=1,nz1
  3   z0(i1,k1)=z0(1,k1)
c
cc  map horizontal grid positions:
      do 11 i1=1,nx1
      xxh=xx1(i1)
      do 12 i=1,nx
      ii=i
      if(xx(i).ge.xxh) go to 16
 12   continue
 16   iim=max0(1,ii-1)
      x0(i1,1)=float(iim)+(xxh-xx(iim))/(xx(ii)-xx(iim)+1.e-6)
 11   continue
      do 13 i1=1,nx1
      do 13 k1=2,nz1
 13   x0(i1,k1)=x0(i1,1)
cc
cc  call Piotr's interpolation routine
      call inter2(ain,x0,z0)
cc
      return
      end
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INTER2(XF,XD1,XD2)
C IOR=ORDER OF ACCURACY/2; ONLY EVEN ORDER TRMBACK SCHEMES ARE CONSIDERED

      PARAMETER(IOR=2)
      PARAMETER(LINER=0)
      DIMENSION XF(*),XD1(*),XD2(*)
CC  N1 - HORIZONTAL INDEX, N2 - VERTICAL INDEX
      include 'param.grid'
      PARAMETER (N1=nx, N2=nz,NN=N1*N2)
      DIMENSION Z(NN,-IOR:IOR)
      DATA  EP/ 1.E-10/
      PARAMETER(NONOS=1)
      REAL      MX,MN
      PARAMETER(IBC=0)
      COMMON // IG0(NN),JG0(NN),X(-IOR+1:N1+IOR,-IOR+1:N2+IOR)
C  next is for shavano: 
C      DONOR(Y1,Y2,A)=CVMGM(Y2,Y1,A)*A
C  next is for workstation:
      DONOR(Y1,Y2,A)=AMAX1(0.,A)*Y1 + AMIN1(0.,A)*Y2
      TR2(Y1,Y2,A)=A*.5*(Y1+Y2)-A**2*.5*(Y2-Y1)
      TR4(YM1,Y0,YP1,YP2,A)=A/12.*(7.*(YP1+Y0)-(YP2+YM1))
     1 -A**2/24.*(15.*(YP1-Y0)-(YP2-YM1))-A**3/12.*((YP1+Y0)
     2 -(YP2+YM1))+A**4/24.*(3.*(YP1-Y0)-(YP2-YM1))
      TR6(YM2,YM1,Y0,YP1,YP2,YP3,A)=-A/60.*(-YM2+8.*YM1-37.*Y0
     1                                     -37.*YP1+8.*YP2-YP3)
     2-A**2/360.*(-2.*YM2+25.*YM1-245.*Y0+245.*YP1-25.*YP2+2.*YP3)
     3-A**3/48.*(YM2-7.*YM1+6.*Y0+6.*YP1-7.*YP2+YP3)
     4-A**4/144.*(YM2-11.*YM1+28.*Y0-28.*YP1+11.*YP2-YP3)
     5-A**5/240.*(-YM2+3.*YM1-2.*Y0-2.*YP1+3.*YP2-YP3)
     6-A**6/720.*(-YM2+5.*YM1-10.*Y0+10.*YP1-5.*YP2+YP3)
      PP(XI)=AMAX1(0.,XI)
      PN(XI)=AMIN1(0.,XI)
C
CC CHECK COSISTENCY OF THE DATA:
      IF(NX.NE.N1.OR.NZ.NE.N2) THEN
      PRINT 777
 777  FORMAT(2X,'!!! CALLS TO INTER2 WITH NON-MATCHING DIMENSIONS.'
     1  ,' STOP.')
      STOP
      ENDIF
CC
      DO 1 K=1,NN
      IG0(K)=NINT(XD1(K))
    1 JG0(K)=NINT(XD2(K))

C  GRID EXTENSION FOR BC REMOVAL 
      DO 508 I=1,N1      
      DO 509 J=1,N2
      II=(J-1)*N1+I
  509 X(I,J)=XF(II) 
      DO 5091 IS=1-IOR,0
C     II=(1-1)*N1+I
 5091 X(I,IS)=XF(I)
      DO 5092 IS=1,IOR
      II=(N2-1)*N1+I
 5092 X(I,N2+IS)=XF(II)
  508 CONTINUE
      DO 507 J=-IOR+1,N2+IOR
      DO 5071 IS=-IOR+1,1
 5071 X(IS,J)=X(1,J)*(1-IBC)+IBC*X(N1+IS-1,J)
      DO 5072 IS=0,IOR
 5072 X(N1+IS,J)=X(N1,J)*(1-IBC)+IBC*X(1+IS,J)
  507 CONTINUE
C  END OF GRID EXTENSION
C                     
C
C  HERE STARTS REZIDUAL ADVECTION
C
                     DO 50 J=-IOR,IOR
C
      IF(LINER.EQ.1) THEN
      DO 211 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
  211 Z(II,J)=Y0-(FL1-FL0) 
      GO TO 50
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 311 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  311 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 321 II=1,NN
      U=IG0(II)-XD1(II)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      F0=TR2(YM1, Y0,U)
      F1=TR2(Y0 ,YP1,U)
  321 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 312 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  312 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 322 II=1,NN
      U=IG0(II)-XD1(II)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      F0=TR4(YM2,YM1,Y0 ,YP1,U)
      F1=TR4(YM1,Y0 ,YP1,YP2,U)
  322 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 313 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
      FL0=DONOR(YM1, Y0,U)
      FL1=DONOR(Y0 ,YP1,U)
      W=Y0-(FL1-FL0) 
      MX=AMAX1(YM1,Y0,YP1,W)
      MN=AMIN1(YM1,Y0,YP1,W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
  313 Z(II,J)=W-(F1-F0) 
        ELSE
      DO 323 II=1,NN
      U=IG0(II)-XD1(II)
      YM3=X(IG0(II)-3,JG0(II)+J)
      YM2=X(IG0(II)-2,JG0(II)+J)
      YM1=X(IG0(II)-1,JG0(II)+J)
      Y0 =X(IG0(II)  ,JG0(II)+J)
      YP1=X(IG0(II)+1,JG0(II)+J)
      YP2=X(IG0(II)+2,JG0(II)+J)
      YP3=X(IG0(II)+2,JG0(II)+J)
      F0=TR6(YM3,YM2,YM1,Y0 ,YP1,YP2,U)
      F1=TR6(YM2,YM1,Y0 ,YP1,YP2,YP3,U)
  323 Z(II,J)=Y0-(F1-F0) 
        ENDIF
      ENDIF
C
C
   50 CONTINUE
C  
      IF(LINER.EQ.1) THEN
      DO 212 II=1,NN
      U=JG0(II)-XD2(II)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
  212 XF(II)=Z(II,0)-(FL1-FL0) 
      RETURN
      ENDIF
C
      IF(IOR.EQ.1) THEN
        IF(NONOS.EQ.1) THEN
      DO 411 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  411 CONTINUE
        ELSE
      DO 421 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR2(Z(II,-1),Z(II,0),U)
      F1=TR2(Z(II, 0),Z(II,1),U)
  421 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.2) THEN
        IF(NONOS.EQ.1) THEN
      DO 412 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  412 CONTINUE
        ELSE
      DO 422 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR4(Z(II,-2),Z(II,-1),Z(II,0),Z(II,1),U)
      F1=TR4(Z(II,-1),Z(II, 0),Z(II,1),Z(II,2),U)
  422 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF

      IF(IOR.EQ.3) THEN
        IF(NONOS.EQ.1) THEN
      DO 413 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
      FL0=DONOR(Z(II,-1),Z(II,0),U)
      FL1=DONOR(Z(II, 0),Z(II,1),U)
      W=Z(II,0)-(FL1-FL0) 
      MX=AMAX1(Z(II,-1),Z(II,0),Z(II,1),W)
      MN=AMIN1(Z(II,-1),Z(II,0),Z(II,1),W)
      F0=F0-FL0
      F1=F1-FL1
      OV=(MX-W)/(-PN(F1)+PP(F0)+EP)
      UN=(W-MN)/( PP(F1)-PN(F0)+EP)
      OV=AMIN1(1.,OV)
      UN=AMIN1(1.,UN)
      F0=PP(F0)*OV+PN(F0)*UN
      F1=PP(F1)*UN+PN(F1)*OV
      XF(II)=W-(F1-F0) 
  413 CONTINUE
        ELSE
      DO 423 II=1,NN
      U=JG0(II)-XD2(II)
      F0=TR6(Z(II,-3),Z(II,-2),Z(II,-1),Z(II,0),
     1                     Z(II, 1),Z(II, 2),U)
      F1=TR6(Z(II,-2),Z(II,-1),Z(II, 0),Z(II,1),
     1                     Z(II, 2),Z(II, 3),U)
  423 XF(II)=Z(II,0)-(F1-F0) 
        ENDIF
      ENDIF
      RETURN
      END   

      subroutine sd_adv(ux,uxp,uz,uzp,xp,zp,npx)
      include 'param.grid'
      include 'param.sds'
      dimension ux(nx,nz),uxp(nx,nz),uz(nx,nz),uzp(nx,nz)
      dimension xp(np),zp(np)

      common/grid/ time,dt,dx,dz,dti,dxi,dzi

           do ip=1,npx
cc centered-in-time trajectory
cc predictor:
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1

          uup=.5*(ux(i+1,k)+ux(i,k))
          uum=.5*(ux(i-1,k)+ux(i,k))
          delx=(xp(ip) + dx/2. - (i-1)*dx )/dx
c           print*,'-delx: ',delx
           if(delx.lt.0. .or. delx.gt.1.) stop '111'
          uxc=(1.-delx)*uum + delx*uup
          xxf=xp(ip)+uxc*dt

          wwp=.5*(uz(i,k)+uz(i,k+1))
          wwm=.5*(uz(i,k)+uz(i,k-1))
          delz=(zp(ip) + dz/2. - (k-1)*dz )/dz
           if(delz.lt.0. .or. delz.gt.1.) stop '222'
          uzc=(1.-delz)*wwm + delz*wwp
          zzf=zp(ip)+uzc*dt

cc corrector:
c            do iter=1,0  ! no iteration
            do iter=1,1  ! one iteration
c            do iter=1,2  ! two iterations

          i=(xxf+dx/2.)/dx + 1
          k=(zzf+dz/2.)/dz + 1
          uup=.5*(ux(i+1,k)+ux(i,k))
          uum=.5*(ux(i-1,k)+ux(i,k))
          uup1=.5*(uxp(i+1,k)+uxp(i,k))
          uum1=.5*(uxp(i-1,k)+uxp(i,k))
          uup2=2.*uup-uup1
          uum2=2.*uum-uum1
          wwp=.5*(uz(i,k)+uz(i,k+1))
          wwm=.5*(uz(i,k)+uz(i,k-1))
          wwp1=.5*(uzp(i,k)+uzp(i,k+1))
          wwm1=.5*(uzp(i,k)+uzp(i,k-1))
          wwp2=2.*wwp-wwp1
          wwm2=2.*wwm-wwm1

          delx=(xxf + dx/2. - (i-1)*dx )/dx
           if(delx.lt.0. .or. delx.gt.1.) stop '333'
          uxc2=(1.-delx)*uum2 + delx*uup2
          xxf=xp(ip)+(uxc+uxc2)*dt/2.

          delz=(zzf + dz/2. - (k-1)*dz )/dz
           if(delz.lt.0. .or. delz.gt.1.) stop '444'
          uzc2=(1.-delz)*wwm2 + delz*wwp2
          zzf=zp(ip)+(uzc+uzc2)*dt/2.

                enddo ! iterative corrector

            xp(ip)=xxf
            zp(ip)=zzf

         enddo   !     do ip=1,npx

             return
             end

      subroutine adjust(th,qv,qc,sup,nsd,xp,zp,qcp,plic,radp)
      include 'param.grid'
      include 'param.sds'
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
c SDs:
      dimension xp(np),zp(np)
      dimension qcp(np),plic(np),radp(np)

c SDs parameters assigned to gridboxes:
      dimension th(nx,nz),qv(nx,nz),sup(nx,nz),qc(nx,nz)
      dimension beta(nx,nz),epsilon(nx,nz),epsilon1(nx,nz)

      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlatv,hlats

      a=rg/rv
      c=hlatv/cp
      b=hlatv/(rv*tt0)
      d=hlatv/rv
      e=-cp/rg
      xxlv=hlatv

       do i=1,nxz
       epsilon(i,1) = 0.
       epsilon1(i,1) = 0.
       enddo

       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       do i=1,nx
         thi=1./th(i,k)
         y=b*thetme*tt0*thi
         ees=ee0*exp(b-y)
         qvs=a*ees/(prek-ees)

        abss=sup(i,k)

        tt=th(i,k)/thetme
c        dqsdt = xxlv*qvs/(rv*tt**2)  ! dqs/dT
        dqsdt=xxlv*qvs/(rv*tt**2) * prek/(prek-ees) ! dqs/dT
        ab = 1.+dqsdt*xxlv/cp        ! correction to drop growth rate due to latent heating

c adjust to current supersaturation
         epsilon(i,k) = (qv(i,k)-qvs-abss)/ab

c limit adjustment to available water
         epsilon(i,k)=max(epsilon(i,k),-qc(i,k))

c don't adjust upward if subsaturated
         if (qv(i,k)/qvs.lt.1.) epsilon(i,k)=amin1(0.,epsilon(i,k))

c do not adjust if no cloud water:
         if (epsilon(i,k).gt.0. .and. qc(i,k).lt.1.e-10) epsilon(i,k)=0.

       enddo
       enddo

cc calculate beta, Eq. 12 in Grabowki and Morrison MWR 2008
       do k=1,nz
       do i=1,nx
       beta(i,k)=0.
       enddo
       enddo
             do ip=1,nsd
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
           beta(i,k)=beta(i,k)+plic(ip)*radp(ip) 
             enddo

cc adjustment to cloud water:
        pi=4.*atan(1.)
        rho_w=1.e3
        fact=4./3.*pi*rho_w
             do ip=1,nsd
           qcpbef=qcp(ip)
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
        if(qc(i,k).gt.1.e-8) then
         epsi=epsilon(i,k) * plic(ip)*radp(ip)/beta(i,k)
           epsi=max(epsi,-qcp(ip))
             epsilon1(i,k)=epsilon1(i,k)+epsi 
           qcp(ip) = qcp(ip)+epsi
           radp(ip)=(qcp(ip)/fact/plic(ip)*rho0(k))**.33333
        endif
             enddo

       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       do i=1,nx
c now do the adjustment
         qv(i,k)=qv(i,k)-epsilon1(i,k)
         th(i,k)=th(i,k)+epsilon1(i,k)*thetme*xxlv/cp
       enddo
       enddo

      return
      end

      subroutine sd_phy(theta,qv,sup,fsup,ww,ftt,fqv,
     1            nsd,xp,zp,qcp,plic,radp,iadjust)
      include 'param.grid'
      include 'param.sds'
      common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)
      common /prof_m/ qv_e(nz),tm_e(nz)
      dimension sup(nx,nz),fsup(nx,nz),ftt(nx,nz),fqv(nx,nz),ww(nx,nz)
      dimension theta(nx,nz),qv(nx,nz)
c SDs:
      dimension xp(np),zp(np)
      dimension qcp(np),plic(np),radp(np)

c SDs parameters assigned to gridboxes:
      real npxz(nx,nz),coxz(nx,nz),tausqe(nx,nz)

      common/grid/ time,dt,dx,dz,dti,dxi,dzi
      common /const/ gg,cp,rg,rv
      common/reference/ tt0,ee0
      common/latent/hlat,hlats

c grid in supers to get the same mulitplicity
      dimension sup_sdb(nppg),rad_sdb(nppg)
c parameters for activation:
      common /activ/ sup_sdb,rad_sdb,dcon

c flag for activation:
      dimension iflg(nx,nz)
      common /p_deriv/ dpdz(nz)

       bl(f00,f10,f01,f11,x,y) =
     1   f00*(1.-x)*(1.-y) + f10*x*(1.-y) + f01*(1.-x)*y + f11*x*y

cc activation and condensation growth:
       aa1=0.9152e-10 ! constant in droplet growth equation
       r_sh=1.86e-6   ! including kinetic effects
c       r_sh=0.        ! excluding kinetic effects

        pi=4.*atan(1.)
        rho_w=1.e3
        fact=4./3.*pi*rho_w
      a=rg/rv
      b=hlat/(rv*tt0)
      c=hlat/cp
      d=hlat/rv
      e=-cp/rg
      xxlv=hlat

cc assign number of SDs and local concentartion to gridboxes
       do i=1,nxz
       npxz(i,1)=0.
       coxz(i,1)=0.
       tausqe(i,1)=0.
       enddo
           do ip=1,nsd
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
          npxz(i,k)=npxz(i,k)+1.
          con=plic(ip)
          coxz(i,k)=coxz(i,k)+con
          tausqe(i,k)=tausqe(i,k) + radp(ip)*con
           enddo
cccc constant coefficient:
c       do i=1,nxz
c       tausqe(i,1)=max(tausqe(i,1),.1)
c       tausqe(i,1)=1./(2.8e-4*tausqe(i,1))
c       enddo
ccc variable coefficient:
       do k=1,nz
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       do i=1,nx
       tt=theta(i,k)/thetme
       dv = 8.794e-5*tt**1.81/prek
       coef=4.*pi*dv
       tausqe(i,k)=max(tausqe(i,k),.1)
       tausqe(i,k)=1./(coef*tausqe(i,k))
       enddo
       enddo

           if(iadjust.eq.1) then
cccc quasi-analytic approach for supersaturation:
         do i=1,nx 
         do k=1,nz 
cc place when mixing will be added...
          fqv(i,k)=0.
          ftt(i,k)=0.
        thetme=th_e(k)/tm_e(k)
        prek=1.e5*thetme**e
        thi=1./theta(i,k)
        tt=theta(i,k)/thetme
        y=b*thetme*tt0*thi
        ees=ee0*exp(b-y)
        qvs=a*ees/(prek-ees)
        dqsdt=xxlv*qvs/(rv*tt**2) * prek/(prek-ees) ! dqs/dT
c my:
c        ah=qvs/(prek-ees) * 
c     1    (dpdz(k) + xxlv*prek/(rv*tt**2)*gg/cp) * ww(i,k)
cc my with forces:
        ah=dqsdt*(gg/cp*ww(i,k) - ftt(i,k)/thetme)
     1    + qvs/(prek-ees)*dpdz(k)*ww(i,k) + fqv(i,k)

        xx=1./tausqe(i,k)
        sstemp=ah/xx+(sup(i,k)-ah/xx)*exp(-xx*dt)
        fsup(i,k)=(sstemp-sup(i,k))/dt

ccc recalculate sup as qv/qvs-1 for microphysics:
         sup(i,k)=sup(i,k)/qvs

         enddo 
         enddo 

           endif  ! iadjust

cc activation:
         do i=1,nx 
         do k=1,nz 

cc flag for activation:
         iflg(i,k) = 0

         if(sup(i,k).le.0.) go to 388

           supers=sup(i,k)

cc local concentration and concentration based on supers:
           conc=coxz(i,k)
           call activat(supers,concin,1)

           if(concin .le. conc) go to 388  ! no need to activate...
cc new or additional activation:
             ipin=(concin+dcon/10.)/dcon
             ipin=min(ipin,nppg)
             ipl =conc /dcon
            if(ipl .eq. ipin) go to 388  ! still no need to activate...

       do ip=ipl+1,ipin  ! creating new SDs

         nsd=nsd+1
         ipp=nsd
cccc random position at activation
           rand1=rand()-.5
           rand2=rand()-.5
cc
        xp(ipp)=(i-1)*dx + rand1*dx
        zp(ipp)=(k-1)*dz + rand2*dz

        plic(ipp)=dcon*rho0(k)
        radp(ipp)=rad_sdb(ip)
        qcp(ipp)=fact * radp(ipp)**3.*dcon
        dqc=qcp(ipp)/dt
        fqv(i,k)=fqv(i,k)-dqc
        ftt(i,k)=ftt(i,k)+dqc*hlat/cp*th_e(k)/tm_e(k)
       enddo

cc flag for activation:
         iflg(i,k) = 1

388          continue ! no need to activate...
             if(nsd.gt.np) then
               print*,'after activation',i,k,np,nsd
               STOP 'stop nsd.gt.np'
              endif

         enddo  ! k loop
         enddo  ! i loop

cc reposition droplets after activation:
       call sd_rep_act(xp,zp,nsd,iflg)


cc condensation/evaporation:
cc interp=0/1 - no/yes interpolation to SD position
           interp=0
c           interp=1
           do ip=1,nsd
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
           if(interp.eq.0) then
cc no interpolation
           supers=sup(i,k)
           else
cc interpolation
       ii=xp(ip)/dx + 1
       kk=zp(ip)/dz + 1

        u00=sup(ii,kk)
        u10=sup(ii+1,kk)
        u01=sup(ii,kk+1)
        u11=sup(ii+1,kk+1)

        xsd=(xp(ip) - (ii-1)*dx )/dx
        zsd=(zp(ip) - (kk-1)*dz )/dz

c        if(xsd.lt.0. .or. xsd.gt.1.) stop '111'
c        if(zsd.lt.0. .or. zsd.gt.1.) stop '222'

       supers=bl(u00,u10,u01,u11,xsd,zsd)
           endif

            rad_o=radp(ip)
cc constant:
c       aaa=aa1
cc variable
       thetme=th_e(k)/tm_e(k)
       prek=1.e5*thetme**e
       tt=theta(i,k)/thetme
       dv = 8.794e-5*tt**1.81/prek
       thi=1./theta(i,k)
       y=b*thetme*tt0*thi
       ees=ee0*exp(b-y)
       qvs=a*ees/(prek-ees)
         dqsdt=xxlv*qvs/(rv*tt**2) * prek/(prek-ees) ! dqs/dT
c         dqsdt = xxlv*qvs/(rv*tt**2)  ! dqs/dT
         ab = 1.+dqsdt*xxlv/cp
       aaa=qvs*dv/rho_w/ab
        
      radn2=amax1(0., (rad_o+r_sh)**2 + 2.*aaa*supers*dt)
      rad_n=sqrt(radn2) - r_sh
         radp(ip)=max(0.,rad_n)
         qcp(ip)=fact * radp(ip)**3. *plic(ip)/rho0(k)
          dqc=fact * (radp(ip)**3.-rad_o**3.)
     1             * plic(ip)/rho0(k) / dt
           fqv(i,k)=fqv(i,k)-dqc
           ftt(i,k)=ftt(i,k)+dqc*hlat/cp*th_e(k)/tm_e(k)

cc logic for total evaporation:
         if(rad_n.lt.1.e-8) then
           radp(ip)=0.
           qcp(ip)=0.
           plic(ip)=0.
         endif
           enddo   ! ip loop

cc eliminate empty SDs
       call remove(xp,zp,qcp,plic,radp,nsd)

           return
           end

       subroutine remove(xp,zp,qcp,plic,radp,nsd)
      include 'param.grid'
      include 'param.sds'
cc this routine removes evaporated super-droplets
      dimension xp(np),zp(np)
      dimension qcp(np),plic(np),radp(np)
      dimension xpn(np),zpn(np)
      dimension qcpn(np),plicn(np),radpn(np)

cc zero new locations:
           do ip=1,np
           xpn(ip)=0.
           zpn(ip)=0.
           qcpn(ip)=0.
           plicn(ip)=0.
           radpn(ip)=0.
           enddo

             nsdn=0
             do ip=1,nsd
             if(radp(ip).gt.1.e-7) then
             nsdn=nsdn+1
               xpn(nsdn)=xp(ip)
               zpn(nsdn)=zp(ip)
               qcpn(nsdn)=qcp(ip)
               plicn(nsdn)=plic(ip)
               radpn(nsdn)=radp(ip)
             endif
             enddo

              if(nsdn.gt.np) stop 'np too small'

             nsd=nsdn
             do ip=1,nsd
               xp(ip)=xpn(ip)
               zp(ip)=zpn(ip)
               qcp(ip)=qcpn(ip)
               plic(ip)=plicn(ip)
               radp(ip)=radpn(ip)
             enddo
             do ip=nsd+1,np
               xp(ip)=0.
               zp(ip)=0.
               qcp(ip)=0.
               plic(ip)=0.
               radp(ip)=0.
             enddo

      return
      end

      subroutine sd_rep_act(xp,zp,npx,iflg)
      include 'param.grid'
      include 'param.sds'
      dimension iflg(nx,nz)
      dimension xp(np),zp(np)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
           do ip=1,npx
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
           if(iflg(i,k).eq.1) then
           rand1=rand()-.5
           rand2=rand()-.5
          xp(ip)=(i-1)*dx +  rand1*dx
          zp(ip)=(k-1)*dz +  rand2*dz
           endif
           enddo
          return
          end

      subroutine sd_rep(xp,zp,npx)
      include 'param.grid'
      include 'param.sds'
      dimension xp(np),zp(np)
      common/grid/ time,dt,dx,dz,dti,dxi,dzi
           do ip=1,npx
          i=(xp(ip)+dx/2.)/dx + 1
          k=(zp(ip)+dz/2.)/dz + 1
           rand1=rand()-.5
           rand2=rand()-.5
          xp(ip)=(i-1)*dx +  rand1*dx
          zp(ip)=(k-1)*dz +  rand2*dz
           enddo
          return
          end

       subroutine activat(sup_inp,activ,ifl)
c ifl=0   read data
c ifl=1   get N for input S
c ifl=2   get S for input N
cc A. Ziemniak data, see in ~/Documents/AZIMNIAK
       parameter(nv=271) 
       dimension sup(nv),conc(nv)
       dimension suppl(nv),concpl(nv)
       common /activationAZ/ sup,conc
cc needed for conversion to mixing ratio on input
       include 'param.grid'
       common /prof_d/ rho0(nz),th0(nz),th_e(nz),ux_e(nz),uy_e(nz)

       if(ifl.eq.0) then
cc read N-S data:
       open(14,file='N-S_VOCALS.txt',status='old')

       read(14,*)
       do iv=1,nv
       read (14,*) sup(iv),conc(iv)
       suppl(iv)=sup(iv)*100.
       conc(iv)=conc(iv)/rho0(1) ! converted to m.r. using surface rho
       concpl(iv)=conc(iv)*1.e-6
       enddo

       close(14)

cc plot N-S:
      call setusv('LW',2000)
      call set(.15,.95,.1,.7,0.01,10.,0.,120.,3)
      call labmod('(f5.2)','(f4.0)',5,4,2,2,20,20,0)
      call periml(1,9,3,4)
      call dashdc('$$$$$$$$$$$$$$$$$$$$',10,12)
      call curved(suppl,concpl,nv)
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      CALL plchhq(.55,0.03, 'supersaturation (%)', 0.016,0.,0)
      CALL plchhq(.03,0.40, 'activated CCN (cm**-3)', 0.016,90.,0)
      call frame

       return
       endif

       if(ifl.eq.1) then
cc return number of activated droplets given the input supersaturation
       sup1=sup_inp
cc deal with small and large values:
       if(sup1.lt.0.) then
       activ=0.
       return
       endif

       if(sup1.le.sup(1)) then
c       activ=conc(1)*sup1/sup(1)
       activ=0.
       return
       endif

       if(sup1.gt.sup(nv)) then
       print*,' something wrong with activation. stop'
       print*,'sup, sup(last): ',sup1,sup(nv)
       stop
       endif

cc interpolation needed:
cc find class:
        do ic=2,nv
        if(sup(ic).gt.sup1) then
        iv=ic
        go to 101
        endif
        enddo
         stop 'cannot find class'
101     continue
cc interpolate:
       coe=(sup1-sup(iv-1))/(sup(iv)-sup(iv-1))
       activ=coe*conc(iv) + (1.-coe)*conc(iv-1)
       return
       endif

       if(ifl.eq.2) then
cc return supersaturation given the input number of activated droplets 

       activ1=activ
cc deal with small and large values:
       if(activ1.le.conc(1)) then
c       sup_inp=activ1*sup(1)/conc(1)
       sup_inp=sup(1)
       return
       endif

       if(activ1.gt.conc(nv)) then
       sup_inp=sup(nv)
       endif

cc interpolation needed:
cc find class:
        do ic=2,nv
        if(conc(ic).gt.activ1) then
        iv=ic
        go to 103
        endif
        enddo
         stop 'cannot find class'
103     continue
cc interpolate:
       coe=(activ1-conc(iv-1))/(conc(iv)-conc(iv-1))
       sup_inp=coe*sup(iv) + (1.-coe)*sup(iv-1)
       return
       endif

       end


