      program linden
      implicit double precision (a-h,o-z)

      parameter(maxl=500)
      parameter(taunit=1.4387752224d+00,bh2=59.3220,
     +     pi=3.14159265358979323846d+00,eps=1.d-16)
      character argum*30
c POR density matrix in terms cos(deltaphi)
c     write(6,'(''read in T, Nslice,bconst,npt,iodevn'')')
c     read(5,*)temprt,nslice,bconst,npt,iodevn
      call getarg(1,argum)
      read(argum,*)temprt
      call getarg(2,argum)
      read(argum,*)nslice
      call getarg(3,argum)
      read(argum,*)bconst
      call getarg(4,argum)
      read(argum,*)npt
      call getarg(5,argum)
      read(argum,*)iodevn
c     ... calculate tau
      tau=taunit/(temprt*nslice)
      write(6,'(''tau='',f10.5)')tau

c     ... judge the maximam l quantum number: l means m hereforth
      do l=0,maxl
         if(exp(-tau*bconst*l*(l)).lt.eps)then
            write(6,*)l,exp(-tau*bconst*l*(l)),-tau*bconst*l*(l)
            lmax=l
            goto 20
         endif
      enddo
      write(6,'(''!!! Warning: maxl is reached'')')
      lmax=maxl
 20   continue
      write(6,'(''lmax='',i4)')lmax

      cstep=2.0/dfloat(npt-1)

c     ... get the norm at cost=1.0
      cost=-1.0

      open(7,file='qpr.out',status='unknown')

      do ic=1,npt
         cost=(ic-1)*cstep-1.d0
         call exarho(cost,lmax,maxl,rho,erot1,tau,bconst,iodevn,
     +        nslice,erotsq)
         write(7,'(1p,7(1x,E15.8))')cost,rho,erot1,erotsq
      enddo
      close(7,status='keep')

c     ... calculate the rotational energy at beta
      cost=1.0
      beta=tau*nslice
      write(6,*)'beta=',beta
      call exarho(cost,lmax,maxl,rho,erot,beta,bconst,iodevn,nslice,
     +     erotsq)
      erot=erot*nslice/rho
      erotsq=erotsq*nslice*nslice/rho
      Cv=(erotsq-erot*erot)/temprt**2.d0
      write(6,'(A,E15.8)')'Erot at Beta:',erot
      write(6,'(A,E15.8)')'Cv at Beta:',Cv


      end

c-----------------------------------------------------------------------
      subroutine exarho(cost,lmax,maxl,rho,erot,tau,bconst,iodevn,
     +     nslice,erotsq)
c     ... calculate the exact linear rotor propagator
      implicit double precision(a-h,o-z)
      parameter(pi=3.14159265358979323846d+00,boltz=0.69503476d0)

      rho=0.d0
      erot=0.d0
      erotsq=0.d0

      dt=dacos(cost)
c     ... sum over l: l means m!
      do l=1,lmax
            tmp=exp(-tau*bconst*dfloat(l*(l)))*dcos(dfloat(l)*dt)
            rho=rho+tmp
            erot=erot+tmp*l*(l)*bconst
            erotsq=erotsq+tmp*l*(l)*bconst*l*(l)*bconst
      enddo
c     erot=erot/rho
c     ... convert erot to the unit of K and divide by the number of slice
      erot=2.d0*erot/(nslice*boltz)
c     ... convert erotsq to the unit of K^2 and divide by the square of number of slice
      erotsq=2.d0*erotsq/(nslice*boltz)**2.d0
      rho=2.*rho+1.d0

      rho=rho/(2.0*pi)
      erot=erot/(2.0*pi)
      erotsq=erotsq/(2.0*pi)
      
      return
      end
