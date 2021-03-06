      program linden
      implicit double precision (a-h,o-z)

      parameter(maxl=500)
      parameter(taunit=1.4387752224d+00,bh2=59.3220,
     +          pi=3.14159265358979323846d+00,eps=1.d-16)
      dimension pl(0:maxl)
      character argum*30

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
c ... calculate tau
c.... unit of tau is cm and bconst is cm-1 
      tau=taunit/(temprt*nslice)  
      write(6,'(''tau='',f10.5)')tau

c======================================================================c
c ... judge the maximam l quantum number
      do l=0,maxl
        if(exp(-tau*bconst*l*(l+1)).lt.eps)then
          write(6,*)l,exp(-tau*bconst*l*(l+1)),-tau*bconst*l*(l+1)
          lmax=l
          goto 20
        endif
      enddo
      write(6,'(''!!! Warning: maxl is reached'')')
      lmax=maxl
   20 continue
      write(6,'(''lmax='',i4)')lmax

c======================================================================c
      open(7,file='linden.out',status='unknown')

      cstep=2.0/dfloat(npt-1)
      do ic=1,npt
        cost=(ic-1)*cstep-1.d0
        call exarho(cost,lmax,maxl,pl,rho,erot1,tau,bconst,iodevn,
     +       nslice,erotsq)
        erotf = erot1/rho
        write(7,'(1p,7(1x,E15.8))')cost,rho,erotf,erotsq
      enddo
      close(7,status='keep')

c ... calculate the rotational energy at beta
      cost=1.0
      beta=tau*nslice
      write(6,*)'beta=',beta
      call exarho(cost,lmax,maxl,pl,rho,erot,beta,bconst,iodevn,nslice,
     +            erotsq)
      erot=erot*nslice/rho
      erotsq=erotsq*nslice*nslice/rho
      Cv=(erotsq-erot*erot)/temprt**2.d0
      write(6,'(A,E15.8)')'Erot at Beta:',erot
      write(6,'(A,E15.8)')'Cv at Beta:',Cv
      end

c======================================================================c
      subroutine exarho(cost,lmax,maxl,pl,rho,erot,tau,bconst,iodevn,
     +                  nslice,erotsq)
      implicit double precision(a-h,o-z)
      parameter(pi=3.14159265358979323846d+00,boltz=0.69503476d0)

      dimension pl(0:maxl)

      call lgnd(lmax,cost,pl)

      rho=0.d0
      erot=0.d0
      erotsq=0.d0

      do l=0,lmax
        if((mod(l,2).eq.iodevn).or.iodevn.eq.-1) then
          tmp=dfloat(2*l+1)*pl(l)
          tmp=tmp*exp(-tau*bconst*dfloat(l*(l+1)))
          rho=rho+tmp
          erot=erot+tmp*l*(l+1)*bconst
          erotsq=erotsq+tmp*l*(l+1)*bconst*l*(l+1)*bconst
        endif
      enddo
c     erot=erot/rho
c ... convert erot to the unit of K and divide by the number of slice
      erot=erot/(nslice*boltz)
c ... convert erotsq to the unit of K^2 and divide by the square of number of slice
      erotsq=erotsq/(nslice*boltz)**2.d0
      rho=rho/(4.0*Pi)
      erot=erot/(4.0*pi)
      erotsq=erotsq/(4.0*pi)
      return
      end
c-----------------------------------------------------------------------
      subroutine ratrho(cost,tau,bconst,rho,erot,iodevn,nslice)
      implicit double precision(a-h,o-z)

      parameter(boltz=0.69503476d0)

      if(iodevn.eq.-1)rho=exp((cost-1.0)/(2.0*bconst*tau))
      if(iodevn.eq.0)
     +rho=exp((cost-1.0)/(2.0*bconst*tau))+
     +exp((-cost-1.0)/(2.0*bconst*tau))
      if(iodevn.eq.1)
     +rho=exp((cost-1.0)/(2.0*bconst*tau))-
     +exp((-cost-1.0)/(2.0*bconst*tau))

c     if(iodevn.eq.-1)rho=exp(cost/(2.0*bconst*tau))
c     if(iodevn.eq.0)
c    +rho=exp(cost/(2.0*bconst*tau))+
c    +exp(-cost/(2.0*bconst*tau))
c     if(iodevn.eq.1)
c    +rho=exp(cost/(2.0*bconst*tau))-
c    +exp(-cost/(2.0*bconst*tau))

c ... calculate the rotational energy estimator
      if(iodevn.eq.-1) then
        erot=1.0d0/tau-0.5d0*(1.d0-cost)/(bconst*tau*tau)
      else
        fact1=cost-1.d0
        fact2=-cost-1.d0
        erot=exp(fact1/(2.d0*bconst*tau))*fact1+((-1)**(iodevn))*
     +       exp(fact2/(2.d0*bconst*tau))*fact2
        erot=erot/(
     +       exp(fact1/(2.d0*bconst*tau))+((-1)**(iodevn))*
     +       exp(fact2/(2.d0*bconst*tau))
     +       )
        erot=1.0d0/tau+erot/(2.d0*bconst*tau*tau)
      endif

c ... convert erot to the unit of K and divide by the number of slice
      erot=erot/(nslice*boltz)

      return
      end

ccccccccccccccccccccccccc     Program 5.4     cccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Please Note:                                                         c
c                                                                      c
c (1) This computer program is part of the book, "An Introduction to   c
c     Computational Physics," written by Tao Pang and published and    c
c     copyrighted by Cambridge University Press in 1997.               c
c                                                                      c
c (2) No warranties, express or implied, are made for this program.    c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE LGND(LMAX,X,P)
      implicit double precision(a-h,o-z)
C
C Subroutine to generate Legendre polynomials P_L(X)
C for L = 0,1,...,LMAX with given X.
C
      DIMENSION P(0:LMAX)
      P(0) = 1.
      P(1) = X
      DO 100 L = 1, LMAX-1
        P(L+1) = ((2.0*L+1)*X*P(L)-L*P(L-1))/(L+1)
  100 CONTINUE
      RETURN
      END
