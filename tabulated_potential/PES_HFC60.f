       program main
       integer i1,i2,j1,j2,i0
       double precision EulangJ(2),EulangL(2),EHFC60,R

       do i0=1,1
       
       R=0.11d0+(i0-1)*0.01d0

       do j1=1,1
       do i1=1,181
       do j2=1,1
       do i2=1,181

!minimum:    0.1100  148.3000  180.0000   31.7000    0.0000       -6.92393606
!            0.1100   90.0000  180.0000   90.0000    0.0000       -6.92393608
 
       EulangL(1)=0.d0+(i1-1)*1.d0
       EulangL(2)=180.d0+(j1-1)*2.d0

       EulangJ(1)=0.d0+(i2-1)*1.d0
       EulangJ(2)=0.d0+(j2-1)*2.d0
      
       call enHFC60(R,EulangL,EulangJ,EHFC60)
       write(*,70) R,EulangL,EulangJ,EHFC60

       enddo
       enddo
       enddo
       enddo
       enddo

   70  format(5F10.4,1F18.8)

       end program 



      subroutine enHFC60(R,EulangL,EulangJ,EHFC60)
c
c     input: coordinates of the center of mass of HF molecule are defined by 
c            R, AngleL(1)=th, AngleL(2)=ph; rotation of HF in the centre of mass
c            coordinate system is defined by AngleJ(1)=thp and AngleJ(2)=php;
c            R in Angstrom and Angles in radians
c     output: energy in kcal/mol  
c     This potential is based on the DF-LMP2/cc-pVTZ Counterpoise calculations 

      integer i,j,k,n,p1,p2,p,mp
      parameter(n=12,kcal2k = 503.219565d0)
      double precision EulangL(2),EulangJ(2),EHFC60,x(n),t,R,R2,R3,R4
      double precision th,thp,ph,php 
      dimension p1(n),p2(n),p(n),mp(n)
 
      th=EulangL(1)
      ph=EulangL(2)
      thp=EulangJ(1)
      php=EulangJ(2)

         	p1(1)=0    
         	p1(2)=0    
         	p1(3)=0    
         	p1(4)=0    
         	p1(5)=1    
         	p1(6)=1    
         	p1(7)=1    
         	p1(8)=1    
         	p1(9)=2    
         	p1(10)=2   
         	p1(11)=2   
         	p1(12)=2   
                           
         	p2(1)=0    
         	p2(2)=6    
         	p2(3)=6    
         	p2(4)=6    
         	p2(5)=1    
         	p2(6)=5    
         	p2(7)=5    
         	p2(8)=5    
         	p2(9)=2    
         	p2(10)=4   
         	p2(11)=4   
         	p2(12)=4   
                           
         	p(1)=0     
         	p(2)=6     
         	p(3)=6     
         	p(4)=6     
         	p(5)=0     
         	p(6)=6     
         	p(7)=6     
         	p(8)=6     
         	p(9)=0     
         	p(10)=6    
         	p(11)=6    
         	p(12)=6   
                           
         	mp(1)=0    
         	mp(2)=0    
         	mp(3)=5    
         	mp(4)=-5   
         	mp(5)=0    
         	mp(6)=0    
         	mp(7)=5    
         	mp(8)=-5   
         	mp(9)=0    
         	mp(10)=0   
         	mp(11)=5   
         	mp(12)=-5  
                           

              
       R2=R*R
       R4=R2*R2
       R3=R*R2

       x(1)=  5.825731d0*R4 + 3.505364d0*R2 - 6.871599d0   
       x(2)=  2.927476d-3*R4+1.101096d-4*R2 + 9.487742d-4
       x(3)=-sqrt(7.d0/11.d0)*x(2)
       x(4)=sqrt(7.d0/11.d0)*x(2)


       x(5)=  3.229325d-1*R3-1.619958d-1*R
       x(6)=  -1.172073d-2*R3-2.621524d-3*R             
       x(7)=-sqrt(7.d0/11.d0)*x(6)
       x(8)=sqrt(7.d0/11.d0)*x(6)


       x(9)= 8.895730d-2*R4+1.076448d-2*R2            
       x(10)= -6.372600d-3*R4+5.922033d-3*R2            
       x(11)=-sqrt(7.d0/11.d0)*x(10)
       x(12)=sqrt(7.d0/11.d0)*x(10)
      

       EHFC60=0.d0
       do i=1,n
       EHFC60=EHFC60+x(i)*t(p1(i),p2(i),p(i),mp(i),th,ph,thp,php)
       enddo
       
      end

      function t(p1,p2,p,mp,th,ph,thp,php)
      implicit none
      integer p1,q1,p2,p,mp,r1,r2,r
      real*8 t,th,ph,thp,php,PLMx,
     $     xp1,xp2,xp,xr1,xr2,xr,thrj,pi
      parameter (Pi=3.14159265358979324d0)

      t=0
      do r1=-p1,p1
         do r2=-p2,p2

            r=mp
            if (abs(r).gt.p) cycle


            xp1=p1; xp2=p2; xp=p
            xr1=r1; xr2=r2; xr=r


            if (thrj(xp1,xp2,xp,xr1,xr2,-xr).eq.0.0) cycle


            t = t + thrj(xp1,xp2,xp,xr1,xr2,-xr)
     $           * PLMx(p2,r2,cos(thp*Pi/180.0d0))
     $           * PLMx(p1,r1,cos(th*Pi/180.0d0))
     $           *(cos((r2*php+r1*ph)*Pi/180.0d0))
     $           *sqrt((2.d0*xp1+1.d0)*(2.d0*xp2+1.d0))
     $           *(-1.d0)**(xp1-xp2+xr)*(sqrt(2.d0*xp+1.d0))*2.

         enddo
      enddo
      end



      function PLMx(L,M,C)
      implicit none
      integer L, M
      real*8 PLMx, PLM, C, farity
      PLMx=PLM(L,M,C)
      if (M.lt.0) then
         PLMx=PLMx*farity(M)
      endif
      end

      FUNCTION PLM(LIN,MIN,COSTH)
C
C     COMPUTES NORMALIZED ASSOC. LEGENDRE POLYNOMIALS BY RECURSION.
C     THE VALUES RETURNED ARE NORMALIZED FOR INTEGRATION OVER X
C     (I.E. INTEGRATION OVER COS THETA BUT NOT PHI).
C     NOTE THAT THE NORMALIZATION GIVES
C           PLM(L,0,1)=SQRT(L+0.5)
C           PLM(L,0,X)=SQRT(L+0.5) P(L,X)
C     FOR M.NE.0, THE VALUE RETURNED DIFFERS FROM THE USUAL
C           DEFINITION OF THE ASSOCIATED LEGENDRE POLYNOMIAL
C           (E.G. EDMONDS PAGES 23-24)
C           BY A FACTOR OF (-1)**M*SQRT(L+0.5)*SQRT((L-M)!/(L+M)!)
C     THUS THE SPHERICAL HARMONICS ARE
C          CLM = PLM * EXP(I*M*PHI) / SQRT(L+0.5)
C          YLM = PLM * EXP(I*M*PHI) / SQRT(2*PI)
C     THIS ROUTINE ALWAYS RETURNS THE VALUE FOR ABS(MIN); NOTE THAT
C          FOR MIN.LT.0 THIS VALUE SHOULD BE MULTIPLIED BY FARITY(MIN)
C
C     FUNCTION PM1(LIN,MIN,COSTH)
C     This routine appears to be much more stable for large l, m than
C       the routine from Nerf/ modified according to R.T Pack
C     It was obtained:
C     From: Marie-Lise Dubernet <mld@ipp-garching.mpg.de>
C     Date: Mon, 19 Jun 1995 12:48:11 +0200 (MET DST)
C     Some mods 27-28 June 95 by SG for speed and to accord w/ MOLSCAT 
C     Bugs fixed 21 Sept 95 (SG)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     CHECK FOR ABS(COSTH).LE.1.D0 ... 
      IF (ABS(COSTH).GT.1.D0) THEN
        WRITE(6,*) ' *** ILLEGAL ARGUMENT TO PLM. X =',COSTH
        STOP
      ENDIF
C     SAVE ARGUMENTS IN LOCAL VARIABLES
      L=LIN
      M=ABS(MIN)
      X=COSTH
C
C  IF M>L PLM=0 !
      IF(M.GT.L) THEN
        PLM=0.D0
        RETURN
      ENDIF
      LMAX=L
C
      IF (M.GT.0) GO TO 5
C  HERE FOR REGULAR LEGENDRE POLYNOMIALS
      PLM=1.D0
      PM2=0.D0
      XL=0.D0
      DO 2 L=1,LMAX
      XL=XL+1.D0
      PP=((2.D0*XL-1.D0)*X*PLM-(XL-1.D0)*PM2)/XL
      PM2=PLM
    2 PLM=PP
      GO TO 9000
C
C  HERE FOR ALEXANDER-LEGENDRE POLYNOMIALS
C
    5 IMAX=2*M
      RAT=1.D0
      AI=0.D0
      DO 6 I=2,IMAX,2
      AI=AI+2.D0
    6 RAT=RAT*((AI-1.D0)/AI)
C     Y=SIN(THETA)
      Y=SQRT(1.D0-X*X)
      PLM=SQRT(RAT)*(Y**M)
      PM2=0.D0
      LOW=M+1
      XL=LOW-1
      DO 10 L=LOW,LMAX
      XL=XL+1.D0
      AL=DBLE((L+M)*(L-M))
      AL=1.D0/AL
      AL2=(DBLE((L+M-1)*(L-M-1)))*AL
      AL=SQRT(AL)
      AL2=SQRT(AL2)
      PP=(2.D0*XL-1.D0)*X*PLM*AL-PM2*AL2
      PM2=PLM
   10 PLM=PP
      PLM=PLM*FARITY(MIN)
C
C     CONVERT TO MOLSCAT'S IDIOSYNCRATIC NORMALIZATION
9000  PLM=PLM*SQRT(XL+0.5D0)
      RETURN
      END

      FUNCTION THRJ(F1,F2,F3,G1,G2,G3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SMALL CHANGES 31 JUL 95 (SG)
      SAVE MUNG,X,Y
      PARAMETER (MXIX=302)
      DIMENSION X(MXIX),Y(MXIX)
      DATA MUNG/0/
      IF (MUNG.EQ.21) GO TO 69
      MUNG = 21
      X(1) = 0.D0
      DO 100 I = 1, MXIX-1
      A = I
      X(I+1) = LOG(A) +X(I)
      Y(I+1) = LOG(A)
  100 CONTINUE
   69 IF(F1-ABS(G1)) 1,13,13
   13 IF(F2-ABS(G2))1,14,14
   14 IF(F3-ABS(G3))1,15,15
   15 SUM=F1+F2+F3
      NSUM=SUM+.001D0
      IF(SUM-NSUM)2,2,1
    1 THRJ=0.D0
      RETURN
    2 IF(ABS(G1+G2+G3)-1.D-08)3,3,1
    3 IF(F1+F2-F3)1,4,4
    4 IF(F1+F3-F2)1,5,5
    5 IF(F2+F3-F1)1,6,6
    6 J1=2.D0*F3+2.001D0
      J2=F1+F2-F3+1.001D0
      J3=F1-F2+F3+1.001D0
      J4=-F1+F2+F3+1.001D0
      J5=F1+F2+F3+2.001D0
      J6=F1+G1+1.001D0
      J7=F1-G1+1.001D0
      J8=F2+G2+1.001D0
      J9=F2-G2+1.001D0
      J10=F3+G3+1.001D0
      J11=F3-G3+1.001D0
      IF(J5.GT.MXIX) THEN
        WRITE(6,601) J5,MXIX
  601   FORMAT(' *** DIMENSION ERROR IN THRJ - INDEX.GT.MXIX',2I5)
        STOP
      ENDIF
      R=0.5D0*(Y(J1)+X(J2)+X(J3)+X(J4)-X(J5)
     1+X(J6)+X(J7)+X(J8)+X(J9)+X(J10)+X(J11))
      SUM=0.D0
      F=-1
      KZ=-1
    7 KZ=KZ+1
      F=-F
      J1=KZ+1
      J2=F1+F2-F3-KZ+1.001D0
      IF(J2)20,20,8
    8 J3=F1-G1-KZ+1.001D0
      IF(J3)20,20,9
    9 J4=F2+G2-KZ+1.001D0
      IF(J4)20,20,10
   10 J5=F3-F2+G1+KZ+1.001D0
      IF(J5)7,7,11
   11 J6=F3-F1-G2+KZ+1.001D0
      IF(J6)7,7,12
   12 JMAX=MAX(J1,J2,J3,J4,J5,J6)
      IF(JMAX.GT.MXIX) THEN
        WRITE(6,601) JMAX,MXIX
        STOP
      ENDIF
      S=-(X(J1)+X(J2)+X(J3)+X(J4)+X(J5)+X(J6))
      SUM=SUM+F*EXP(R+S)
      GO TO 7
   20 INT=ABS(F1-F2-G3)+0.0001D0
      VAL=((-1.D0)**INT)*SUM/SQRT(2.D0*F3+1.D0)
      IF(ABS(VAL).LE.1.D-6) VAL=0.D0
      THRJ=VAL
      RETURN
      END

      FUNCTION FARITY(I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FARITY=1.D0
      IF((I/2)*2-I.NE.0) FARITY=-1.D0
      RETURN
      END


