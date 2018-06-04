c     ==========================================================
      subroutine cluster(com_1,com_2,Eulang_1,Eulang_2,E_12)
c     ==========================================================
c     Initially:
c     __________________________________________________
c     this subroutine calculates pair potential between
c     rigid water molecules encapsulated in C60 given the
c     coordinates of centre of mass of water molecules and
c     the Euler angles describing the rotation of each H2O
c     input: com_1, com_2, Eulang_1, Eulang_2 
c     output: E_12
c     phi=Eulang(1)
c     theta=Eulang(2)
c     chi=Eulang(3)
c     __________________________________________________
c     rotmat, computed within the code with Eulang
c     ROwf, RH1wf, RH2wf, RCwf put as data   
      implicit double precision(a-h,o-z)
      dimension delta(3,3),  
     +          com_1(3), com_2(3),d(3),d1(3),d2(3),
     +          Eulang_1(3),Eulang_2(3),R12(3),
     +          rotmat_1(3,3),rotmat_2(3,3),
     +          ROwf(3),RH1wf(3),RH2wf(3),
     +          RO_1_sf(3),RO_2_sf(3),
     +          RH1_1_sf(3),RH1_2_sf(3),
     +          RH2_1_sf(3),RH2_2_sf(3),
     +          alfa_c60(3,3),alfa_c60_1(3,3),alfa_c60_2(3,3),
     +          alfa_h2oc60(3,3),alfa_h2oc60_1(3,3),alfa_h2oc60_2(3,3),
     +          quad(3,3),q1(3,3),q2(3,3),
     +          om(3,3,3),om1(3,3,3),om2(3,3,3),
     +          F(3,3,3,3),F1(3,3,3,3),F2(3,3,3,3)
      parameter (kcal2k=503.218978939d0,polautoA3=0.14819d0,
     +          quadtoDebA=(2.54177d0*0.529177d0),a0=0.529177d0,
     +          am0=2.54175d0)
      autocm=219474.63d0 !from hartree to cm-1
      autoK=315775.13d0 !hartree to Kelvin
      autokcal=627.503d0!hartree to kcal/mol
      autokJ=2625.5d0!hartree to kJ/mol
c permanent dipole moment of h2oc60  
      data d/-0.0011d0,0.0d0,0.2031d0/ !in ea0 pbe0/def2-tzvpp 
c polarizability of c60  pbe0/def2-tzvpp
      data axx_1/546.07d0/,ayy_1/546.07d0/,azz_1/546.07d0/!in a.u.
c polarizability of h2oc60  pbe0/def2-tzvpp
      data axx/547.15d0/,ayy/547.11d0/,azz/547.18d0/!in a.u.
c permanent quadrupole moment of h2oc60  pbe0/def2-tzvpp
      data qxx/0.602/,qyy/-0.603/,qzz/0.001/,qxz/-0.011/!in a.u.
c permanent octupole moment of h2oc60 pbe0/def2-tzvpp
      data oxxx/-0.382d0/,oxyy/0.367d0/,oxxz/1.648d0/,oyyy/-0.002d0/,
     1     oyyz/-0.822d0/,oxzz/0.016d0/,ozzz/-0.826d0/
c permanent hexadecapole moment of h2oc60 pbe0/def2-tzvpp
      data Fxxxx/0.244d0/,Fyyyy/2.384d0/,Fzzzz/-3.612d0/,
     1     Fxxyy/-3.120d0/,Fxxzz/2.876d0/,Fyyzz/0.736d0/



      c6_c60=84589.40d0
      c6_h2oc60=88467.16d0
      alpha_h2oc60=(axx+ayy+azz)/3d0
      alpha_c60=(axx_1+ayy_1+azz_1)/3d0 


      om(1:3,1:3,1:3)=0.d0
      om(1,1,1)=oxxx
      om(1,2,2)=oxyy
      om(1,1,3)=oxxz
      om(2,2,3)=oyyz
      om(1,3,3)=oxzz
      om(3,3,3)=ozzz
      om(2,2,2)=oyyy
      om(1,3,1)=om(1,1,3)
      om(3,1,1)=om(1,1,3)
      om(2,1,2)=om(1,2,2)
      om(2,2,1)=om(1,2,2)
      om(2,3,2)=om(2,2,3)
      om(3,2,2)=om(2,2,3)
      om(3,1,3)=om(1,3,3)
      om(3,3,1)=om(1,3,3)


 
      F(1:3,1:3,1:3,1:3)=0.d0
      F(1,1,1,1)=Fxxxx
      F(2,2,2,2)=Fyyyy
      F(3,3,3,3)=Fzzzz
      F(1,1,2,2)=Fxxyy
      F(1,1,3,3)=Fxxzz
      F(2,2,3,3)=Fyyzz
      F(2,3,2,3)=F(2,2,3,3)
      F(2,3,3,2)=F(2,2,3,3)
      F(3,2,2,3)=F(2,2,3,3)
      F(3,3,2,2)=F(2,2,3,3)
      F(3,2,3,2)=F(2,2,3,3)
      F(1,3,1,3)=F(1,1,3,3)
      F(1,3,3,1)=F(1,1,3,3)
      F(3,1,1,3)=F(1,1,3,3)
      F(3,3,1,1)=F(1,1,3,3)
      F(3,1,3,1)=F(1,1,3,3)
      F(1,2,1,2)=F(1,1,2,2)
      F(1,2,2,1)=F(1,1,2,2)
      F(2,1,1,2)=F(1,1,2,2)
      F(2,2,1,1)=F(1,1,2,2)
      F(2,1,2,1)=F(1,1,2,2)



c  coordinates of water in body fixed frame in AA
      data ROwf/0.d0,0.d0,-0.06563807d0/,  
     +     RH1wf/0.7575d0,0.d0,0.52086193d0/,
     +     RH2wf/-0.7575d0,0.d0,0.52086193d0/

c delta function     
      do i=1,3
        do j=1,3
          if (i.eq.j) then
            delta(i,j)=1.0d0
          else
            delta(i,j)=0.0d0
          endif
        enddo
      enddo

      E_12=0.d0
c polarizability of water molecule in AA**3      
      do i=1,3
       do j=1,3
         alfa_h2oc60(i,j)=0.d0
         alfa_c60(i,j)=0.d0
         quad(i,j)=0.d0
       enddo 
      enddo
      alfa_h2oc60(1,1)=axx
      alfa_h2oc60(2,2)=ayy
      alfa_h2oc60(3,3)=azz
      alfa_c60(1,1)=axx_1
      alfa_c60(2,2)=ayy_1
      alfa_c60(3,3)=azz_1
      quad(1,1)=qxx
      quad(2,2)=qyy
      quad(3,3)=qzz
      quad(1,3)=qxz
      quad(3,1)=quad(1,3)
c      write(*,*)quad(1,1),quad(2,2),quad(3,3)

c     prepare rotational matrix for water1 
c     and coordinates of atoms in space sixed frame
      call matpre(Eulang_1, rotmat_1)
      do i=1,3
         RO_1_sf(i)=0.d0
         RH1_1_sf(i)=0.d0
         RH2_1_sf(i)=0.d0
      enddo
      call rottrn(rotmat_1, ROwf, RO_1_sf, com_1)
      call rottrn(rotmat_1, RH1wf, RH1_1_sf, com_1)
      call rottrn(rotmat_1, RH2wf, RH2_1_sf, com_1)
  

c     prepare rotational matrix for water2
c     and coordinates of atoms in space sixed frame
      call matpre(Eulang_2, rotmat_2)
      do i=1,3
         RO_2_sf(i)=0.d0
         RH1_2_sf(i)=0.d0
         RH2_2_sf(i)=0.d0
      enddo
      call rottrn(rotmat_2, ROwf, RO_2_sf, com_2)
      call rottrn(rotmat_2, RH1wf, RH1_2_sf, com_2)
      call rottrn(rotmat_2, RH2wf, RH2_2_sf, com_2)

c dipole moment of each water@c60 molecule in space fixed frame
      d1(1:3)=0.d0
      d2(1:3)=0.d0
      do i=1,3
       do j=1,3
         d1(i)=d1(i)+rotmat_1(i,j)*d(j)
         d2(i)=d2(i)+rotmat_2(i,j)*d(j)
       enddo
      enddo

c polarizabilty tesors of each water molecule in space fixed frame
      alfa_h2oc60_1(1:3,1:3)=0.d0
      alfa_h2oc60_2(1:3,1:3)=0.d0
      alfa_c60_1(1:3,1:3)=0.d0
      alfa_c60_2(1:3,1:3)=0.d0
      q1(1:3,1:3)=0.d0
      q2(1:3,1:3)=0.d0
      do i1=1,3
       do i2=1,3
        do j1=1,3
         do j2=1,3
         alfa_h2oc60_1(i2,j2)=alfa_h2oc60_1(i2,j2)+
     1                rotmat_1(i2,i1)*rotmat_1(j2,j1)*alfa_h2oc60(i1,j1)
         alfa_h2oc60_2(i2,j2)=alfa_h2oc60_2(i2,j2)+
     1                rotmat_2(i2,i1)*rotmat_2(j2,j1)*alfa_h2oc60(i1,j1)
         alfa_c60_1(i2,j2)=alfa_c60_1(i2,j2)+
     1                rotmat_1(i2,i1)*rotmat_1(j2,j1)*alfa_c60(i1,j1)
         alfa_c60_2(i2,j2)=alfa_c60_2(i2,j2)+
     1                rotmat_2(i2,i1)*rotmat_2(j2,j1)*alfa_c60(i1,j1)
         q1(i2,j2)=q1(i2,j2)+rotmat_1(i2,i1)*rotmat_1(j2,j1)*quad(i1,j1)
         q2(i2,j2)=q2(i2,j2)+rotmat_2(i2,i1)*rotmat_2(j2,j1)*quad(i1,j1)
         enddo
        enddo
       enddo
      enddo

c polarizabilty tesors of each water molecule in space fixed frame
      om1(1:3,1:3,1:3)=0.d0
      om2(1:3,1:3,1:3)=0.d0
      do i1=1,3
       do i2=1,3
        do j1=1,3
         do j2=1,3
          do k1=1,3
           do k2=1,3  
         om1(i2,j2,k2)=om1(i2,j2,k2)+rotmat_1(i2,i1)*rotmat_1(j2,j1)
     1                 *rotmat_1(k2,k1)*om(i1,j1,k1)
         om2(i2,j2,k2)=om2(i2,j2,k2)+rotmat_2(i2,i1)*rotmat_2(j2,j1)
     1                 *rotmat_2(k2,k1)*om(i1,j1,k1)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      F1(1:3,1:3,1:3,1:3)=0.d0
      F2(1:3,1:3,1:3,1:3)=0.d0
      do i1=1,3
       do i2=1,3
        do j1=1,3
         do j2=1,3
          do k1=1,3
           do k2=1,3
            do l1=1,3
             do l2=1,3
         F1(i2,j2,k2,l2)=F1(i2,j2,k2,l2)+rotmat_1(i2,i1)*rotmat_1(j2,j1)
     1                 *rotmat_1(k2,k1)*rotmat_1(l2,l1)*F(i1,j1,k1,l1)
         F2(i2,j2,k2,l2)=F2(i2,j2,k2,l2)+rotmat_2(i2,i1)*rotmat_2(j2,j1)
     1                 *rotmat_2(k2,k1)*rotmat_2(l2,l1)*F(i1,j1,k1,l1)
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo



      R12(1:3)=0.0
      do i=1,3
         R12(i)=com_1(i)-com_2(i)
      enddo
      R=dsqrt(R12(1)*R12(1)+R12(2)*R12(2)+R12(3)*R12(3))

c ... pair interaction potential
c dipole-dipole interaction
      v3_elec=0.d0
      do i=1,3
       do j=1,3
        v3_elec=v3_elec - d1(i)*d2(j)*R**(-5.d0)*(3.d0*R12(i)*R12(j)-
     1 R**2.d0*delta(i,j))   
       enddo
      enddo
c dipole-quadrupole interaction
      v4_elec=0.d0
      do i=1,3
       do j=1,3
        do k=1,3
       Tijk=-3.d0*R**(-7.d0)*(5.d0*R12(i)*R12(j)*R12(k)-
     1 R**2.d0*(R12(i)*delta(j,k)+R12(j)*delta(i,k)+R12(k)*delta(i,j)))
        v4_elec=v4_elec-(-d1(i)*q2(j,k)+d2(i)*q1(j,k))*Tijk/3d0
        enddo
       enddo
      enddo

c dipole-octupole interaction
      v5_elec=0.d0
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
        Tijkl=3d0*(35d0*R12(i)*R12(j)*R12(k)*R12(l)-5d0*R**2d0*(R12(i)*
     1  R12(j)*delta(k, l)+R12(i)*R12(k)*delta(j, l)+R12(i)*R12(l)*
     1  delta(j, k)+R12(j)*R12(k)*delta(i, l)+R12(j)*R12(l)*delta(i, k)+
     1  R12(k)*R12(l)*delta(i, j))+R**4d0*(delta(i, j)*delta(k, l)+
     1  delta(i, k)*delta(j, l)+delta(i, l)*delta(j, k)))/R**9d0

        v5_elec=v5_elec-(d1(i)*om2(j,k,l)+d2(i)*om1(j,k,l))*Tijkl/15d0
     1                 +(q1(i,j)*q2(k,l))*Tijkl/9d0
         enddo
        enddo
       enddo
      enddo

c dipole-hexadecapole interaction
      v6_elec=0.d0
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
           do m=1,3
        Tijklm=-15d0*(63d0*R12(i)*R12(j)*R12(k)*R12(l)*R12(m)
     1  -7d0*R**2d0*(R12(j)*R12(k)*R12(l)*delta(i, m)+R12(i)*R12(k)*
     1  R12(l)*delta(j, m)+R12(i)*R12(j)*R12(l)*delta(k, m)+R12(i)*
     1  R12(j)*R12(k)*delta(l, m)+R12(m)*R12(j)*R12(k)*delta(i, l)+
     1  R12(m)*R12(i)*R12(k)*delta(j, l)+R12(m)*R12(i)*R12(j)*
     1  delta(k, l)+R12(i)*R12(l)*R12(m)*delta(j, k)+R12(j)*R12(l)*
     1  R12(m)*delta(i, k)+R12(k)*R12(l)*R12(m)*delta(i, j))+R**4d0*
     1  (R12(k)*delta(i, l)*delta(j, m)+R12(j)*delta(i, l)*delta(k, m)+
     1  R12(k)*delta(i, m)*delta(j, l)+R12(i)*delta(k, m)*delta(j, l)+
     1  R12(j)*delta(i, m)*delta(k, l)+R12(i)*delta(j, m)*delta(k, l)+
     1  R12(i)*delta(l, m)*delta(j, k)+R12(l)*delta(i, m)*delta(j, k)+
     1  R12(j)*delta(l, m)*delta(i, k)+R12(l)*delta(j, m)*delta(i, k)+
     1  R12(k)*delta(l, m)*delta(i, j)+R12(l)*delta(k, m)*delta(i, j)+
     1  R12(m)*delta(i, l)*delta(j, k)+R12(m)*delta(j, l)*delta(i, k)+
     1  R12(m)*delta(k, l)*delta(i, j)))/R**11d0

        v6_elec=v6_elec-(d1(i)*F2(j,k,l,m)+d2(i)*F1(j,k,l,m))*
     1          Tijklm/105d0   
     1            -(q1(i,j)*om2(k,l,m)-q2(i,j)*om1(k,l,m))*Tijklm/45d0
          enddo
         enddo
        enddo
       enddo
      enddo
      
c induction(v5 ~R^-6)  and  dispersion (v6 ~R^-6) interactions
      v6_ind=0.d0
      v6_h2oc60=0.d0
      v6_c60=0d0
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
       Tik=R**(-5.d0)*(3.d0*R12(i)*R12(k)-R**2.d0*delta(i,k))
       Tjl=R**(-5.d0)*(3.d0*R12(j)*R12(l)-R**2.d0*delta(j,l))
       v6_ind=v6_ind-0.5d0*(alfa_h2oc60_1(i,j)*d2(k)*d2(l)
     1              +alfa_h2oc60_2(i,j)*d1(k)*d1(l))*Tik*Tjl
       v6_h2oc60=v6_h2oc60-C6_h2oc60/(6.d0*alpha_h2oc60**2d0)*
     1           alfa_h2oc60_1(i,j)*alfa_h2oc60_2(k,l)*Tik*Tjl
       v6_c60=v6_c60-C6_c60/(6.d0*alpha_c60**2d0)*
     1           alfa_c60_1(i,j)*alfa_c60_2(k,l)*Tik*Tjl
         enddo
        enddo
       enddo
      enddo

      v3_elec=v3_elec*a0**3.d0  
      v4_elec=v4_elec*a0**4.d0
      v5_elec=v5_elec*a0**5.d0
      v6_elec=v6_elec*a0**6.d0
      v6_ind=v6_ind*a0**6.d0
      v6_disp=(v6_h2oc60-v6_c60)*a0**6.d0 !in a.u.



c      write(*,*)"dd(R-3)=   ",v3_elec*autoK,
c     1          "  dq(R-4)=   ",v4_elec*autoK
c      write(*,*)"do+qq(R-5)=",v5_elec*autoK,
c     1          "  dF+qo(R-6)=",v6_elec*autoK
c      write(*,*)"add(R-6)=  ",v6_ind*autoK   ,
c     1          "  aa(R-6)=   ",v6_disp*autoK
      
c      open(22,file="ENERGYY",access="append")
      edisp=v6_disp
      eelec=v3_elec+v4_elec+v5_elec+v6_elec
      eind=v6_ind
c      write(*,*) "Eelec=",eelec*autoK,"K   Eind=",eind*autoK,
c     1           "K   Edisp=",edisp*autoK,"K"
c      write(22,*) "Eelec=",eelec*autoK," Eind=",eind*autoK,
c     1           " Edisp=",edisp*autoK," Etot=",(edisp+eelec+eind)*autoK

c      close(22)

      E_12_au=v3_elec+v4_elec+v5_elec+v6_disp+v6_ind+v6_elec  !in Hartree
      E_12=E_12_au*autoK !in Kelvin
      E_12_kcal=E_12_au*autokcal !in Kcal/mol  
      E_12_kJ=E_12_au*autokJ !in kJ/mol
      E_12_cm=E_12_au*autocm !in cm-1      
 
c      write(*,*) E_12
c      write(44,72) RO_1_sf,RH1_1_sf,RH2_1_sf      
c      write(44,73) RO_2_sf,RH1_2_sf,RH2_2_sf
c      write(44,71) com_1(1:3),com_2(1:3),R12(1:3),R
c      write(44,74) E_12,E_12_cm
c      write(44,*)  "            "
   71 format("R1=(",3F8.4,")   R2=(",3F8.4,")   R12=(",3F8.4,")  R=",
     1      1F8.4, "   Angstrom")
   72 format("1st mol  O=(",3F12.8,")  H=(",3F12.8,")  H=(",3F12.8,")")
   73 format("2nd mol  O=(",3F12.8,")  H=(",3F12.8,")  H=(",3F12.8,")")
   74 format("EKel = ",1F10.4,"  Ecm-1 = ",1F10.4)
      return
      end
