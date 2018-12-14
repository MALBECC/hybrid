c This subroutine computes forces exerted by electrons and nucleus
c on the point charges
c----------------------------------------------------------------------------
      subroutine mmforce(na,nac,natot,nesp,ntpl,nspin,ntm,isa,list,r,
     .                   f,pc,ucell,rcorte,dvol,Rho)

      use neutralatom
      implicit real*8 (a-h,o-z)
      integer na,nac,natot,nesp,isa(na),list(nac)
      integer ntpl,nspin,ntm(3)
      real*8 rcorte,vna,grvna(3),fel(3),rr(3)
      real*8 r(3,natot),ucell(3,3),f(3,natot),dvol,pc(nac)
      real Rho(ntpl)
        f=0.0d0
        rcorte2=rcorte**2
c --------------------------------------------------------------------------
c loop over MM atoms
        do js=na+1,natot

c Forces only for ListQMMM solvent atoms 
      if(list(js-na).eq.0) then
        xq=r(1,js)
        yq=r(2,js)
        zq=r(3,js)
     
c loop over the mesh points 
      do iz=0,ntm(3)-1
      do iy=0,ntm(2)-1
      do ix=0,ntm(1)-1

      xm = ucell(1,1) * ix/ntm(1) + ucell(1,2) * iy/ntm(2)
     .   + ucell(1,3) * iz/ntm(3)
      ym = ucell(2,1) * ix/ntm(1) + ucell(2,2) * iy/ntm(2)
     .   + ucell(2,3) * iz/ntm(3)
      zm = ucell(3,1) * ix/ntm(1) + ucell(3,2) * iy/ntm(2)
     .   + ucell(3,3) * iz/ntm(3)

      dx = xq - xm
      dy = yq - ym
      dz = zq - zm

      d2 = dx**2 + dy**2 + dz**2       

      imesh = 1 + ix + ntm(1)*iy + ntm(1)*ntm(2)*iz

C Forces exerted on point charges by electrons 
      if (d2.gt.rcorte2) then
      d=sqrt(d2)
      De = -2.*dvol*Rho(imesh)*pc(js-na)/d**3
      else
      De = 0.0 
      endif
 
       f(1,js) = f(1,js) + De*dx
       f(2,js) = f(2,js) + De*dy
       f(3,js) = f(3,js) + De*dz
       
       enddo !!grid
       enddo
       enddo

C loop over QM atoms
        j2=js
      do  j1=1,na
        is = isa(j1)
        ra = rcut(is)
       tx=r(1,j1)-r(1,j2)
       ty=r(2,j1)-r(2,j2)
       tz=r(3,j1)-r(3,j2)
        rr(1)=tx
        rr(2)=ty
        rr(3)=tz
       dd2=tx**2 + ty**2 + tz**2
       dd=sqrt(dd2)

C Forces exerted on point charges by nucleus
        fel=0.0
        if (dd.lt.ra) then
          call vna_sub(is,rr,vna,grvna)
          fel(1:3)=grvna(1:3)*pc(j2-na)
        endif

        f(1,j2)=f(1,j2) - fel(1) 
        f(2,j2)=f(2,j2) - fel(2)
        f(3,j2)=f(3,j2) - fel(3)

        enddo !! at st

        endif !! list

        enddo  !! at sv

       return
       end

