c This subroutine calculates the potential due to solvent
c at each point of the mesh and keeps it in Vqm        
c----------------------------------------------------------------------------
      subroutine pcpot(na,nac,natot,ntpl,ntm,list,r,ucell,pc,rcorte,Vqm)

      implicit real*8 (a-h,o-z)
      integer na,nac,natot,list(nac)
      integer ntpl,ntm(3)
      real*8 r(3,natot),ucell(3,3),pc(nac),rcorte
      real Vqm(ntpl)
      Vqm=0.0
      rcorte2=rcorte**2
c----------------------------------------------------------------------------
c loop over MM atoms
      do js=na+1,natot

c Vqm only for ListQMMM solvent atoms
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

c calculation of the external potential due to point charges
       if (d2.gt.rcorte2) then
       d=sqrt(d2)
       Vqm(imesh)=Vqm(imesh)+pc(js-na)/d 
       else
       Vqm(imesh)=Vqm(imesh)+pc(js-na)/rcorte      
       endif

       enddo !!grid
       enddo
       enddo

        endif !!list

        enddo  !!sv atoms

c change units
       Vqm=-2.*Vqm

       return
       end

