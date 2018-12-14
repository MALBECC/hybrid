        module neutralatom 

        integer, parameter       :: nemax = 100  ! max. number of chem. species
        integer, parameter       :: nmax = 10000 ! max. number of table lines

        double precision, save   :: delta(nemax) ! Dist between points in table
        double precision, save   :: rcut(nemax)  ! Max radius of Vna
        double precision, save   :: vnatable(nmax,nemax) ! Table with Vna 
        integer,save             :: nvna(nemax)  ! Number of entries in table

        private nemax,delta,vnatable,nvna,nmax
        public rcut, vna_sub, read_na

        contains

        subroutine read_na(ne,atsim)

        use ionew
        implicit none
        integer ne
        character atsim(ne)*2, filevna*24, paste*24
        external paste

C internal variables
        integer unit, is, i
        real*8 x(nmax)

C Reads the neutral atom potentials from files 
        if (ne .gt. nemax) stop 'Increase max number of species!!'
        call io_assign(unit)

        do is=1,ne
          call io_assign(unit)
          filevna = paste('VNA.', atsim(is))
          open(unit,file=filevna,status='old',form='formatted')
          read(unit,*)
          do i=1,nmax
            nvna(is)=i
            if (i .gt. nmax) stop 'Increase size of VNA matrix!!'
            read(unit,*,err=100,end=100) x(i),vnatable(i,is)
          enddo
100       continue
          delta(is) = x(2)-x(1)
          rcut(is)=x(nvna(is)-1)
          call io_close(unit)
        enddo

        return
        end subroutine read_na

c..............................................................................

         subroutine vna_sub(is,rr,vna,grvna)

         implicit none
         integer is
         double precision rr(3),vna,dvna,dist,grvna(3),yp1,ypn,y2(nmax)

C  initialize splines tables
         yp1 = huge(1.d0)
         ypn = huge(1.d0)
         call spline(delta(is),vnatable(1,is),nvna(is),yp1,ypn,y2)
C  calculate distance
         dist = dsqrt(rr(1)**2 + rr(2)**2 + rr(3)**2)

C  interpolate
         call splint(delta(is),vnatable(1,is),y2,nvna(is),dist,vna,dvna)

C  calculate gradients
         grvna(1:3) = dvna * rr(1:3)/dist

         return
         end subroutine vna_sub

         end module
