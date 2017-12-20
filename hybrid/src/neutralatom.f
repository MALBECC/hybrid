        module neutralatom 

        integer, parameter       :: nemax = 100  ! max. number of chem. species
        integer, parameter       :: nmax = 10000 ! max. number of table lines

        double precision, save   :: delta(nemax) ! Dist between points in table
        double precision, save   :: rcut(nemax)  ! Max radius of Vna
        double precision, save   :: vnatable(nmax,nemax) ! Table with Vna 
        integer,save             :: nvna(nemax)  ! Number of entries in table

        private nemax,delta,vnatable,nvna,nmax
        public rcut, vna_sub 

        contains
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
