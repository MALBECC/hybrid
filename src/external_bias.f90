    subroutine external_bias(bias_type,natot,rclas,fdummy,E)
    use scarlett, only : Ang
    implicit none
    integer, intent(in) :: bias_type, natot
    double precision, dimension(3,natot), intent(in) :: rclas
    double precision, dimension(3,natot), intent(inout) :: fdummy
    double precision, intent(inout) :: E
    integer :: i
    if (bias_type .eq. 1) then !Müller-Brown Potential
      do i=1, natot
        call MB(rclas(1,i)/Ang,rclas(2,i)/Ang,fdummy(1,i), fdummy(2,i),E)
      end do
    elseif (bias_type .eq. 2) then !Estrin Potential
      do i=1, natot
        call EDA(rclas(1,i)/Ang,fdummy(1,i),E)
      end do
    end if
    
    end subroutine external_bias


    subroutine EDA(x,fx,E)
    use scarlett, only : eV,kcal
    implicit none
    double precision, dimension(4) :: a,b,c,af
    double precision, intent(in) :: x
    double precision, intent(out) :: E, fx
    integer :: i
    a=(/15.d0,-30.d0, -50.d0, 25.d0/)
    af=(/-6.d0,6.d0, -200.d0, 200.d0/)
    b=(/0.4d0, 0.2d0, -4.d0 , -8.d0/)
    c=(/-2.4d0, -1.2d0, 8.d0, 16.d0/)
    do i=1,4
      E=E+eval_exp(x,a(i),b(i),c(i))
      fx=fx+eval_exp(x,af(i),b(i),c(i))
    end do
    E=E*5.12d-2
    fx=fx*5.12d-2

    contains
    double precision function eval_exp(x,a,b,c)
    implicit none
    double precision, intent(in) :: x,a,b,c
    eval_exp=0.d0
    eval_exp=a*2.71d0**(b*x+c)
    return
    end function eval_exp
    end subroutine EDA


    subroutine MB(x,y,fx,fy,E)
    use scarlett, only : eV,kcal
    implicit none
    double precision, intent(in) :: x,y
    double precision, intent(inout) :: fx,fy,E
    double precision, dimension(4) :: Ap,a,b,c,x0,y0
    integer :: i
    Ap= (/-200d0,-100d0,-170d0,15d0 /)
    Ap=Ap*eV/kcal
    a=(/-1d0,-1d0,-6.5d0,0.7d0/)
    b=(/0.d0,0.d0,11.d0,0.6d0/)
    c=(/-10d0,-10d0,-6.5d0,0.7d0/)
    x0=(/1.d0,0.d0,-0.5d0,-1.d0/)
    y0=(/0.d0,0.5d0,1.5d0,1.d0/)
    do i=1,4
      E=E+energy(x,y,Ap(i),a(i),b(i),c(i),x0(i),y0(i))
      fx=fx+f_x(x,y,Ap(i),a(i),b(i),c(i),x0(i),y0(i))
      fy=fy+f_y(x,y,Ap(i),a(i),b(i),c(i),x0(i),y0(i))
    end do
    
    contains
    
    double precision function energy(x,y,Ap,a,b,c,x0,y0)
    implicit none
    double precision, intent(in) :: x,y,Ap,a,b,c,x0,y0
    energy=0.d0
    energy=Ap*exp(a*(x-x0)**2+b*(x-x0)*(y-y0)+c*(y-y0)**2)
    return
    end function energy
    
    double precision function f_x(x,y,Ap,a,b,c,x0,y0)
    implicit none
    double precision, intent(in) :: x,y,Ap,a,b,c,x0,y0
    f_x=0.d0
    f_x=Ap*exp(a*(x-x0)**2+b*(x-x0)*(y-y0)+c*(y-y0)**2)
    f_x=-f_x*(2.d0*a*(x-x0)+b*(y-y0))
    return
    end function f_x
    
    double precision function f_y(x,y,Ap,a,b,c,x0,y0)
    implicit none
    double precision, intent(in) :: x,y,Ap,a,b,c,x0,y0
    f_y=0.d0
    f_y=Ap*exp(a*(x-x0)**2+b*(x-x0)*(y-y0)+c*(y-y0)**2)
    f_y=-f_y*(2.d0*c*(y-y0)+b*(x-x0))
    return
    end function f_y
    
    end subroutine MB
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
