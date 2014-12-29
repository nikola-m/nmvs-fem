  module shape_functions
  use prec
  implicit none

  contains

  function shape_function(xi,eta,n)
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: xi,eta

  real(dp),dimension(n) :: shape_function
  real(dp),dimension(n) :: shapefun

  select case(n)
     case(4)  ! Quad4 element
         shapefun(1)=0.25*(1-xi)*(1-eta)
         shapefun(2)=0.25*(1+xi)*(1-eta)
         shapefun(3)=0.25*(1+xi)*(1+eta)
         shapefun(4)=0.25*(1-xi)*(1+eta)
     case(8)   ! Quad8 Serendipity element
         shapefun(1)=-0.25*(1-xi)*(1-eta)*(1+xi+eta)
         shapefun(2)=  0.5*(1-xi)*(1+xi)*(1-eta)
         shapefun(3)=-0.25*(1+xi)*(1-eta)*(1-xi+eta)
         shapefun(4)=  0.5*(1+xi)*(1+eta)*(1-eta)
         shapefun(5)=-0.25*(1+xi)*(1+eta)*(1-xi-eta)
         shapefun(6)=  0.5*(1-xi)*(1+xi)*(1+eta)
         shapefun(7)=-0.25*(1-xi)*(1+eta)*(1+xi-eta)
         shapefun(8)=  0.5*(1-xi)*(1+eta)*(1-eta)
     case default
        write(*,*) "shape_function ERROR: Non existing shape function type!"
  end select
  shape_function=shapefun
  end function shape_function

  function shape_function_der(xi,eta,n,d)
  implicit none
  integer, intent(in) :: d,n
  real(dp), intent(in) :: xi,eta

  real(dp),dimension(d,n) :: shape_function_der
  real(dp),dimension(d,n) :: dhdxi
  select case(d)
    case(2)
      select case(n)
        case(4) ! Quad4 element
            dhdxi(1,1)=-0.25*(1-eta)
            dhdxi(1,2)=0.25*(1-eta)
            dhdxi(1,3)=0.25*(1+eta)
            dhdxi(1,4)=-0.25*(1+eta)

            dhdxi(2,1)=-0.25*(1-xi)
            dhdxi(2,2)=-0.25*(1+xi)
            dhdxi(2,3)=0.25*(1+xi)
            dhdxi(2,4)=0.25*(1-xi)
        case(8) ! Quad8 Serendipity element
            dhdxi(1,1)=-0.25*(-1+eta)*(2*xi+eta)
            dhdxi(1,2)=   xi*(-1+eta)
            dhdxi(1,3)= 0.25*(-1+eta)*(eta-2*xi)
            dhdxi(1,4)=  -0.5*(1+eta)*(-1+eta)
            dhdxi(1,5)=  0.25*(1+eta)*(2*xi+eta)
            dhdxi(1,6)=   -xi*(1+eta)
            dhdxi(1,7)= -0.25*(1+eta)*(eta-2*xi)
            dhdxi(1,8)=   0.5*(1+eta)*(-1+eta)

            dhdxi(2,1)=-0.25*(-1+xi)*(xi+2*eta)
            dhdxi(2,2)=   0.5*(1+xi)*(-1+xi)
            dhdxi(2,3)=  0.25*(1+xi)*(2*eta-xi)
            dhdxi(2,4)=  -eta*(1+xi)
            dhdxi(2,5)=  0.25*(1+xi)*(xi+2*eta)
            dhdxi(2,6)=  -0.5*(1+xi)*(-1+xi)
            dhdxi(2,7)=-0.25*(-1+xi)*(2*eta-xi)
            dhdxi(2,8)=  eta*(-1+xi)
        case default
            write(*,*) "shape_function_der ERROR: Non existing shape function type!"
      end select   
    case default
      write(*,*) "shape_function_der ERROR: Elements of that dimension not implemented yet!"
    end select
  shape_function_der=dhdxi
  end function shape_function_der


  end module shape_functions
