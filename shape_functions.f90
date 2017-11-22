  module shape_functions
  use prec
  implicit none


  interface shape_function
    module procedure shape_function_2d
    module procedure shape_function_3d
  end interface

  interface shape_function_der
    module procedure shape_function_der_2d
    module procedure shape_function_der_3d
  end interface

  private

  public :: shape_function, shape_function_der

  contains

  function shape_function_2d(xi,eta,n)
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: xi,eta

  real(dp),dimension(n) :: shape_function_2d
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
  shape_function_2d=shapefun
  end function shape_function_2d

  function shape_function_3d(xi,eta,zeta,n)
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: xi,eta,zeta

  real(dp),dimension(n) :: shape_function_3d
  real(dp),dimension(n) :: shapefun
  real(dp) :: xi4

  select case(n)
     case(10)  ! Tet10 element
        xi4 = 1.0_dp-xi-eta-zeta
        shapefun(1) = (2.*xi-1.0_dp)*xi
        shapefun(2) = (2.*eta-1.0_dp)*eta
        shapefun(3) = (2.*zeta-1.0_dp)*zeta
        shapefun(4) = (2.*xi4-1.0_dp)*xi4
        shapefun(5) = 4.*xi*eta
        shapefun(6) = 4.*eta*zeta
        shapefun(7) = 4.*zeta*xi
        shapefun(8) = 4.*xi*xi4
        shapefun(9) = 4.*eta*xi4
        shapefun(10) = 4.*zeta*xi4
     case(20)   ! Brick20 element
        shapefun(1) = (1.0_dp-xi)*(1.0_dp-eta)*(1.0_dp-zeta)*(-xi-eta-zeta-2.)*0.125_dp
        shapefun(2) = (1.0_dp+xi)*(1.0_dp-eta)*(1.0_dp-zeta)*(xi-eta-zeta-2.)*0.125_dp
        shapefun(3) = (1.0_dp+xi)*(1.0_dp+eta)*(1.0_dp-zeta)*(xi+eta-zeta-2.)*0.125_dp
        shapefun(4) = (1.0_dp-xi)*(1.0_dp+eta)*(1.0_dp-zeta)*(-xi+eta-zeta-2.)*0.125_dp
        shapefun(5) = (1.0_dp-xi)*(1.0_dp-eta)*(1.0_dp+zeta)*(-xi-eta+zeta-2.)*0.125_dp
        shapefun(6) = (1.0_dp+xi)*(1.0_dp-eta)*(1.0_dp+zeta)*(xi-eta+zeta-2.)*0.125_dp
        shapefun(7) = (1.0_dp+xi)*(1.0_dp+eta)*(1.0_dp+zeta)*(xi+eta+zeta-2.)*0.125_dp
        shapefun(8) = (1.0_dp-xi)*(1.0_dp+eta)*(1.0_dp+zeta)*(-xi+eta+zeta-2.)*0.125_dp
        shapefun(9)  = (1.0_dp-xi**2)*(1.0_dp-eta)*(1.0_dp-zeta)*0.25_dp
        shapefun(10) = (1.0_dp+xi)*(1.0_dp-eta**2)*(1.0_dp-zeta)*0.25_dp
        shapefun(11) = (1.0_dp-xi**2)*(1.0_dp+eta)*(1.0_dp-zeta)*0.25_dp
        shapefun(12) = (1.0_dp-xi)*(1.0_dp-eta**2)*(1.0_dp-zeta)*0.25_dp
        shapefun(13) = (1.0_dp-xi**2)*(1.0_dp-eta)*(1.0_dp+zeta)*0.25_dp
        shapefun(14) = (1.0_dp+xi)*(1.0_dp-eta**2)*(1.0_dp+zeta)*0.25_dp
        shapefun(15) = (1.0_dp-xi**2)*(1.0_dp+eta)*(1.0_dp+zeta)*0.25_dp
        shapefun(16) = (1.0_dp-xi)*(1.0_dp-eta**2)*(1.0_dp+zeta)*0.25_dp
        shapefun(17) = (1.0_dp-xi)*(1.0_dp-eta)*(1.0_dp-zeta**2)*0.25_dp
        shapefun(18) = (1.0_dp+xi)*(1.0_dp-eta)*(1.0_dp-zeta**2)*0.25_dp
        shapefun(19) = (1.0_dp+xi)*(1.0_dp+eta)*(1.0_dp-zeta**2)*0.25_dp
        shapefun(20) = (1.0_dp-xi)*(1.0_dp+eta)*(1.0_dp-zeta**2)*0.25_dp
     case default
        write(*,*) "shape_function ERROR: Non existing shape function type!"
  end select
  shape_function_3d=shapefun
  end function shape_function_3d

  function shape_function_der_2d(xi,eta,n,d)
  implicit none
  integer, intent(in) :: d,n
  real(dp), intent(in) :: xi,eta

  real(dp),dimension(d,n) :: shape_function_der_2d
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
  shape_function_der_2d=dhdxi
  end function shape_function_der_2d

  function shape_function_der_3d(xi,eta,zeta,n,d)
  implicit none
  integer, intent(in) :: d,n
  real(dp), intent(in) :: xi,eta,zeta

  real(dp),dimension(d,n) :: shape_function_der_3d
  real(dp),dimension(d,n) :: dhdxi
  real(dp) :: xi4

  select case(d)
    case(3)
      select case(n)
        case(10) ! Tet10 element
          xi4 = 1.-xi-eta-zeta
          dhdxi(1,1) = (4.*xi-1.)
          dhdxi(2,2) = (4.*eta-1.)
          dhdxi(3,3) = (4.*zeta-1.)
          dhdxi(1,4) = -(4.*xi4-1.)
          dhdxi(2,4) = -(4.*xi4-1.)
          dhdxi(3,4) = -(4.*xi4-1.)
          dhdxi(1,5) = 4.*eta
          dhdxi(2,5) = 4.*xi
          dhdxi(2,6) = 4.*zeta
          dhdxi(3,6) = 4.*eta
          dhdxi(1,7) = 4.*zeta
          dhdxi(3,7) = 4.*xi 
          dhdxi(1,8) = 4.*(xi4-xi)
          dhdxi(2,8) = -4.*xi
          dhdxi(3,8) = -4.*xi
          dhdxi(1,9) = -4.*eta
          dhdxi(2,9) = 4.*(xi4-eta)
          dhdxi(3,9) = -4.*eta
          dhdxi(1,10) = -4.*zeta*xi4
          dhdxi(2,10) = -4.*zeta
          dhdxi(3,10) = 4.*(xi4-zeta)
        case(20) ! Brick20 element
          dhdxi(1,1) = (-(1.-eta)*(1.-zeta)*(-xi-eta-zeta-2.)-(1.-xi)*(1.-eta)*(1.-zeta))*0.125_dp
          dhdxi(2,1) = (-(1.-xi)*(1.-zeta)*(-xi-eta-zeta-2.)-(1.-xi)*(1.-eta)*(1.-zeta))*0.125_dp
          dhdxi(3,1) = (-(1.-xi)*(1.-eta)*(-xi-eta-zeta-2.)-(1.-xi)*(1.-eta)*(1.-zeta))*0.125_dp

          dhdxi(1,2) = ((1.-eta)*(1.-zeta)*(xi-eta-zeta-2.)+(1.+xi)*(1.-eta)*(1.-zeta))*0.125_dp
          dhdxi(2,2) = (-(1.+xi)*(1.-zeta)*(xi-eta-zeta-2.)-(1.+xi)*(1.-eta)*(1.-zeta))*0.125_dp
          dhdxi(3,2) = (-(1.+xi)*(1.-eta)*(xi-eta-zeta-2.)-(1.+xi)*(1.-eta)*(1.-zeta))*0.125_dp

          dhdxi(1,3) = ((1.+eta)*(1.-zeta)*(xi+eta-zeta-2.)+(1.+xi)*(1.+eta)*(1.-zeta))*0.125_dp
          dhdxi(2,3) = ((1.+xi)*(1.-zeta)*(xi+eta-zeta-2.)+(1.+xi)*(1.+eta)*(1.-zeta))*0.125_dp
          dhdxi(3,3) = (-(1.+xi)*(1.+eta)*(xi+eta-zeta-2.)-(1.+xi)*(1.+eta)*(1.-zeta))*0.125_dp

          dhdxi(1,4) = (-(1.+eta)*(1.-zeta)*(-xi+eta-zeta-2.)-(1.-xi)*(1.+eta)*(1.-zeta))*0.125_dp
          dhdxi(2,4) = ((1.-xi)*(1.-zeta)*(-xi+eta-zeta-2.)+(1.-xi)*(1.+eta)*(1.-zeta))*0.125_dp
          dhdxi(3,4) = (-(1.-xi)*(1.+eta)*(-xi+eta-zeta-2.)-(1.-xi)*(1.+eta)*(1.-zeta))*0.125_dp

          dhdxi(1,5) = (-(1.-eta)*(1.+zeta)*(-xi-eta+zeta-2.)-(1.-xi)*(1.-eta)*(1.+zeta))*0.125_dp
          dhdxi(2,5) = (-(1.-xi)*(1.+zeta)*(-xi-eta+zeta-2.)-(1.-xi)*(1.-eta)*(1.+zeta))*0.125_dp
          dhdxi(3,5) = ((1.-xi)*(1.-eta)*(-xi-eta+zeta-2.)+(1.-xi)*(1.-eta)*(1.+zeta))*0.125_dp

          dhdxi(1,6) = ((1.-eta)*(1.+zeta)*(xi-eta+zeta-2.)+(1.+xi)*(1.-eta)*(1.+zeta))*0.125_dp
          dhdxi(2,6) = (-(1.+xi)*(1.+zeta)*(xi-eta+zeta-2.)-(1.+xi)*(1.-eta)*(1.+zeta))*0.125_dp
          dhdxi(3,6) = ((1.+xi)*(1.-eta)*(xi-eta+zeta-2.)+(1.+xi)*(1.-eta)*(1.+zeta))*0.125_dp

          dhdxi(1,7) = ((1.+eta)*(1.+zeta)*(xi+eta+zeta-2.)+(1.+xi)*(1.+eta)*(1.+zeta))*0.125_dp
          dhdxi(2,7) = ((1.+xi)*(1.+zeta)*(xi+eta+zeta-2.)+(1.+xi)*(1.+eta)*(1.+zeta))*0.125_dp
          dhdxi(3,7) = ((1.+xi)*(1.+eta)*(xi+eta+zeta-2.)+(1.+xi)*(1.+eta)*(1.+zeta))*0.125_dp

          dhdxi(1,8) = (-(1.+eta)*(1.+zeta)*(-xi+eta+zeta-2.)-(1.-xi)*(1.+eta)*(1.+zeta))*0.125_dp
          dhdxi(2,8) = ((1.-xi)*(1.+zeta)*(-xi+eta+zeta-2.)+(1.-xi)*(1.+eta)*(1.+zeta))*0.125_dp
          dhdxi(3,8) = ((1.-xi)*(1.+eta)*(-xi+eta+zeta-2.)+(1.-xi)*(1.+eta)*(1.+zeta))*0.125_dp

          dhdxi(1,9)  = -2.*xi*(1.-eta)*(1.-zeta)*0.25_dp
          dhdxi(2,9)  = -(1.-xi**2)*(1.-zeta)*0.25_dp
          dhdxi(3,9)  = -(1.-xi**2)*(1.-eta)*0.25_dp

          dhdxi(1,10)  = (1.-eta**2)*(1.-zeta)*0.25_dp
          dhdxi(2,10)  = -2.*eta*(1.+xi)*(1.-zeta)*0.25_dp
          dhdxi(3,10)  = -(1.-eta**2)*(1.+xi)*0.25_dp

          dhdxi(1,11)  = -2.*xi*(1.+eta)*(1.-zeta)*0.25_dp
          dhdxi(2,11)  = (1.-xi**2)*(1.-zeta)*0.25_dp
          dhdxi(3,11)  = -(1.-xi**2)*(1.+eta)*0.25_dp

          dhdxi(1,12)  = -(1.-eta**2)*(1.-zeta)*0.25_dp
          dhdxi(2,12)  = -2.*eta*(1.-xi)*(1.-zeta)*0.25_dp
          dhdxi(3,12)  = -(1.-eta**2)*(1.-xi)*0.25_dp

          dhdxi(1,13)  = -2.*xi*(1.-eta)*(1.+zeta)*0.25_dp
          dhdxi(2,13)  = -(1.-xi**2)*(1.+zeta)*0.25_dp
          dhdxi(3,13)  = (1.-xi**2)*(1.-eta)*0.25_dp

          dhdxi(1,14)  = (1.-eta**2)*(1.+zeta)*0.25_dp
          dhdxi(2,14)  = -2.*eta*(1.+xi)*(1.+zeta)*0.25_dp
          dhdxi(3,14)  = (1.-eta**2)*(1.+xi)*0.25_dp

          dhdxi(1,15)  = -2.*xi*(1.+eta)*(1.+zeta)*0.25_dp
          dhdxi(2,15)  =  (1.-xi**2)*(1.+zeta)*0.25_dp
          dhdxi(3,15)  = (1.-xi**2)*(1.+eta)*0.25_dp

          dhdxi(1,16)  = -(1.-eta**2)*(1.+zeta)*0.25_dp
          dhdxi(2,16)  = -2.*eta*(1.-xi)*(1.+zeta)*0.25_dp
          dhdxi(3,16)  = (1.-eta**2)*(1.-xi)*0.25_dp

          dhdxi(1,17) = -(1.-eta)*(1.-zeta**2)*0.25_dp
          dhdxi(2,17) = -(1.-xi)*(1.-zeta**2)*0.25_dp
          dhdxi(3,17) = -zeta*(1.-xi)*(1.-eta)*0.5_dp

          dhdxi(1,18) = (1.-eta)*(1.-zeta**2)*0.25_dp
          dhdxi(2,18) = -(1.+xi)*(1.-zeta**2)*0.25_dp
          dhdxi(3,18) = -zeta*(1.+xi)*(1.-eta)*0.5_dp

          dhdxi(1,19) = (1.+eta)*(1.-zeta**2)*0.25_dp
          dhdxi(2,19) = (1.+xi)*(1.-zeta**2)*0.25_dp
          dhdxi(3,19) = -zeta*(1.+xi)*(1.+eta)*0.5_dp

          dhdxi(1,20) = -(1.+eta)*(1.-zeta**2)*0.25_dp
          dhdxi(2,20) = (1.-xi)*(1.-zeta**2)*0.25_dp
          dhdxi(3,20) = -zeta*(1.-xi)*(1.+eta)*0.5_dp
        case default
            write(*,*) "shape_function_der ERROR: Non existing shape function type!"
      end select   
    case default
      write(*,*) "shape_function_der ERROR: Elements of that dimension not implemented yet!"
    end select
  shape_function_der_3d=dhdxi
  end function shape_function_der_3d


  end module shape_functions
