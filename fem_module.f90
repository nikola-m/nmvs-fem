   module fem_module
!
!  Functions for basic matrix calculations
!
   use prec

   implicit none

!..Podaci vezani za sam FEM program
   integer :: d     ! Dimension of the problem, d=2 for 2d problems, d=3 for 3d.
   integer :: ndof  ! Number of DOFs per node      
   integer :: edof  ! nonoel*ndof No. of degrees of freedom (DOF), per element. TODO Mislim da je ovo takodje broj kolona u 'data interpolation matrix'.
   integer :: sdof ! nonome*ndof
   integer, dimension(:), allocatable :: bcdof ! Constrained global DOFs. SIZE(uslovno nonobo)
   real(dp), dimension(:), allocatable :: bcval ! Values at constrained global DOFS

!  Sparse matrix data structure
   integer ( kind = 4 ) :: nz_num                          ! Number of non- zero element in sparse matrix
   real(dp), allocatable, dimension (:) :: a               ! System matrix values in COO format.
   integer ( kind = 4 ), allocatable, dimension (:) :: ja  ! Columns in COO format.
   integer ( kind = 4 ), allocatable, dimension (:) :: ia  ! Rows in COO format.
   integer ( kind = 4 ), allocatable :: diag(:)            ! Index in 'a' array where the values from the main diagonal are. MAIN DIAGONAL=a(diag)
   real(dp), allocatable, dimension (:) :: f               ! RHS vector
   real(dp), allocatable, dimension (:) :: node_u          ! Solution vector
   integer ( kind = 4 ), allocatable, dimension (:) :: adj_col ! Adjecency info - important for finding nz_num and setup of IA, JA arrays

!..END: Podaci vezani za sam FEM program

!..Podaci vezani za mrezu.
   integer :: nel               ! No. of elements in the MESH. Every element in the mesh gets one number between 1 and nel.
   integer :: nonoel            ! NO. of NOdes in the ELement. Serves for the local numeration of nodes in the element. Moze i kao parametar. Npr 4 za Q4, 8 za Q8
   integer :: nonome            ! NO. of NOdes in the MEsh. Every node in the mesh gets one number between 1 and nonome.
   integer :: nofa              ! Total number of Faces in the mesh. Denotes global address of a specific face in the mesh
   integer :: nogra             ! Broj nodova koji svaki granicni element "naslanja" na granicu, za tri elemente 2; za quad 2 za tet 3; za brick 4.
   integer :: nonobo            ! NO. of nodes that lay on the boundary of the mesh. Denotes global address of a specific boundary node in the mesh
   integer :: nofabo            ! No. of element faces that are on the boundary of the mesh. Denotes global address of a specific boundary face in the mesh.

   

   real(dp), dimension(:,:), allocatable :: coord   ! Coordinates for every node in the mesh. coord(1,43) is x-coordinate of node 43.SIZE (d,nonome). 
   integer, dimension(:,:), allocatable :: glono    ! Mapira lokalni broj DOF na globalni. Treba za asamblazu. SIZE (nonoel,nel)
   integer, dimension(:), allocatable :: map        ! Mapira lokalni broj DOF na globalni. Treba sa asamblazu. SIZE (edof=nonoel*ndof)
   integer, dimension(:,:), allocatable :: element_neighbor  ! Broj susedne celije, sa kojom se deli lice definisano i-j lokolanim nodovima. (4,nonome)
   logical, allocatable, dimension(:) :: node_boundary        ! True/False
   integer, allocatable, dimension(:) :: node_condition       ! Integer denoting condition for a specific node, takes values among 0,1,2,3. 
                                                              ! 0 - don't write eq. for this node.
                                                              ! 1 - write FEM equation for this node.
                                                              ! 2 - Dirichlet BC is setup at this node.
                                                              ! 3 - Neumann BC
                                                                   
!   integer, dimension(:,:), allocatable :: globono  ! Globalni indeks n-tog lokalnog noda na m-toj granicnom licu (ili licu za 3d). SIZE (nogra,nofabo)
!   integer, dimension(:), allocatable :: i_cond_lim ! Vezano za specificni FEM proble. Na svakom granicnom licu je definisan broj 1,2 ili 3 koji oznacava da li se postavlja Dirichlet, Neumann ili Robin BC uslov. SIZE (1,nofabo)
!..END: Podaci vezani za mrezu.  

!..Podaci potrebni za kvadrature vezane za REFERENTNI element
   real(dp), dimension(:,:), allocatable :: psi          ! Vrednost i-te bazne funckije (i=1,nonoel) u j-toj Gausovoj tacki (j=1,lq). SIZE (nonoel,lq)
   real(dp), dimension(:,:,:), allocatable :: dpsi_dxi   ! Vrednost DIFERENCIJALA u k-tom pravcu (k=1,d), i-te bazne funckije (i=1,nonoel) u j-toj Gausovoj tacki (j=1,lq). 
                                                         ! Drugim recima ovo su komponente GRAD vectora bazne funkcije u toj tacki. SIZE (d,nonoel,lq).

   real(dp), dimension(:,:), allocatable :: theta        ! Values of nf basis functions at Gauss quadrature points (\xi_1,...\xi_{lq}). SIZE (nf,lq)
   real(dp), dimension(:,:,:), allocatable :: grad_theta ! Values of d directional derivatives of nf basis functions at Gauss quadrature points (\xi_1,...\xi_{lq}). SIZE (d,nf,lq)

   integer :: lq                                         ! number of quadrature points for specific Gauss Quadrature rule.
   real(dp), dimension(:), allocatable :: weight         ! The l-th Gauss node's (l=1,lq) weight in reference element SIZE (lq).
   real(dp), dimension(:), allocatable :: xi             ! The l-th Gauss nodes. SIZE (lq)
!..END: Podaci potrebni za kvadrature vezane za REFERENTNI element 

 
   contains
!=======================================================================
   function elementdof(glono,nonoel,ndof,iel) result(map)
   implicit none
   integer, dimension(nonoel*ndof) :: map
   
   integer, intent(in) :: nonoel,ndof,iel
   integer, dimension(nonoel,nel), intent(in) :: glono
   integer :: k,i,j,start
   k=0
   do i=1,nonoel
     start = (glono(i,iel)-1)*ndof
       do j=1,ndof
         k=k+1
         map(k)=start+j
       end do
   end do
   end function elementdof
!=======================================================================
   function ref_ksi(l,m) result(ksi)
!
!  returns koordinate, l-te Gausove tacke u m-tom elementu
!
!  Result
   real(dp), dimension(d) :: ksi 
!  Input
   integer, intent(in) :: l,m
!  Locals
   integer :: k,i

   do k=1,d
     do i=1,nonoel
       ksi(k) = ksi(k) + coord(k,glono(i,m))*psi(i,l)
     enddo
   enddo

   end function ref_ksi
!=======================================================================
   function jacobian(l,m) result(jac)
!
!  Jacobian matrix at l-th Gauss quad point of m-th element
!
   implicit none
!  Result
   real(dp), dimension(d,d) :: jac
!  Input
   integer, intent(in) :: l,m
!  Locals
   integer :: k1,k2,n

!  Initialize
   jac=0.0_dp

!  Fill-in
   do k1=1,d!f
      do k2=1,d!c
       do n=1,nonoel!a
           jac(k1,k2)=jac(k1,k2)+dpsi_dxi(k1,n,l)*coord(k2,glono(n,m))
         enddo
     enddo
   enddo

   end function jacobian
!=======================================================================

   function jweight(l,m)
!
! A Jacobian value in m-th element, at a Gauss node l, times the l-th weight: \omega_l*det(J_{K_m}(\xi_l))
! poids[l,m]  u Francuskoj knjizi
!
   use matrix_module
   implicit none

!  Result
   real(dp) :: jweight
!  Input
   integer, intent(in) :: l,m

   jweight = weight(l)*det(jacobian(l,m),d)

   end function jweight

!=======================================================================

  function shape_fun_at_gauss_pts_referent_element(quad_rule,xi) result(psi)
  use shape_functions
  implicit none
  integer, intent(in) :: quad_rule
  real(dp), dimension(nonoel,quad_rule*quad_rule) :: psi
  real(dp), dimension(quad_rule),intent(in) :: xi  ! 1d-nodes and weights
! Locals
  integer :: i,j,l
  real(dp) :: ksi,eta
! Loop over Integration points
  do i=1,quad_rule
    ksi=xi(i)
    do j=1,quad_rule
      eta=xi(j)
      l=(i-1)+j
      psi(:,l) = shape_function(ksi,eta,nonoel)
    end do
  end do
  end function

!=======================================================================

  function shape_fun_derivs_at_gauss_pts_referent_element(quad_rule,xi) result(dpsi_dxi)
  use shape_functions
  implicit none
  integer, intent(in) :: quad_rule
  real(dp), dimension(d,nonoel,quad_rule*quad_rule) :: dpsi_dxi
  real(dp), dimension(quad_rule),intent(in) :: xi ! 1d-nodes
! Locals
  integer :: i,j,l
  real(dp) :: ksi,eta
! Loop over Integration points
  do i=1,quad_rule
    ksi=xi(i)
    do j=1,quad_rule
      eta=xi(j)
      l=(i-1)+j
      dpsi_dxi(:,:,l) = shape_function_der(ksi,eta,nonoel,d) 
    end do
  end do 
  end function

!########################################################################

  function element_siff_mat(iel,dpsi_dxi, &
                            quad_rule,xi,weight, &
                            E,nu,t, &
                            d,nonoel,edof) result(kel)
!
! Returns element stiffness matrix of size edof*edof
!------------------------------------------------------------------------
!
  use matrix_module 
  implicit none
! Input
  integer, intent(in) :: quad_rule,d,nonoel,edof
  integer, intent(in) :: iel ! element no.
!  real(dp), dimension(nonoel,quad_rule*quad_rule) :: psi
  real(dp), dimension(d,nonoel,quad_rule*quad_rule) :: dpsi_dxi
  real(dp), dimension(quad_rule),intent(in) :: xi,weight ! 1d-nodes and weights
  real(dp), intent(in) :: E, nu ! Coefs for constitutive relation
  real(dp), intent(in) :: t ! Thickness of the material

! Result
  real(dp), dimension(edof,edof) :: kel

! Locals
  integer :: i,j,ii,jj,l
  real(dp) :: m1,m2 ! Coefs to simplify constitutive relation
  real(dp) :: ksi,eta,wtx,wty
  real(dp), dimension(d,d) :: jac, invjac
  real(dp) :: detjac,poid
  real(dp), dimension(d,nonoel) :: dphi_dx
  real(dp), dimension(8,8) :: temp


  m1=(1-nu)/2; m2=E/(1-nu**2) ! Coefs for Stress-strain matrix

  kel=0.0_dp ! Initialize
! Loop over Integration points
  do i=1,quad_rule
    ksi=xi(i)
    wtx=weight(i)
    do j=1,quad_rule
      eta=xi(j)
      wty=weight(j)

      l=(i-1)+j ! No. of integration point (l=1,lq)

! Jacobian computation for every element (need it to differentiate in physical coordinates)
      jac=jacobian(l,iel)
      detjac = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
      invjac=0.0_dp
      invjac(1,1)= jac(2,2)/detjac
      invjac(1,2)=-jac(1,2)/detjac
      invjac(2,1)=-jac(2,1)/detjac
      invjac(2,2)= jac(1,1)/detjac

! shape function derivatives in physical coordinates
      do ii=1,nonoel ! Number of shape functions, or unknown nodal values in element
       do jj=1,d     ! Directions
          dphi_dx(jj,ii)=invjac(jj,1)*dpsi_dxi(1,ii,l)+invjac(jj,2)*dpsi_dxi(2,ii,l)
        end do
      end do

! dphi_dx(2,3)==N3*ddy, dakle prvo je pravac diferenciranja
      temp(1,1) = dphi_dx(1,1)*dphi_dx(1,1) + dphi_dx(2,1)*m1*dphi_dx(2,1)
      temp(1,2) = dphi_dx(1,1)*m1*dphi_dx(2,1) + dphi_dx(2,1)*nu*dphi_dx(1,1)
      temp(1,3) = dphi_dx(1,2)*dphi_dx(1,1) + dphi_dx(2,2)*m1*dphi_dx(2,1)
      temp(1,4) = dphi_dx(1,2)*m1*dphi_dx(2,1) + dphi_dx(2,2)*nu*dphi_dx(1,1)
      temp(1,5) = dphi_dx(1,3)*dphi_dx(1,1) + dphi_dx(2,3)*m1*dphi_dx(2,1)
      temp(1,6) = dphi_dx(1,3)*m1*dphi_dx(2,1) + dphi_dx(2,3)*nu*dphi_dx(1,1)
      temp(1,7) = dphi_dx(1,4)*dphi_dx(1,1) + dphi_dx(2,4)*m1*dphi_dx(2,1)
      temp(1,8) = dphi_dx(1,4)*m1*dphi_dx(2,1) + dphi_dx(2,4)*nu*dphi_dx(1,1)


      temp(2,1) = dphi_dx(2,1)*m1*dphi_dx(1,1) + dphi_dx(1,1)*nu*dphi_dx(2,1)
      temp(2,2) = dphi_dx(2,1)*dphi_dx(2,1) + dphi_dx(1,1)*m1*dphi_dx(1,1)
      temp(2,3) = dphi_dx(2,2)*m1*dphi_dx(1,1) + dphi_dx(1,2)*nu*dphi_dx(2,1)
      temp(2,4) = dphi_dx(2,2)*dphi_dx(2,1) + dphi_dx(1,2)*m1*dphi_dx(1,1)
      temp(2,5) = dphi_dx(2,3)*m1*dphi_dx(1,1) + dphi_dx(1,3)*nu*dphi_dx(2,1)
      temp(2,6) = dphi_dx(2,3)*dphi_dx(2,1) + dphi_dx(1,3)*m1*dphi_dx(1,1)
      temp(2,7) = dphi_dx(2,4)*m1*dphi_dx(1,1) + dphi_dx(1,4)*nu*dphi_dx(2,1)
      temp(2,8) = dphi_dx(2,4)*dphi_dx(2,1) + dphi_dx(1,4)*m1*dphi_dx(1,1)

      temp(3,1) = dphi_dx(1,1)*dphi_dx(1,2) + dphi_dx(2,1)*m1*dphi_dx(2,2)
      temp(3,2) = dphi_dx(1,1)*m1*dphi_dx(2,2) + dphi_dx(2,1)*nu*dphi_dx(1,2)
      temp(3,3) = dphi_dx(1,2)*dphi_dx(1,2) + dphi_dx(2,2)*m1*dphi_dx(2,2)
      temp(3,4) = dphi_dx(1,2)*m1*dphi_dx(2,2) + dphi_dx(2,2)*nu*dphi_dx(1,2)
      temp(3,5) = dphi_dx(1,3)*dphi_dx(1,2) + dphi_dx(2,3)*m1*dphi_dx(2,2)
      temp(3,6) = dphi_dx(1,3)*m1*dphi_dx(2,2) + dphi_dx(2,3)*nu*dphi_dx(1,2)
      temp(3,7) = dphi_dx(1,4)*dphi_dx(1,2) + dphi_dx(2,4)*m1*dphi_dx(2,2)
      temp(3,8) = dphi_dx(1,4)*m1*dphi_dx(2,2) + dphi_dx(2,4)*nu*dphi_dx(1,2)

      temp(4,1) = dphi_dx(2,1)*m1*dphi_dx(1,2) + dphi_dx(1,1)*nu*dphi_dx(2,2)
      temp(4,2) = dphi_dx(2,1)*dphi_dx(2,2) + dphi_dx(1,1)*m1*dphi_dx(1,2)
      temp(4,3) = dphi_dx(2,2)*m1*dphi_dx(1,2) + dphi_dx(1,2)*nu*dphi_dx(2,2)
      temp(4,4) = dphi_dx(2,2)*dphi_dx(2,2) + dphi_dx(1,2)*m1*dphi_dx(1,2)
      temp(4,5) = dphi_dx(2,3)*m1*dphi_dx(1,2) + dphi_dx(1,3)*nu*dphi_dx(2,2)
      temp(4,6) = dphi_dx(2,3)*dphi_dx(2,2) + dphi_dx(1,3)*m1*dphi_dx(1,2)
      temp(4,7) = dphi_dx(2,4)*m1*dphi_dx(1,2) + dphi_dx(1,4)*nu*dphi_dx(2,2)
      temp(4,8) = dphi_dx(2,4)*dphi_dx(2,2) + dphi_dx(1,4)*m1*dphi_dx(1,2)

      temp(5,1) = dphi_dx(1,1)*dphi_dx(1,3) + dphi_dx(2,1)*m1*dphi_dx(2,3)
      temp(5,2) = dphi_dx(1,1)*m1*dphi_dx(2,3) + dphi_dx(2,1)*nu*dphi_dx(1,3)
      temp(5,3) = dphi_dx(1,2)*dphi_dx(1,3) + dphi_dx(2,2)*m1*dphi_dx(2,3)
      temp(5,4) = dphi_dx(1,2)*m1*dphi_dx(2,3) + dphi_dx(2,2)*nu*dphi_dx(1,3)
      temp(5,5) = dphi_dx(1,3)*dphi_dx(1,3) + dphi_dx(2,3)*m1*dphi_dx(2,3)
      temp(5,6) = dphi_dx(1,3)*m1*dphi_dx(2,3) + dphi_dx(2,3)*nu*dphi_dx(1,3)
      temp(5,7) = dphi_dx(1,4)*dphi_dx(1,3) + dphi_dx(2,4)*m1*dphi_dx(2,3)
      temp(5,8) = dphi_dx(1,4)*m1*dphi_dx(2,3) + dphi_dx(2,4)*nu*dphi_dx(1,3)

      temp(6,1) = dphi_dx(2,1)*m1*dphi_dx(1,3) + dphi_dx(1,1)*nu*dphi_dx(2,3)
      temp(6,2) = dphi_dx(2,1)*dphi_dx(2,3) + dphi_dx(1,1)*m1*dphi_dx(1,3)
      temp(6,3) = dphi_dx(2,2)*m1*dphi_dx(1,3) + dphi_dx(1,2)*nu*dphi_dx(2,3)
      temp(6,4) = dphi_dx(2,2)*dphi_dx(2,3) + dphi_dx(1,2)*m1*dphi_dx(1,3)
      temp(6,5) = dphi_dx(2,3)*m1*dphi_dx(1,3) + dphi_dx(1,3)*nu*dphi_dx(2,3)
      temp(6,6) = dphi_dx(2,3)*dphi_dx(2,3) + dphi_dx(1,3)*m1*dphi_dx(1,3)
      temp(6,7) = dphi_dx(2,4)*m1*dphi_dx(1,3) + dphi_dx(1,4)*nu*dphi_dx(2,3)
      temp(6,8) = dphi_dx(2,4)*dphi_dx(2,3) + dphi_dx(1,4)*m1*dphi_dx(1,3)
 
      temp(7,1) = dphi_dx(1,1)*dphi_dx(1,4) + dphi_dx(2,1)*m1*dphi_dx(2,4)
      temp(7,2) = dphi_dx(1,1)*m1*dphi_dx(2,4) + dphi_dx(2,1)*nu*dphi_dx(1,4)
      temp(7,3) = dphi_dx(1,2)*dphi_dx(1,4) + dphi_dx(2,2)*m1*dphi_dx(2,4)
      temp(7,4) = dphi_dx(1,2)*m1*dphi_dx(2,4) + dphi_dx(2,2)*nu*dphi_dx(1,4)
      temp(7,5) = dphi_dx(1,3)*dphi_dx(1,4) + dphi_dx(2,3)*m1*dphi_dx(2,4)
      temp(7,6) = dphi_dx(1,3)*m1*dphi_dx(2,4) + dphi_dx(2,3)*nu*dphi_dx(1,4)
      temp(7,7) = dphi_dx(1,4)*dphi_dx(1,4) + dphi_dx(2,4)*m1*dphi_dx(2,4)
      temp(7,8) = dphi_dx(1,4)*m1*dphi_dx(2,4) + dphi_dx(2,4)*nu*dphi_dx(1,4)
 
      temp(8,1) = dphi_dx(2,1)*m1*dphi_dx(1,4) + dphi_dx(1,1)*nu*dphi_dx(2,4)
      temp(8,2) = dphi_dx(2,1)*dphi_dx(2,4) + dphi_dx(1,1)*m1*dphi_dx(1,4)
      temp(8,3) = dphi_dx(2,2)*m1*dphi_dx(1,4) + dphi_dx(1,2)*nu*dphi_dx(2,4)
      temp(8,4) = dphi_dx(2,2)*dphi_dx(2,4) + dphi_dx(1,2)*m1*dphi_dx(1,4)
      temp(8,5) = dphi_dx(2,3)*m1*dphi_dx(1,4) + dphi_dx(1,3)*nu*dphi_dx(2,4)
      temp(8,6) = dphi_dx(2,3)*dphi_dx(2,4) + dphi_dx(1,3)*m1*dphi_dx(1,4)
      temp(8,7) = dphi_dx(2,4)*m1*dphi_dx(1,4) + dphi_dx(1,4)*nu*dphi_dx(2,4)
      temp(8,8) = dphi_dx(2,4)*dphi_dx(2,4) + dphi_dx(1,4)*m1*dphi_dx(1,4)

      poid=m2*wtx*wty*detjac
      kel = kel + poid*temp

    end do
  end do
! End: Loop over Integration points
  end function element_siff_mat
!########################################################################

   function result_of_rhs_assembly(n) result(rhs) 
!
!  Assemble RHS vector. Treba: ncel, lq(number of quad points), nbf(number of basis funct)
!
   implicit none
!  Result
   real(dp), dimension(n) :: rhs
!  Input
   integer, intent(in) :: n
!  Locals
   integer :: iel,l,ni,i

   rhs=0.0_dp

!  Loop over elements
   do iel=1,nel
   !  Loop over Gauss points
      do l=1,lq
      ! Loop over basis functions
         do ni=1,nonoel
           i=glono(ni,iel)
           ! Accumulate
           !rhs(i)=rhs(i)+weight(l,iel)*phys_rhs(lq)*bi
         end do
      enddo
   enddo
   end function result_of_rhs_assembly


  end module fem_module
