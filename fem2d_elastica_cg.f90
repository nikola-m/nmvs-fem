  program fem2d_elastica_cg
!
! elastica pefix > output.txt; where prefix is the common filename.
!
  use prec
  use quadrature
  use fem_module
  use fem_mesh
  use output
  
  implicit none

  integer :: iel,inode,idof
  integer :: i,j,k,ni,nj
  integer :: quad_rule,lq2

  real(dp) :: E,nu,t ! Young's modulus E; Poisson coefficient nu; models' thickness t

  integer :: it_max  ! Maximum no. of iterations
  real(dp) :: tol    ! Error tolerance - criteria for termination of iterations

  integer :: num_constr

  integer :: nrdof                               ! Number of nodes where traction forces are applied
  integer, dimension(:),allocatable :: rightdof  ! DOF's where we make changes to RHS vector - traction force components, size (ndof*nrdof)
  real(dp), dimension(:),allocatable :: rhsval   ! Values of traction force components

  real(dp), dimension(:,:),allocatable :: kel  ! Element stiffness matrix

  logical, parameter :: debug = .true.

  character ( len = 255 ) element_filename
  character ( len = 255 ) node_filename
  character ( len = 255 ) input_filename
  character ( len = 255 ) solution_filename
  character ( len = 255 ) prefix

  integer :: iarg
  integer :: iargc
  integer :: ios

  integer :: num_arg


!######## Section 1: INTRO ########################################################
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  write ( *, '(a)' ) 'FEM2D_ELASTICA_CG'
  write ( *, '(a)' ) '  A version of NMVS_ELASTICA using sparse storage'
  write ( *, '(a)' ) '  and a Preconditioned Conjugate Gradient solver.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution of the Plane Elasticity Problems'
  write ( *, '(a)' ) '  in an arbitrary region, in 2 dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The finite element method,'
  write ( *, '(a)' ) '  with QUAD4 elements is used.'
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the filename prefix:'
    read ( *, '(a)', iostat = ios ) prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_ELASTICA_CG - Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, prefix )

  end if
!
!  Create the filenames.
!
  input_filename = trim ( prefix ) // '_input.dat'
  node_filename = trim ( prefix ) // '_nodes.dat'
  element_filename = trim ( prefix ) // '_elements.dat'
  solution_filename = trim ( prefix ) // '_values.plt'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node file is "' // trim ( node_filename ) // '".'
  write ( *, '(a)' ) '  Element file is "' // trim ( element_filename ) // '".'

!######## END Section 1: INTRO ####################################################

!######## Section 2: INPUT ########################################################

  open(unit=3,file=input_filename,status='old')
  rewind 3
!
! Material properties (consider reading them from input file).
!
  read(3,*) E   ! Young's modulus
  read(3,*) nu  ! Poisson's coeficient
  read(3,*) t   ! Thickness
!
! Related to FEM solution
!
  d=2                ! d=2 fr 2D problems
  nonoel=4           ! if QAUD4 nonel=4; if QUAD8 nonoel=8
  ndof=2             ! DOF's per node
  edof = nonoel*ndof ! Number of DOFs per element
  lq=2               ! Gauss rule in 1D; lq=2 for two-point Gauss quad rule
  quad_rule=lq       ! Alias for lq
  lq2=lq*lq          ! Sqare to get no. of quad points in 2D.
!
! For CG solver
!
  read(3,*) it_max ! Maximum number of iterations for iterative linear solver
  read(3,*) tol    ! Error tolerance, if it is reached, the iterations are stopped.
!
! Constrained DOFs and their constrained values.
!
  read(3,*) num_constr ! Number of constrained DOFs

  allocate ( bcdof(num_constr) )
  allocate ( bcval(num_constr) )  

  do i=1,num_constr
    read(3,*) inode,idof,bcval(i)
    bcdof(i) = (inode-1)*ndof+idof ! Global DOF number of a constrained (via BCs) DOF 
  end do
!  bcdof = reshape([1, 2, 241, 241],shape(bcdof))
!  bcval = reshape([0., 0., 0., 0.],shape(bcval))

! Traction force components
  read(3,*) nrdof
  allocate ( rightdof(nrdof*ndof) )
  allocate ( rhsval(nrdof*ndof) )

  write(*,*) 
  write(*,*) " Node v.s. Traction force components:"
  write(*,*) 
  do i=1,nrdof
    idof = (i-1)*ndof+1
    read(3,*) inode,rhsval(idof),rhsval(idof+1)
    write(*,'(i4,2x,e13.7,2x,e13.7)') inode,rhsval(idof), rhsval(idof+1)
    rightdof(idof)= (inode-1)*ndof+1
    rightdof(idof+1)= (inode-1)*ndof+2
  end do

!
! Finished reading input.
!
  rewind 3
  close(3)
!
! Some initial memory allocations
!
  allocate ( xi(quad_rule) )
  allocate ( weight(quad_rule) )
  allocate ( psi(nonoel,lq2 ) )
  allocate ( dpsi_dxi(d,nonoel,lq2) )
  allocate ( kel(edof,edof) )
  allocate ( map(edof) )

!######## END Section 2: INPUT #####################################################

!######## Section 3: MESH  #########################################################
!
! Read mesh 
! Prilagodjeno za AUTOMESH2D mesh generator. Kod JB se to drugacije radi
! Fajl se procita i subrutina izbroji redove i kolone i na taj nacin nadje 'nonome' i 'd'
!
  open(unit=4,file=node_filename,status='old')
  rewind 4
  read(4,*) nonome

  sdof = nonome*ndof
  allocate ( node_boundary(nonome) )
  allocate ( node_condition(sdof) )
  allocate ( coord(d,nonome) )
  allocate ( glono(nonoel,nonome) )

  do inode=1,nonome
    read(4,*) k,coord(1,inode), coord(2,inode)
  enddo
! Comment out next three lines if mesh is in a single .msh file
  close(4)
  open(unit=5,file=element_filename,status='old')
  rewind 5
! Replace '5' by '4' in next few lines, if mesh is in a single .msh file
  read(5,*) nel
  do iel=1,nel
      read(5,*) k,(glono(inode,iel) , inode=1,nonoel)
  end do 
  rewind 5
  close(5) 

!
! Write a short report about a mesh
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes =          ', nonome
  call r8mat_transpose_print_some ( d, nonome, coord, 1, 1, &
    d, 10, '  First 10 nodes' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Element order =            ', nonoel
  write ( *, '(a,i8)' ) '  Number of elements =       ', nel
  call i4mat_transpose_print_some ( 4, nonome, glono, &
    1, 1, 4, 10, '  First 10 elements' )
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Quadrature order =          ', quad_rule

!######## END Section 3: MESH ######################################################


!######## Section 4: BOUNDARY CONDITIONS  ##########################################

!!  Determine which nodes are boundary nodes and which have a
!!  finite element unknown.  Then set the boundary values.
!!
!  call triangulation_order3_boundary_node ( node_num, element_num, &
!    element_node, node_boundary )
!!
!!  Determine the node conditions.
!!  For now, we'll just assume all boundary nodes are Dirichlet.
!!
  node_condition(1:sdof) = 1

!  do node = 1, node_num
!    if ( node_boundary(node) ) then
!      node_condition(node) = 2
!    end if
!  end do
  node_condition(bcdof) = 2 

!######## END Section 4: BOUNDARY CONDITIONS  ######################################

!######## Section 5: SPARSE MATRIX ALLOCATION  #####################################

!
!  Determine the element neighbor array, just so we can estimate
!  the nonzeros.
!
  allocate ( element_neighbor(4,nel) )

  call tesselation_order4_neighbor_quads ( nel, glono, &
    element_neighbor )

  if ( debug ) then
  call i4mat_transpose_print_some ( 4, nel, element_neighbor, &
    1, 1, 4, 20, '  First 20 elements of element_neighbor:' )
  end if
!
!  Count the number of nonzeros.
!
  allocate ( adj_col(1:nonome*ndof+1) )

  call tesselation_order4_adj_count ( nonome, nel, ndof, glono, &
    element_neighbor, nz_num, adj_col )

  if ( debug ) then
    call i4vec_print ( 20, adj_col, '  First 20 lines of the Adjacency vector:' )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzero coefficients NZ_NUM = ', nz_num
!
!  Set up the sparse row and column index vectors.
!
  allocate ( ia(1:nz_num) )
  allocate ( ja(1:nz_num) )

  call tesselation_order4_adj_set2 ( nonome, nel, ndof, glono, &
    element_neighbor, nz_num, adj_col, ia, ja )

  if ( debug ) then
    call i4vec_print2 ( 100, ia, ja, '  First 100 lines of the IA and JA arrays:' )
  end if

  deallocate ( adj_col )
  deallocate ( element_neighbor )
!
!  Index the diagonal elements for use by the CG solver.
!
  allocate ( diag(1:sdof) ) !sdof=nonome*ndof
  k=nonome*ndof
  call diag_index ( nz_num, ia, ja, k, diag )

  if ( debug ) then
    call i4vec_print ( 10, diag, '  First 10 lines of Diagonal adjacency vector:' )
  end if
!
!  Allocate space for the coefficient matrix A and right hand side F.
!
  allocate ( a(nz_num) )
  allocate ( f(sdof) )
  allocate ( node_u(sdof) )

!######## END Section 5: SPARSE MATRIX ALLOCATION  #################################

!######## Section 6: MAIN FEM COMPUTATION  #########################################

!  
! Gaus quadrature nodes and weights - depending on quad_rule
!
  xi=gauss_pts(quad_rule)
  weight=gauss_wts(quad_rule)
!
! Shape functions and their derivatives on referent element
!
  lq=quad_rule
  psi=shape_fun_at_gauss_pts_referent_element(lq,xi)
  dpsi_dxi=shape_fun_derivs_at_gauss_pts_referent_element(lq,xi)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Assemble the finite element coefficient matrix A and the right-hand side F.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  call assemble_poisson_dsp ( node_num, node_xy, element_num, &
!    element_node, quad_num, nz_num, ia, ja, a, f )

  f(1:sdof) = 0.0_dp
  a(1:nz_num) = 0.0_dp
  do iel=1,nel   ! Loop over elements
    ! Compute element stiffness matrix
    kel=element_siff_mat(iel,dpsi_dxi,           &
                            quad_rule,xi,weight, &
                            E,nu,t,              &
                            d,nonoel,edof)
    ! Map local DOF to global
    map = elementdof(glono,nonoel,ndof,iel)
    ! Assemby
    do i=1,edof
      ni=map(i) 
      do j=1,edof
        nj=map(j)
        call dsp_ij_to_k ( nz_num, ia, ja, ni, nj, k )
        a(k) = a(k) + kel(i,j)
      end do
    end do
  enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!
! Traction force components.
! 

! rightdof=reshape([30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40],shape(rightdof))
!  rightdof = (rightdof-1)*ndof+1
  f(rightdof)=rhsval(:) 

!  rightdof=reshape([30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40],shape(rightdof))
!  rightdof = (rightdof-1)*ndof+1
!  f(rightdof)=1./dble(21)

!  rightdof=reshape([221, 223, 225, 227, 229, 231, 233, 235, 237, 239, 241],shape(rightdof))
!  f(rightdof)=-1e+5/dble(11)
!  leftdof=reshape([1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 201],shape(leftdof))
!  f(leftdof)=1e+5/dble(11)
!  f(bcdof)=bcval

!
! Print some elements of matrix A and of vector F.
!
  if ( debug ) then
    call dsp_print_some ( nonome, nonome, nz_num, ia, ja, a, 1, 1, &
      10, 10, '  Part of Finite Element matrix A:' )
    call r8vec_print_some ( nonome, f, 1, 10, &
      '  Part of right hand side vector F:' )
  end if
!
!  Adjust the linear system to account for Dirichlet boundary conditions.
!
  call dirichlet_apply_dsp ( sdof, coord, node_condition, nz_num, &
    ia, ja, a, f )

  if ( debug ) then
    call dsp_print_some ( nonome, nonome, nz_num, ia, ja, a, 1, 1, &
      10, 10, '  Part of A after adjustment for Dirichlet condition:' )
    call r8vec_print_some ( nonome, f, 1, 10, &
      '  Part of F after adjustment for Dirichlet condition:' )
  end if
!
!  Solve the linear system using the conjugate gradient method.
!
  do k = 1, nz_num
    if ( ia(k) < 1 .or. nonome*ndof < ia(k) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_ELASTICA_CG - Fatal error!'
      write ( *, '(a)' ) '  Illegal IA(K)'
      stop
    end if
    if ( ja(k) < 1 .or. nonome*ndof < ja(k) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_ELASTICA_CG - Fatal error!'
      write ( *, '(a)' ) '  Illegal JA(K)'
      stop
    end if
  end do

  call solve_cg ( sdof, diag, nz_num, ia, ja, a, f, node_u, it_max, tol )

  if ( debug .or. .true. ) then
    call r8vec_print_some ( sdof, node_u, 1, 10, &
      '  Part of the solution vector U:' )
  end if

!
! Write output file in TECPLOT, PARAVIEW format
!
  call write_tecplot(solution_filename,element_filename,nonome,nonoel,coord,node_u)

!######## END Section 6: MAIN FEM COMPUTATION  #####################################
  
!######## Section 7: FINISHING PROGAM EXECUTION ####################################

!
! Deallocate arrays
!
  deallocate ( a )
  deallocate ( diag )
  deallocate ( f )
  deallocate ( glono )
  deallocate ( ia )
  deallocate ( ja )
  deallocate ( node_boundary )
  deallocate ( node_condition )
  deallocate ( node_u )
  deallocate ( coord )
  deallocate ( bcdof )
  deallocate ( bcval )
  deallocate ( rightdof ) 
  deallocate ( rhsval ) 
  deallocate ( xi )
  deallocate ( weight )
  deallocate ( psi )
  deallocate ( dpsi_dxi )
  deallocate ( kel )
  deallocate ( map )

!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM2D_ELASTICA_CG:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

!######## END Section 7: FINISHING PROGAM EXECUTION ################################
  stop
  end program
