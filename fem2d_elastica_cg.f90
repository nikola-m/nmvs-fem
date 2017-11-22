  program fem2d_elastica_cg
!
! elastica pefix > output.txt; where prefix is the common filename.
!
  use prec
  use quadrature
  use fem_module
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

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7,7x)') i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end

subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4mat.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 10
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7)') i
    end do

    write ( *, '(''  Row '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end

subroutine tesselation_order4_neighbor_quads ( quad_num, &
  quad_node, quad_neighbor )

!*****************************************************************************80
!
!! TESSELATION_ORDER4_NEIGHBOR_QUADS determines quad neighbors.
!
!  Discussion:
!
!    A tesselation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each quad.  However, in some cases, it is necessary to know
!    quad adjacency information, that is, which quad, if any,
!    is adjacent to a given quad on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 4* quad_NUM
!    data items.
!
!    Note that ROW is a work array allocated dynamically inside this
!    routine.  It is possible, for very large values of quad_NUM,
!    that the necessary amount of memory will not be accessible, and the
!    routine will fail.  This is a limitation of the implementation of
!    dynamic arrays in FORTRAN90.  One way to get around this would be
!    to require the user to declare ROW in the calling routine
!    as an allocatable array, get the necessary memory explicitly with
!    an ALLOCATE statement, and then pass ROW into this routine.
!
!    Of course, the point of dynamic arrays was to make it easy to
!    hide these sorts of temporary work arrays from the poor user!
!
!    This routine was revised to store the edge data in a column
!    array rather than a row array.
!
!  Example:
!
!    The input information from TRIANGLE_NODE:
!
!    quad   Nodes
!    --------   ---------------
!     1         3      4      1
!     2         3      1      2
!     3         3      2      8
!     4         2      1      5
!     5         8      2     13
!     6         8     13      9
!     7         3      8      9
!     8        13      2      5
!     9         9     13      7
!    10         7     13      5
!    11         6      7      5
!    12         9      7      6
!    13        10      9      6
!    14         6      5     12
!    15        11      6     12
!    16        10      6     11
!
!    The output information in TRIANGLE_NEIGHBOR:
!
!    Triangle Neighboring Triangles
!    --------  ---------------------
!
!     1        -1     -1      2
!     2         1      4      3
!     3         2      5      7
!     4         2     -1      8
!     5         3      8      6
!     6         5      9      7
!     7         3      6     -1
!     8         5      4     10
!     9         6     10     12
!    10         9      8     11
!    11        12     10     14
!    12         9     11     13
!    13        -1     12     16
!    14        11     -1     15
!    15        16     14     -1
!    16        13     15     -1
!
!  Licensing:
!
!    This code is disquadbuted under the GNU LGPL license.
!
!  Modified:
!
!    14 February 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) quad_NUM, the number of quads.
!
!    Input, integer ( kind = 4 ) quad_NODE(4,quad_NUM), the nodes
!    that make up each quad.
!
!    Output, integer ( kind = 4 ) quad_NEIGHBOR(4,quad_NUM), the four
!    quads that are direct neighbors of a given quad.
!    quad_NEIGHBOR(1,I) is the index of the quad which touches side 1,
!    defined by nodes 2 and 3, and so on.  quad_NEIGHBOR(1,I) is negative
!    if there is no neighbor on that side.  In this case, that side of the
!    quad lies on the boundary of the tesselation.
!
  implicit none

  integer ( kind = 4 ) quad_num
  integer ( kind = 4 ), parameter :: quad_order = 4

  integer ( kind = 4 ) col(4,4*quad_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) side1
  integer ( kind = 4 ) side2
  integer ( kind = 4 ) quad_neighbor(4,quad_num)
  integer ( kind = 4 ) quad
  integer ( kind = 4 ) quad_node(quad_order,quad_num)
  integer ( kind = 4 ) quad1
  integer ( kind = 4 ) quad2
!
!  Step 1.
!  From the list of nodes for quad T, of the form: (I,J,K)
!  construct the FOUR neighbor relations:
!
!    (I,J,1,T) or (J,I,1,T),
!    (J,K,2,T) or (K,J,2,T),
!    (K,L,3,T) or (L,K,3,T)
!    (L,I,4,T) or (I,L,4,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  do quad = 1, quad_num

    i = quad_node(1,quad)
    j = quad_node(2,quad)
    k = quad_node(3,quad)
    l = quad_node(4,quad)

    if ( i < j ) then
      col(1:4,4*(quad-1)+1) = (/ i, j, 1, quad /)
    else
      col(1:4,4*(quad-1)+1) = (/ j, i, 1, quad /)
    end if

    if ( j < k ) then
      col(1:4,4*(quad-1)+2) = (/ j, k, 2, quad /)
    else
      col(1:4,4*(quad-1)+2) = (/ k, j, 2, quad /)
    end if

    if ( k < l ) then
      col(1:4,4*(quad-1)+3) = (/ k, l, 3, quad /)
    else
      col(1:4,4*(quad-1)+3) = (/ l, k, 3, quad /)
    end if

    if ( l < i ) then
      col(1:4,4*(quad-1)+4) = (/ l, i, 4, quad /)
    else
      col(1:4,4*(quad-1)+4) = (/ i, l, 4, quad /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on rows 1 and 2; the routine we call here
!  sorts on rows 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two quads share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
!  we make sure that these two columns occur consecutively.  That will
!  make it easy to notice that the quads are neighbors.
!
  call i4col_sort_a ( 4, 4*quad_num, col )
!
!  Step 3. Neighboring quads show up as consecutive columns with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in quad_NEIGHBOR.
!
  quad_neighbor(1:4,1:quad_num) = -1

  icol = 1

  do

    if ( 4 * quad_num <= icol ) then
      exit
    end if

    if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
      icol = icol + 1
      cycle
    end if

    side1 = col(3,icol)
    quad1 = col(4,icol)
    side2 = col(3,icol+1)
    quad2 = col(4,icol+1)


    quad_neighbor(side1,quad1) = quad2
    quad_neighbor(side2,quad2) = quad1

    icol = icol + 2

  end do

  return
end

subroutine tesselation_order4_adj_count ( node_num, quad_num, ndof, &
  quad_node, quad_neighbor, adj_num, adj_col )

!*****************************************************************************80
!
!! TESSELATION_ORDER4_ADJ_COUNT counts adjacencies in a tesselation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The tesselation is assumed to involve 4-node quads.
!
!    Two nodes are "adjacent" if they are both nodes in some quad.
!    Also, a node is considered to be adjacent to itself.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid.
!
!
!    Below, we have a chart that summarizes the adjacency relationships
!    in the sample grid.  On the left, we list the node, and its neighbors,
!    with an asterisk to indicate the adjacency of the node to itself
!    (in some cases, you want to count this self adjacency and in some
!    you don't).  On the right, we list the number of adjancencies to
!    lower-indexed nodes, to the node itself, to higher-indexed nodes,
!    the total number of adjacencies for this node, and the location
!    of the first and last entries required to list this set of adjacencies
!    in a single list of all the adjacencies.
!
!    N   Adjacencies                Below  Self   Above   Total First  Last
!
!   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
!    1:  *  2  6                        0     1       2       3     1     3
!    2:  1  *  3  6  7                  1     1       3       5     4     8
!    3:  2  *  4  7  8                  1     1       3       5     9    13
!    4:  3  *  5  8  9                  1     1       3       5    14    18
!    5:  4  *  9 10                     1     1       2       4    19    22
!    6:  1  2  *  7 11                  2     1       2       5    23    27
!    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
!    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
!    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
!   10:  5  9  * 14 15                  2     1       2       5    49    53
!   11:  6  7  * 12 16                  2     1       2       5    54    58
!   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
!   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
!   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
!   15: 10 14  * 19 20                  2     1       2       5    80    84
!   16: 11 12  * 17 21                  2     1       2       5    85    89
!   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
!   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
!   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
!   20: 15 19  * 24 25                  2     1       2       5   111   115
!   21: 16 17  * 22                     2     1       1       4   116   119
!   22: 17 18 21  * 23                  3     1       1       5   120   124
!   23: 18 19 22  * 24                  3     1       1       5   125   129
!   24: 19 20 23  * 25                  3     1       1       5   130   134
!   25: 20 24  *                        2     1       0       3   135   137
!   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 August 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) quad_NUM, the number of quads.
!
!    Input, integer ( kind = 4 ) NDOF, the number of DOFs at node.
!
!    Input, integer ( kind = 4 ) quad_NODE(4,quad_NUM), lists the
!    nodes that make up each quad, in counterclockwise order.
!
!    Input, integer ( kind = 4 ) quad_NEIGHBOR(4,quad_NUM), for each
!    side of a quad, lists the neighboring quad, or -1 if there is
!    no neighbor.
!
!    Output, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Output, integer ( kind = 4 ) ADJ_COL(NODE_NUM+1).  Information about
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) quad_num
  integer ( kind = 4 ) ndof
  integer ( kind = 4 ), parameter :: quad_order = 4

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) adj_col(node_num+1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) n7
  integer ( kind = 4 ) n8
  integer ( kind = 4 ) n_adj
  integer ( kind = 4 ) quad
  integer ( kind = 4 ) quad2
  integer ( kind = 4 ) quad_neighbor(4,quad_num)
  integer ( kind = 4 ) quad_node(quad_order,quad_num)

  adj_num = 0
  n_adj = ndof
!
!  Set every node to be adjacent to itself.
!
  adj_col(1:node_num*ndof) = n_adj
!
!  Examine each quad.
!
  do quad = 1, quad_num
    ! Ovi n-ovi n1-n8 su ne brojevi nodova, vec globalnih dof, znaci ne IE nego LM
    n1 = (quad_node(1,quad)-1)*ndof + 1
    n2 = (quad_node(2,quad)-1)*ndof + 1
    n3 = (quad_node(3,quad)-1)*ndof + 1
    n4 = (quad_node(4,quad)-1)*ndof + 1

    n5 = (quad_node(1,quad)-1)*ndof + 2
    n6 = (quad_node(2,quad)-1)*ndof + 2
    n7 = (quad_node(3,quad)-1)*ndof + 2
    n8 = (quad_node(4,quad)-1)*ndof + 2
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (quad2 <= 0)
!  or if this quad is the first of the pair in which the edge
!  occurs (quad < quad2).
!
    quad2 = quad_neighbor(1,quad) ! Neighbour 1 == we share edge 1-2

    if ( quad2 < 0 .or. quad < quad2 ) then
      adj_col(n1) = adj_col(n1) + n_adj
      adj_col(n2) = adj_col(n2) + n_adj

      adj_col(n5) = adj_col(n5) + n_adj
      adj_col(n6) = adj_col(n6) + n_adj
    end if
!
!  Add edge (2,3).
!
    quad2 = quad_neighbor(2,quad) ! Neighbour 2 == we share edge 2-3

    if ( quad2 < 0 .or. quad < quad2 ) then
      adj_col(n2) = adj_col(n2) + n_adj
      adj_col(n3) = adj_col(n3) + n_adj

      adj_col(n6) = adj_col(n6) + n_adj
      adj_col(n7) = adj_col(n7) + n_adj
    end if
!
!  Add edge (3,4).
!
    quad2 = quad_neighbor(3,quad) ! Neighbour 3 == we share edge 3-4

    if ( quad2 < 0 .or. quad < quad2 ) then
      adj_col(n3) = adj_col(n3) + n_adj
      adj_col(n4) = adj_col(n4) + n_adj

      adj_col(n7) = adj_col(n7) + n_adj
      adj_col(n8) = adj_col(n8) + n_adj
    end if

!
!  Add edge (4.1).
!
    quad2 = quad_neighbor(4,quad) ! Neighbour 4 == we share edge 4-1

    if ( quad2 < 0 .or. quad < quad2 ) then
      adj_col(n1) = adj_col(n1) + n_adj
      adj_col(n4) = adj_col(n4) + n_adj

      adj_col(n5) = adj_col(n5) + n_adj
      adj_col(n8) = adj_col(n8) + n_adj
    end if

!
! Links along diagonal of the quad, every dof is linked with two other dof's
! diagonally from it on quad.
!
      adj_col(n1) = adj_col(n1) + n_adj
      adj_col(n2) = adj_col(n2) + n_adj
      adj_col(n3) = adj_col(n3) + n_adj
      adj_col(n4) = adj_col(n4) + n_adj
      adj_col(n5) = adj_col(n5) + n_adj
      adj_col(n6) = adj_col(n6) + n_adj
      adj_col(n7) = adj_col(n7) + n_adj
      adj_col(n8) = adj_col(n8) + n_adj
  end do
!
!  We used ADJ_COL to count the number of entries in each column.
!  Convert it to pointers into the ADJ array.
!
  adj_col(2:node_num*ndof+1) = adj_col(1:node_num*ndof)

  adj_col(1) = 1
  do i = 2, node_num*ndof+1
    adj_col(i) = adj_col(i-1) + adj_col(i)
  end do
!
! Number of non-zeros
!
  adj_num = adj_col(node_num*ndof+1) - 1

  return
end

subroutine tesselation_order4_adj_set2 ( node_num, quad_num, ndof, &
  quad_node, quad_neighbor, adj_num, adj_col, ia, ja )

!*****************************************************************************80
!
!! TESSELATION_ORDER4_ADJ_SET2 sets adjacencies in a tesselation.
!
!  Discussion:
!
!    This routine is called to set up the arrays IA and JA that
!    record which nodes are adjacent in a tesselation.
!
!    The tesselation is assumed to involve 4-node quads.
!
!    Two nodes are "adjacent" if they are both nodes in some quad.
!    Also, a node is considered to be adjacent to itself.
!
!    This routine can be used to set up the sparse triplet storage
!    for a linear quad finite element discretization of Poisson's
!    equation in two dimensions.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid
!
!
!    Below, we have a chart that summarizes the adjacency relationships
!    in the sample grid.  On the left, we list the node, and its neighbors,
!    with an asterisk to indicate the adjacency of the node to itself
!    (in some cases, you want to count this self adjacency and in some
!    you don't).  On the right, we list the number of adjancencies to
!    lower-indexed nodes, to the node itself, to higher-indexed nodes,
!    the total number of adjacencies for this node, and the location
!    of the first and last entries required to list this set of adjacencies
!    in a single list of all the adjacencies.
!
!    N   Adjacencies                Below  Self    Above  Total First  Last
!
!   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
!    1:  *  2  6                        0     1       2       3     1     3
!    2:  1  *  3  6  7                  1     1       3       5     4     8
!    3:  2  *  4  7  8                  1     1       3       5     9    13
!    4:  3  *  5  8  9                  1     1       3       5    14    18
!    5:  4  *  9 10                     1     1       2       4    19    22
!    6:  1  2  *  7 11                  2     1       2       5    23    27
!    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
!    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
!    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
!   10:  5  9  * 14 15                  2     1       2       5    49    53
!   11:  6  7  * 12 16                  2     1       2       5    54    58
!   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
!   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
!   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
!   15: 10 14  * 19 20                  2     1       2       5    80    84
!   16: 11 12  * 17 21                  2     1       2       5    85    89
!   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
!   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
!   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
!   20: 15 19  * 24 25                  2     1       2       5   111   115
!   21: 16 17  * 22                     2     1       1       4   116   119
!   22: 17 18 21  * 23                  3     1       1       5   120   124
!   23: 18 19 22  * 24                  3     1       1       5   125   129
!   24: 19 20 23  * 25                  3     1       1       5   130   134
!   25: 20 24  *                        2     1       0       3   135   137
!   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
!
!    For this example, the initial portion of the IA and JA arrays will be:
!
!      (1,1), (1,2), (1,6),
!      (2,1), (2,2), (2,3), (2,6), (2,7),
!      (3,2), (3,3), (3,4), (3,7), (3,8),
!      ...
!      (25,20), (25,24), (25,25)
!
!    for a total of 137 pairs of values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2007(JB); 5 November 2013(NM)
!
!  Author (of the original routine that treated triangle elements):
!
!    John Burkardt
!
!  Adaptation for QUAD elements:
! 
!     Nikola Mirkov
!
!  Parameters
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) quad_NUM, the number of quads.
!
!    Input, integer ( kind = 4 ) quad_NODE(4,quad_NUM), lists the nodes
!    that make up each quad in counterclockwise order.
!
!    Input, integer ( kind = 4 ) quad_NEIGHBOR(4,quad_NUM), for each
!    side of a quad, lists the neighboring quad, or -1 if there is
!    no neighbor.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacencies.
!
!    Input, integer ( kind = 4 ) ADJ_COL(NODE_NUM+1).  Information about
!    column J is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
!
!    Output, integer ( kind = 4 ) IA(ADJ_NUM), JA(ADJ_NUM), the adjacency
!    information.
!
  implicit none

  integer ( kind = 4 ) adj_num
  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) quad_num
  integer ( kind = 4 ) ndof
  integer ( kind = 4 ), parameter :: quad_order = 4

  integer ( kind = 4 ) adj_col(node_num*ndof+1)
  integer ( kind = 4 ) adj_copy(node_num*ndof)
  integer ( kind = 4 ) ia(adj_num)
  integer ( kind = 4 ) ja(adj_num)
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) n2
  integer ( kind = 4 ) n3
  integer ( kind = 4 ) n4
  integer ( kind = 4 ) n5
  integer ( kind = 4 ) n6
  integer ( kind = 4 ) n7
  integer ( kind = 4 ) n8
  integer ( kind = 4 ) n_adj
  integer ( kind = 4 ) node
  integer ( kind = 4 ) quad
  integer ( kind = 4 ) quad2
  integer ( kind = 4 ) quad_neighbor(4,quad_num)
  integer ( kind = 4 ) quad_node(quad_order,quad_num)

  n_adj = ndof ! When we have 'ndof' degrees of freedom at every node.

  ia(1:adj_num) = -1
  ja(1:adj_num) = -1

  adj_copy(1:node_num*ndof) = adj_col(1:node_num*ndof)
!
!  Set every {node} DOF to be adjacent to itself.
!
  do node = 1, node_num*ndof
    ia(adj_copy(node)) = node
    ja(adj_copy(node)) = node
    adj_copy(node) = adj_copy(node) + 1
  end do
!
!  Set every DOF to be adjacent to the DOF it shares the node with!
!
  do node= 1, node_num
    n1 = (node-1)*ndof + 1
    n2 = (node-1)*ndof + 2

    ia(adj_copy(n1)) = n1
    ja(adj_copy(n1)) = n2
    adj_copy(n1) = adj_copy(n1) + 1

    ia(adj_copy(n2)) = n2
    ja(adj_copy(n2)) = n1
    adj_copy(n2) = adj_copy(n2) + 1
  end do

!
!  Examine each quad.
!
  do quad = 1, quad_num

    ! Ovi n-ovi n1-n8 su ne brojevi nodova, vec globalnih dof, znaci ne IE nego LM
    n1 = (quad_node(1,quad)-1)*ndof + 1
    n2 = (quad_node(1,quad)-1)*ndof + 2
    n3 = (quad_node(2,quad)-1)*ndof + 1
    n4 = (quad_node(2,quad)-1)*ndof + 2
    n5 = (quad_node(3,quad)-1)*ndof + 1
    n6 = (quad_node(3,quad)-1)*ndof + 2
    n7 = (quad_node(4,quad)-1)*ndof + 1
    n8 = (quad_node(4,quad)-1)*ndof + 2
!
! U svakom slucaju ce sve diagonale biti povezane
!
!(1-5); (1-6)
      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n5
      adj_copy(n1) = adj_copy(n1) + 1

      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n6
      adj_copy(n1) = adj_copy(n1) + 1
!(2-5); (2-6)
      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n5
      adj_copy(n2) = adj_copy(n2) + 1

      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n6
      adj_copy(n2) = adj_copy(n2) + 1
!(5-1); (5-2)
      ia(adj_copy(n5)) = n5
      ja(adj_copy(n5)) = n1
      adj_copy(n5) = adj_copy(n5) + 1

      ia(adj_copy(n5)) = n5
      ja(adj_copy(n5)) = n2
      adj_copy(n5) = adj_copy(n5) + 1
!(6-1); (6-2)
      ia(adj_copy(n6)) = n6
      ja(adj_copy(n6)) = n1
      adj_copy(n6) = adj_copy(n6) + 1

      ia(adj_copy(n6)) = n6
      ja(adj_copy(n6)) = n2
      adj_copy(n6) = adj_copy(n6) + 1
!(3-7); (3-8)
      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n7
      adj_copy(n3) = adj_copy(n3) + 1

      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n8
      adj_copy(n3) = adj_copy(n3) + 1
!(4-7);(4-8)
      ia(adj_copy(n4)) = n4
      ja(adj_copy(n4)) = n7
      adj_copy(n4) = adj_copy(n4) + 1

      ia(adj_copy(n4)) = n4
      ja(adj_copy(n4)) = n8
      adj_copy(n4) = adj_copy(n4) + 1
!(7-3);(7-4)
      ia(adj_copy(n7)) = n7
      ja(adj_copy(n7)) = n3
      adj_copy(n7) = adj_copy(n7) + 1

      ia(adj_copy(n7)) = n7
      ja(adj_copy(n7)) = n4
      adj_copy(n7) = adj_copy(n7) + 1
!(8-3);(8-4)
      ia(adj_copy(n8)) = n8
      ja(adj_copy(n8)) = n3
      adj_copy(n8) = adj_copy(n8) + 1

      ia(adj_copy(n8)) = n8
      ja(adj_copy(n8)) = n4
      adj_copy(n8) = adj_copy(n8) + 1

!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (quad2 <= 0)
!  or if this quad is the first of the pair in which the edge
!  occurs (quad < quad2).
!
    quad2 = quad_neighbor(1,quad)

    if ( quad2 < 0 .or. quad < quad2 ) then

!      ia(adj_copy(n1)) = n1
!      ja(adj_copy(n1)) = n2
!      adj_copy(n1) = adj_copy(n1) + 1

      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1

      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n4
      adj_copy(n1) = adj_copy(n1) + 1
!--
!      ia(adj_copy(n2)) = n2
!      ja(adj_copy(n2)) = n1
!      adj_copy(n2) = adj_copy(n2) + 1
!--
      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1

      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n4
      adj_copy(n2) = adj_copy(n2) + 1
!---------------------
      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1

      ia(adj_copy(n4)) = n4
      ja(adj_copy(n4)) = n1
      adj_copy(n4) = adj_copy(n4) + 1

      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1

      ia(adj_copy(n4)) = n4
      ja(adj_copy(n4)) = n2
      adj_copy(n4) = adj_copy(n4) + 1

    end if
!
!  Add edge (2,3).
!
    quad2 = quad_neighbor(2,quad)

    if ( quad2 < 0 .or. quad < quad2 ) then

!      ia(adj_copy(n3)) = n3
!      ja(adj_copy(n3)) = n4
!      adj_copy(n3) = adj_copy(n3) + 1

      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n5
      adj_copy(n3) = adj_copy(n3) + 1

      ia(adj_copy(n3)) = n3
      ja(adj_copy(n3)) = n6
      adj_copy(n3) = adj_copy(n3) + 1
!--
!      ia(adj_copy(n4)) = n4
!      ja(adj_copy(n4)) = n3
!      adj_copy(n4) = adj_copy(n4) + 1
!--
      ia(adj_copy(n4)) = n4
      ja(adj_copy(n4)) = n5
      adj_copy(n4) = adj_copy(n4) + 1

      ia(adj_copy(n4)) = n4
      ja(adj_copy(n4)) = n6
      adj_copy(n4) = adj_copy(n4) + 1
!---------------------
      ia(adj_copy(n5)) = n5
      ja(adj_copy(n5)) = n3
      adj_copy(n5) = adj_copy(n5) + 1

      ia(adj_copy(n6)) = n6
      ja(adj_copy(n6)) = n3
      adj_copy(n6) = adj_copy(n6) + 1

      ia(adj_copy(n5)) = n5
      ja(adj_copy(n5)) = n4
      adj_copy(n5) = adj_copy(n5) + 1

      ia(adj_copy(n6)) = n6
      ja(adj_copy(n6)) = n4
      adj_copy(n6) = adj_copy(n6) + 1

    end if
!
!  Add edge (3,4).
!
    quad2 = quad_neighbor(3,quad)

    if ( quad2 < 0 .or. quad < quad2 ) then

!      ia(adj_copy(n5)) = n5
!      ja(adj_copy(n5)) = n6
!      adj_copy(n5) = adj_copy(n5) + 1

      ia(adj_copy(n5)) = n5
      ja(adj_copy(n5)) = n7
      adj_copy(n5) = adj_copy(n5) + 1

      ia(adj_copy(n5)) = n5
      ja(adj_copy(n5)) = n8
      adj_copy(n5) = adj_copy(n5) + 1
!--
!      ia(adj_copy(n6)) = n6
!      ja(adj_copy(n6)) = n5
!      adj_copy(n6) = adj_copy(n6) + 1
!--
      ia(adj_copy(n6)) = n6
      ja(adj_copy(n6)) = n7
      adj_copy(n6) = adj_copy(n6) + 1

      ia(adj_copy(n6)) = n6
      ja(adj_copy(n6)) = n8
      adj_copy(n6) = adj_copy(n6) + 1
!---------------------
      ia(adj_copy(n7)) = n7
      ja(adj_copy(n7)) = n5
      adj_copy(n7) = adj_copy(n7) + 1

      ia(adj_copy(n8)) = n8
      ja(adj_copy(n8)) = n5
      adj_copy(n8) = adj_copy(n8) + 1

      ia(adj_copy(n7)) = n7
      ja(adj_copy(n7)) = n6
      adj_copy(n7) = adj_copy(n7) + 1

      ia(adj_copy(n8)) = n8
      ja(adj_copy(n8)) = n6
      adj_copy(n8) = adj_copy(n8) + 1

    end if

!
!  Add edge (4,1).
!
    quad2 = quad_neighbor(4,quad)

    if ( quad2 < 0 .or. quad < quad2 ) then

!      ia(adj_copy(n7)) = n7
!      ja(adj_copy(n7)) = n8
!      adj_copy(n7) = adj_copy(n7) + 1

      ia(adj_copy(n7)) = n7
      ja(adj_copy(n7)) = n1
      adj_copy(n7) = adj_copy(n7) + 1

      ia(adj_copy(n7)) = n7
      ja(adj_copy(n7)) = n2
      adj_copy(n7) = adj_copy(n7) + 1
!--
!      ia(adj_copy(n8)) = n8
!      ja(adj_copy(n8)) = n7
!      adj_copy(n8) = adj_copy(n8) + 1
!--
      ia(adj_copy(n8)) = n8
      ja(adj_copy(n8)) = n1
      adj_copy(n8) = adj_copy(n8) + 1

      ia(adj_copy(n8)) = n8
      ja(adj_copy(n8)) = n2
      adj_copy(n8) = adj_copy(n8) + 1
!---------------------
      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n7
      adj_copy(n1) = adj_copy(n1) + 1

      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n7
      adj_copy(n2) = adj_copy(n2) + 1

      ia(adj_copy(n1)) = n1
      ja(adj_copy(n1)) = n8
      adj_copy(n1) = adj_copy(n1) + 1

      ia(adj_copy(n2)) = n2
      ja(adj_copy(n2)) = n8
      adj_copy(n2) = adj_copy(n2) + 1

    end if

  end do
!
!  Lexically sort the IA, JA values.
!
  call i4vec2_sort_a ( adj_num, ia, ja )

  return
end

subroutine i4col_sort_a ( m, n, a )

!*****************************************************************************80
!
!! I4COL_SORT_A ascending sorts an I4COL.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two real
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, the first time they differ, X is smaller.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
!    a vector of data.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, integer ( kind = 4 ) A(M,N).
!    On input, the array of N columns of M-vectors.
!    On output, the columns of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  if ( m <= 0 ) then
    return
  end if

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call i4col_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4col_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine i4col_swap ( m, n, a, i, j )

!*****************************************************************************80
!
!! I4COL_SWAP swaps columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      A = (
!        1  4  3  2
!        5  8  7  6
!        9 12 11 10 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns of
!    length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be swapped.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) col(m)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
    write ( *, '(a)' ) '  I or J is out of bounds.'
    write ( *, '(a,i8)' ) '  I =    ', i
    write ( *, '(a,i8)' ) '  J =    ', j
    write ( *, '(a,i8)' ) '  N =    ', n
    stop

  end if

  if ( i == j ) then
    return
  end if

  col(1:m) = a(1:m,i)
  a(1:m,i) = a(1:m,j)
  a(1:m,j) = col(1:m)

  return
end

subroutine i4col_compare ( m, n, a, i, j, isgn )

!*****************************************************************************80
!
!! I4COL_COMPARE compares columns I and J of an I4COL.
!
!  Example:
!
!    Input:
!
!      M = 3, N = 4, I = 2, J = 4
!
!      A = (
!        1  2  3  4
!        5  6  7  8
!        9 10 11 12 )
!
!    Output:
!
!      ISGN = -1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors
!    of length M.
!
!    Input, integer ( kind = 4 ) I, J, the columns to be compared.
!    I and J must be between 1 and N.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, column I < column J,
!     0, column I = column J,
!    +1, column J < column I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
!
!  Check.
!
  if ( i < 1 .or. n < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index I is out of bounds.'
    write ( *, '(a,i8)' ) '  I = ', i
    stop
  end if

  if ( j < 1 .or. n < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Column index J is out of bounds.'
    write ( *, '(a,i8)' ) '  J = ', j
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= m )

    if ( a(k,i) < a(k,j) ) then
      isgn = -1
      return
    else if ( a(k,j) < a(k,i) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end

subroutine i4vec2_sort_a ( n, a1, a2 )

!*****************************************************************************80
!
!! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
!
!  Discussion:
!
!    Each item to be sorted is a pair of integers (I,J), with the I
!    and J values stored in separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items of data.
!
!    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ) temp

  if ( n <= 1 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( n, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      temp  = a1(i)
      a1(i) = a1(j)
      a1(j) = temp

      temp  = a2(i)
      a2(i) = a2(j)
      a2(j) = temp
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call i4vec2_compare ( n, a1, a2, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end

subroutine sort_heap_external ( n, indx, i, j, isgn )

!*****************************************************************************80
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2004
!
!  Author:
!
!    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of items to be sorted.
!
!    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer ( kind = 4 ) I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ( kind = 4 ) ISGN, results of comparison of elements
!    I and J. (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: i_save = 0
  integer ( kind = 4 ) indx
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  integer ( kind = 4 ), save :: j_save = 0
  integer ( kind = 4 ), save :: k = 0
  integer ( kind = 4 ), save :: k1 = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ), save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end

subroutine diag_index ( m, ia, ja, n, diag )

!*****************************************************************************80
!
!! DIAG_INDEX determines where the diagonal matrix entries are stored.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of adjacencies.
!
!    Input, integer ( kind = 4 ) IA(M), JA(M), the row and column indices 
!    of adjacencies.
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Output, integer ( kind = 4 ) DIAG(N), contains for each index 1 <= I <= N, 
!    the unique index J such that IA[J] = JA[J] = I.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(m)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(m)

  diag(1:n) = -1

  do j = 1, m
    if ( ia(j) == ja(j) ) then
      i = ia(j)
      if ( diag(i) /= -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIAG_INDEX - Fatal error!'
        write ( *, '(a)' ) '  Multiple occurrences of diagonal pairs.'
        write ( *, '(a,i6,a,i6,a,i6,a)' ) &
          '  IA(', j, ') = JA(', j, ') = ', ia(j), ' and'
        write ( *, '(a,i6,a,i6,a,i6)'   ) &
          '  IA(', diag(i), ') = JA(', diag(i), ') = ', ia(j)
        stop
      end if
      diag(i) = j
    end if
  end do

  do i = 1,  n
    if ( diag(i) == -1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIAG_INDEX - Fatal error!'
      write ( *, '(a,i6,a)' ) '  DIAG(', i, ') = -1.'
      stop
    end if
  end do

  return
end

subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

!*****************************************************************************80
!
!! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data items.
!
!    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
!    of each item.
!
!    Input, integer ( kind = 4 ) I, J, the items to be compared.
!
!    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
!    -1, item I < item J,
!     0, item I = item J,
!    +1, item J < item I.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a1(n)
  integer ( kind = 4 ) a2(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j

  isgn = 0

       if ( a1(i) < a1(j) ) then

    isgn = -1

  else if ( a1(i) == a1(j) ) then

         if ( a2(i) < a2(j) ) then
      isgn = -1
    else if ( a2(i) < a2(j) ) then
      isgn = 0
    else if ( a2(j) < a2(i) ) then
      isgn = +1
    end if

  else if ( a1(j) < a1(i) ) then

    isgn = +1

  end if

  return
end

subroutine i4vec_print ( n, a, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i12)' ) i, ':', a(i)
  end do

  return
end

subroutine i4vec_print2 ( n, a, b, title )

!*****************************************************************************80
!
!! I4VEC_PRINT prints an I4VEC.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 May 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, integer ( kind = 4 ) A(N), first vector to be printed.
!
!    Input, integer ( kind = 4 ) B(N), second vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(n),b(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,a,2x,i4,2x,i4)' ) i, ':', a(i), b(i)
  end do

  return
end

subroutine dsp_ij_to_k ( nz_num, row, col, i, j, k )

!*****************************************************************************80
!
!! DSP_IJ_TO_K seeks the compressed index of the (I,J) entry of A.
!
!  Discussion:
!
!    If A(I,J) is nonzero, then its value is stored in location K.
!
!    This routine searches the DSP storage structure for the index K
!    corresponding to (I,J), returning -1 if no such entry was found.
!
!    This routine assumes that the data structure has been sorted,
!    so that the entries of ROW are ascending sorted, and that the
!    entries of COL are ascending sorted, within the group of entries
!    that have a common value of ROW.
!
!    The DSP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    The DSP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements in
!    the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and
!    column indices of the nonzero elements.
!
!    Input, integer ( kind = 4 ) I, J, the row and column indices of the
!    matrix entry.
!
!    Output, integer ( kind = 4 ) K, the DSP index of the (I,J) entry.
!
  implicit none

  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) md
  integer ( kind = 4 ) row(nz_num)

  lo = 1
  hi = nz_num

  do

    if ( hi < lo ) then
      k = -1
      exit
    end if

    md = ( lo + hi ) / 2

    if ( row(md) < i .or. ( row(md) == i .and. col(md) < j ) ) then
      lo = md + 1
    else if ( i < row(md) .or. ( row(md) == i .and. j < col(md) ) ) then
      hi = md - 1
    else
      k = md
      exit
    end if

  end do

  return
end

subroutine dsp_print_some ( m, n, nz_num, row, col, a, ilo, jlo, &
  ihi, jhi, title )

!*****************************************************************************80
!
!! DSP_PRINT_SOME prints some of a DSP matrix.
!
!  Discussion:
!
!    This version of DSP_PRINT_SOME has been specifically modified to allow,
!    and correctly handle, the case in which a single matrix location
!    A(I,J) is referenced more than once by the sparse matrix structure.
!    In such cases, the routine prints out the sum of all the values.
!
!    The DSP storage format stores the row, column and value of each nonzero
!    entry of a sparse matrix.
!
!    It is possible that a pair of indices (I,J) may occur more than
!    once.  Presumably, in this case, the intent is that the actual value
!    of A(I,J) is the sum of all such entries.  This is not a good thing
!    to do, but I seem to have come across this in MATLAB.
!
!    The DSP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
!    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 September 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    of the matrix.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero elements
!    in the matrix.
!
!    Input, integer ( kind = 4 ) ROW(NZ_NUM), COL(NZ_NUM), the row and column
!    indices of the nonzero elements.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero elements of the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title to print.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) aij(incx)
  integer ( kind = 4 ) col(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) row(nz_num)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if
!
!  Print the columns of the matrix, in strips of 5.
!
  do j2lo = jlo, jhi, incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '
    write ( *, '(''  Col:  '',5(i7,7x))' ) ( j, j = j2lo, j2hi )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
      aij(1:inc) = 0.0D+00
!
!  Is matrix entry K actually the value of A(I,J), with J2LO <= J <= J2HI?
!  Because MATLAB seems to allow for multiple (I,J,A) entries, we have
!  to sum up what we find.
!
      do k = 1, nz_num

        if ( i == row(k) .and. &
             j2lo <= col(k) .and. &
             col(k) <= j2hi ) then

          j2 = col(k) - j2lo + 1
          aij(j2) = aij(j2) + a(k)

        end if

      end do

      if ( any ( aij(1:inc) /= 0.0D+00 ) ) then
        write ( *, '(i5,1x,5g14.6)' ) i, aij(1:inc)
      end if

    end do

  end do

  return
end

subroutine r8vec_print_some ( n, a, i_lo, i_hi, title )

!*****************************************************************************80
!
!! R8VEC_PRINT_SOME prints "some" of an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 October 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) I_LO, I_HI, the first and last indices
!    to print.  The routine expects 1 <= I_LO <= I_HI <= N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i_hi
  integer ( kind = 4 ) i_lo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = max ( i_lo, 1 ), min ( i_hi, n )
    write ( *, '(2x,i8,2x,g14.8)' ) i, a(i)
  end do

  return
end

subroutine dirichlet_apply_dsp ( node_num, node_xy, node_condition, &
  nz_num, ia, ja, a, f )

!*****************************************************************************80
!
!! DIRICHLET_APPLY_DSP accounts for Dirichlet boundary conditions.
!
!  Discussion:
!
!    It is assumed that the matrix A and right hand side F have already been
!    set up as though there were no boundary conditions.  This routine
!    then modifies A and F, essentially replacing the finite element equation
!    at a boundary node NODE by a trivial equation of the form
!
!      A(NODE,NODE) * U(NODE) = NODE_BC(NODE)
!
!    where A(NODE,NODE) = 1.
!
!    This routine assumes that the coefficient matrix is stored in a
!    sparse triplet format.
!
!    This routine implicitly assumes that the sparse matrix has a storage
!    location for every diagonal element...or at least for those diagonal
!    elements corresponding to boundary nodes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of nodes.
!
!    Input, integer ( kind = 4 ) NODE_CONDITION(NODE_NUM), reports the
!    condition used to set the unknown associated with the node.
!    0, unknown.
!    1, finite element equation.
!    2, Dirichlet condition;
!    3, Neumann condition.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the nonzero entries.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the coefficient matrix,
!    stored in sparse triplet format; on output, the matrix has been adjusted
!    for Dirichlet boundary conditions.
!
!    Input/output, real ( kind = 8 ) F(NODE_NUM), the right hand side.
!    On output, the right hand side has been adjusted for Dirichlet
!    boundary conditions.
!
  implicit none

  integer ( kind = 4 ) node_num
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ), dimension(nz_num) :: a
  integer ( kind = 4 ) column
  integer ( kind = 4 ), parameter :: DIRICHLET = 2
  real ( kind = 8 ), dimension(node_num) :: f
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) node
  real ( kind = 8 ), dimension ( node_num ) :: node_bc
  integer ( kind = 4 ) node_condition(node_num)
  real ( kind = 8 ), dimension(2,node_num) :: node_xy
  integer ( kind = 4 ) nz
!
!  TODO Retrieve the Dirichlet boundary condition value at every node. 
!
!@#  call dirichlet_condition ( node_num, node_xy, node_bc ) ! SUPPLY THIS FUNCTION
!
!  Consider every matrix entry, NZ.
!
!  If the row I corresponds to a boundary node, then
!  zero out all off diagonal matrix entries, set the diagonal to 1,
!  and the right hand side to the Dirichlet boundary condition value.
!
  do nz = 1, nz_num

    node = ia(nz)

    if ( node_condition(node) == DIRICHLET ) then

      column = ja(nz)

      if ( column == node ) then
        a(nz) = 1.0D+00
        f(node) = 0.0D0 ! node_bc(node) FIXME
      else
        a(nz) = 0.0D+00
      end if

    end if

  end do

  return
end

subroutine solve_cg ( n, diag, nz_num, ia, ja, a, b, x, it_max, tol )

!*****************************************************************************80
!
!! SOLVE_CG solves a linear system using the conjugate gradient method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes.
!
!    Input, integer ( kind = 4 ) DIAG(N), contains for each index 1 <= I <= N, 
!    the unique index J such that IA(J) = JA(J) = I.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzero entries.
!
!    Input, integer ( kind = 4 ) IA(NZ_NUM), JA(NZ_NUM), the row and column
!    indices of the nonzero entries.
!
!    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero entries of the matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) bnrm2
  integer ( kind = 4 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nz_num)
  integer ( kind = 4 ) it
  integer ( kind = 4 ) it_max
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) job
  integer ( kind = 4 ) k
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) q(n)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rnrm2
  real ( kind = 8 ) tol
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) z(n)

  it = 0
  bnrm2 = sqrt ( sum ( b(1:n)**2 ) )

  x(1:n) = b(1:n) / a(diag(1:n))

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Step        Residual'
  write ( *, '(a)' ) ''

  job = 1

  do
  
    call cg_rc ( n, b, x, r, z, p, q, job )

    if ( job == 1 ) then
!   Compute q = A * p.

      q(1:n) = 0.0D+00
      do k = 1, nz_num
        i = ia(k)
        j = ja(k)
        q(i) = q(i) + a(k) * p(j)
      end do

    else if ( job == 2 ) then
!   Solve M * z = r. NOTE: HERE WE USE DIAGONAL PRECONDITIONER!
      z(1:n) = r(1:n) / a(diag(1:n))

    else if ( job == 3 ) then
!   Residual update. Compute r = r - A * x. 
      do k = 1, nz_num
        i = ia(k)
        j = ja(k)
        r(i) = r(i) - a(k) * x(j)
      end do

    else if ( job == 4 ) then
!   Perform stopping test.

      rnrm2 = sqrt ( sum ( r(1:n)**2 ) )

      if ( bnrm2 == 0.0D+00 ) then
        if ( rnrm2 <= tol ) then
          exit
        end if
      else
        if ( rnrm2 <= tol * bnrm2 ) then
          exit
        end if
      end if

      it = it + 1
      write ( *, '(2x,i4,2x,g14.6)' ) it, rnrm2

      if ( it_max <= it ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a)' ) '  Terminating early.'
        exit
      end if

    end if

    job = 2

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of iterations was ', it
  write ( *, '(a,g14.6)' ) '  Estimated error is ', rnrm2

  return
end

subroutine cg_rc ( n, b, x, r, z, p, q, job )

!*****************************************************************************80
!
!! CG_RC is a reverse communication conjugate gradient routine.
!
!  Discussion:
!
!    This routine seeks a solution of the linear system A*x=b
!    where b is a given right hand side vector, A is an n by n
!    symmetric positive definite matrix, and x is an unknown vector
!    to be determined.
!
!    Under the assumptions that the matrix A is large and sparse,
!    the conjugate gradient method may provide a solution when
!    a direct approach would be impractical because of excessive
!    requirements of storage or even of time.
!
!    The conjugate gradient method presented here does not require the 
!    user to store the matrix A in a particular way.  Instead, it only 
!    supposes that the user has a way of calculating
!      y = alpha * A * x + b * y
!    and of solving the preconditioned linear system
!      M * x = b
!    where M is some preconditioning matrix, which might be merely
!    the identity matrix, or a diagonal matrix containing the
!    diagonal entries of A.
!
!    This routine was extracted from the "templates" package.
!    There, it was not intended for direct access by a user;
!    instead, a higher routine called "cg()" was called once by
!    the user.  The cg() routine then made repeated calls to 
!    cgrevcom() before returning the result to the user.
!
!    The reverse communication feature of cgrevcom() makes it, by itself,
!    a very powerful function.  It allows the user to handle issues of
!    storage and implementation that would otherwise have to be
!    mediated in a fixed way by the function argument list.  Therefore,
!    this version of cgrecom() has been extracted from the templates
!    library and documented as a stand-alone procedure.
!
!    The user sets the value of JOB to 1 before the first call,
!    indicating the beginning of the computation, and to the value of
!    2 thereafter, indicating a continuation call.  
!    The output value of JOB is set by cgrevcom(), which
!    will return with an output value of JOB that requests a particular
!    new action from the user.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994,
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
!    Input/output, real ( kind = 8 ) X(N).  On first call, the user 
!    should store an initial guess for the solution in X.  On return with
!    JOB = 4, X contains the latest solution estimate.
!
!    Input/output, real ( kind = 8 ) R(N), Z(N), P(N), Q(N),
!    information used by the program during the calculation.  The user
!    does not need to initialize these vectors.  However, specific
!    return values of JOB may require the user to carry out some computation
!    using data in some of these vectors.
!
!    Input/output, integer ( kind = 4 ) JOB, communicates the task to be done.
!    The user needs to set the input value of JOB to 1, before the first call,
!    and then to 2 for every subsequent call for the given problem.
!    The output value of JOB indicates the requested user action.  
!    * JOB = 1, compute Q = A * P;
!    * JOB = 2: solve M*Z=R, where M is the preconditioning matrix;
!    * JOB = 3: compute R = R - A * X;
!    * JOB = 4: check the residual R for convergence.  
!               If satisfactory, terminate the iteration.
!               If too many iterations were taken, terminate the iteration.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(n)
  real ( kind = 8 ) beta
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) job
  real ( kind = 8 ) p(n)
  real ( kind = 8 ) pdotq
  real ( kind = 8 ) q(n)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_old
  integer ( kind = 4 ) rlbl
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) z(n)
!
!  Some local variables must be preserved between calls.
!
  save iter
  save rho
  save rho_old
  save rlbl
!
!  Initialization.
!  Ask the user to compute the initial residual.
!
  if ( job == 1 ) then

    r(1:n) = b(1:n)

    job = 3
    rlbl = 2
!
!  Begin first conjugate gradient loop.
!  Ask the user for a preconditioner solve.
!
  else if ( rlbl == 2 ) then

    iter = 1

    job = 2
    rlbl = 3
!
!  Compute the direction.
!  Ask the user to compute ALPHA.
!  Save A*P to Q.
!
  else if ( rlbl == 3 ) then

    rho = dot_product ( r, z )

    if ( 1 < iter ) then
      beta = rho / rho_old
      z(1:n) = z(1:n) + beta * p(1:n)
    end if

    p(1:n) = z(1:n)

    job = 1
    rlbl = 4
!
!  Compute current solution vector.
!  Ask the user to check the stopping criterion.
!
  else if ( rlbl == 4 ) then

    pdotq = dot_product ( p, q )
    alpha = rho / pdotq
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * q(1:n)

    job = 4
    rlbl = 5
!
!  Begin the next step.
!  Ask for a preconditioner solve.
!
  else if ( rlbl == 5 ) then

    rho_old = rho
    iter = iter + 1

    job = 2
    rlbl = 3

  end if

  return
end

  subroutine write_tecplot(solution_filename,element_filename,num_nodes,element_order,coord,node_u)
  use prec
  implicit none

  character ( len = * ) solution_filename
  character ( len = * ) element_filename
  integer, intent(in) :: num_nodes
  integer, intent(in) :: element_order
  real(dp),dimension(2,num_nodes), intent(in) :: coord
  real(dp), dimension(num_nodes), intent(in) :: node_u

  integer :: k
  integer :: num_elements
  integer :: iel
  integer :: inode
  integer :: offset
  real(dp) :: tmp1,tmp2
  integer, dimension(element_order) :: conectivity

  open(unit=6,file=solution_filename,status='replace')
  rewind 6

  open(unit=5,file=element_filename,status='old')
  rewind 5
  read(5,*) num_elements

  write(6,*) 'TITLE = "FEM2D_ELASTICA_CG_SOLUTION"'
  write(6,*) 'VARIABLES = "X", "Y", "UX", "UY"'
  write(6,*) 'ZONE DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL,N=', num_nodes,',E=',num_elements

  offset=1
  do inode=1,num_nodes
      tmp1=node_u(offset)
      tmp2=node_u(offset+1)
      offset=2*inode+1
      write(6,*) coord(1,inode), coord(2,inode), tmp1, tmp2 
  end do

  do iel=1,num_elements
      read(5,*) k, (conectivity(inode) , inode=1,element_order)
      write(6,*) (conectivity(inode) , inode=1,element_order) 
  end do 
  rewind 5
  close(5) 
  close(6)

  return
  end
