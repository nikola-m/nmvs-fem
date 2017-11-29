module output

!*****************************************************************************80
!
! Module with routines for writing postprocessing files in Tecplot and Paraview.
!
!*****************************************************************************80
!
use prec

implicit none

contains

subroutine write_tecplot(solution_filename,element_filename,num_nodes,element_order,coord,node_u)

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

end module