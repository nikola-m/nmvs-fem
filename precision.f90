module prec
!
! define precision for floating-point real
!
   ! single precision
   integer, parameter :: sp = selected_real_kind(6,37)
   ! double precision    
   integer, parameter :: dp = kind(1.0d0) !selected_real_kind(15,307) 
   ! quadruple precision
   integer, parameter :: qp = selected_real_kind(33,4931)

end module prec
