module precision1
  implicit none
  integer, parameter :: dp = KIND(1.d0)
  integer, parameter :: sp = KIND(1.0)
  integer, parameter :: qp = selected_real_kind(32)
end module precision1
