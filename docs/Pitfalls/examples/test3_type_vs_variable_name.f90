! ===================================================
!  Test 3: Derived Type and Variable share a name
! ===================================================
! Inside "createcar", the dummy argument is named "car",
! the same name as the derived type "car" it is declared
! as. gfortran rejects this with a chain of errors;
! nvfortran compiles and runs it without complaint.
!
! Try compiling with:
!   gfortran test3_type_vs_variable_name.f90 -o test3
!   nvfortran test3_type_vs_variable_name.f90 -o test3

module classmodule
implicit none
type :: car
  character(10) :: color
  real          :: fuel
end type car

contains
subroutine createcar(car, color, fuel)
  type(car), intent(out) :: car
  character(*), intent(in) :: color
  real, intent(in) :: fuel
  car%color = color
  car%fuel = fuel
  print*, color, fuel
end subroutine createcar

end module classmodule

program test
  use classmodule
  implicit none

  type(car) :: mycar
  call createcar(mycar, "red", 50.0)

end program test
