! ===================================================
!  Test 2: Module and Main Program share a name
! ===================================================
! A module named "test" and a program named "test"
! live in the same file. gfortran rejects this;
! nvfortran compiles and runs it without complaint.
!
! Try compiling with:
!   gfortran test2_module_vs_program_name.f90 -o test2
!   nvfortran test2_module_vs_program_name.f90 -o test2

module test
contains
  ! procedures
end module test

program test
  print *, "Running simulation"
end program test
