! ===================================================
!  Test 1: Variable and Main Program share a name
! ===================================================
! A variable named "test" is declared inside a program
! also named "test". This is rejected by every
! standards-conforming compiler.
!
! Try compiling with:
!   gfortran test1_variable_vs_program_name.f90 -o test1
!   nvfortran test1_variable_vs_program_name.f90 -o test1

program test
  character :: test
  print*, "Test"
end program test
