<style>H1{color:Black;}</style>

<style>H2{color:DarkOrange;}</style>

<style>p{color:Black;}</style>

# Fortran Naming Conflicts in Practice: Code Samples and Compiler Behavior

This article is a follow-up to *Understanding Naming Conflicts in Fortran: Program Names vs File Names*. That article covered the theory. Here, we run actual code through two different compilers — `gfortran` and `nvfortran` — and compare what happens.

If you only read one section, make it this: **the same code can compile fine on one compiler and fail on another.** This is the core portability issue every Fortran developer runs into sooner or later.

Each test below is shown inline as a code block, and is also available as a standalone, runnable file in the [`examples/`](./examples) folder, so you can compile and try these yourself without copy-pasting out of markdown.

---

## Setup

Each test below uses its own file, and both compilers are recent versions:

```text
$ gfortran --version
GNU Fortran (Ubuntu 15.2.0-16ubuntu1) 15.2.0

$ nvfortran --version
nvfortran 24.11-0 64-bit target on x86-64 Linux -tp haswell
```

Each test is compiled the same simple way, using its own file name:

```bash
gfortran <filename>.f90 -o <output_name>
nvfortran <filename>.f90 -o <output_name>
```

---

## Test 1: A Variable With the Same Name as the Program

```fortran
program test
  character :: test
  print*, "Test"
end program test
```

*Runnable file: [`examples/test1_variable_vs_program_name.f90`](./examples/test1_variable_vs_program_name.f90)*

Here, the program is called `test`, and inside it we also declare a variable called `test`. This is the kind of thing that looks harmless but touches a real rule in Fortran: a name can't be reused for two different purposes in the same scope.

**gfortran result:**

```text
test1_variable_vs_program_name.f90:13:19:
   13 |   character :: test
      |                   1
Error: Symbol 'test' at (1) cannot have a type
```

**nvfortran result:**

```text
NVFORTRAN-S-0043-Illegal attempt to redefine subprogram or entry name test (test1_variable_vs_program_name.f90: 13)
  0 inform,   0 warnings,   1 severes, 0 fatal for test
```

Both compilers reject this, though the error messages are worded differently. This is a genuine conflict, not a compiler quirk — the program name `test` and the variable name `test` collide inside the same program unit, and every standards-conforming compiler should catch it.

```{note}
When both compilers reject the same code, you can be confident it maybe a real language violation, not a portability gray area.
```

---

## Test 2: A Module and a Program With the Same Name

```fortran
module test
contains
  ! procedures
end module test

program test
  print *, "Running simulation"
end program test
```

*Runnable file: [`examples/test2_module_vs_program_name.f90`](./examples/test2_module_vs_program_name.f90)*

This is where things get interesting.

**gfortran result:**

```text
   12 | module test
      |        2~~~
......
   17 | program test
      |            1
Error: Global name ‘test’ at (1) is already being used as a MODULE at (2)
```

**nvfortran result:**

```text
$ ./test2
Running simulation
```

`nvfortran` compiles and runs this without complaint. `gfortran` refuses to build it at all.

This is a real divergence between compilers, and it's a good example of why "it compiles on my machine" isn't the same as "it's correct code." A module and a program share the same *global namespace* in Fortran, so reusing a name across them is risky even when one particular compiler happens to allow it.

```{warning}
Code that compiles under nvfortran but not gfortran (or vice versa) is not portable, even if it "works." Treat compiler acceptance as a hint, not a guarantee of correctness.
```

---

## Test 3: A Derived Type with the Same Name as the Variable

```fortran
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
```

*Runnable file: [`examples/test3_type_vs_variable_name.f90`](./examples/test3_type_vs_variable_name.f90)*

When Test 3 is compiled on its own with the `car` subroutine argument reusing the derived type name `car` as a variable name, gfortran raises a chain of errors:

```text
test3_type_vs_variable_name.f90:22:12:

   22 |   type(car), intent(out) :: car
      |            1
Error: Derived type ‘car’ at (1) is being used before it is defined
test3_type_vs_variable_name.f90:25:7:

   25 |   car%color = color
      |       1
Error: Symbol ‘car’ at (1) has no IMPLICIT type
test3_type_vs_variable_name.f90:26:7:

   26 |   car%fuel = fuel
      |       1
Error: Symbol ‘car’ at (1) has no IMPLICIT type
test3_type_vs_variable_name.f90:21:22:

   21 | subroutine createcar(car, color, fuel)
      |                      1~~
Error: Symbol ‘car’ at (1) has no IMPLICIT type
test3_type_vs_variable_name.f90:33:7:

   33 |   use classmodule
      |       1
Fatal Error: Cannot open module file ‘classmodule.mod’ for reading at (1): No such file or directory
compilation terminated.
```

**nvfortran result:**

```text
$ ./test3
red    50.00000
```

Again, nvfortran accepts what gfortran rejects. The lesson is the same as Test 2: naming a variable identically to its own derived type is legal in some compilers' interpretation and not in others.

---

## Side-by-Side Comparison

| Test | Code Pattern | gfortran | nvfortran |
|------|-------------|----------|-----------|
| 1 | Variable name same as program name | ❌ Rejected | ❌ Rejected |
| 2 | Module name same as program name | ❌ Rejected | ✅ Compiles and runs |
| 3 | Variable name same as its derived type | ❌ Rejected | ✅ Compiles and runs |

The pattern here is worth remembering: **conflicts within a single program unit (Test 1) are caught everywhere. Conflicts across separate units — module vs. program, type vs. variable — are where compilers start to disagree.**

---

## Why This Matters for Portability

If you're writing code for your own use on a single machine, a compiler that quietly accepts a naming clash might feel like a convenience. But scientific and HPC codebases are rarely built with just one compiler forever. Code gets:

* Ported to a new cluster with a different toolchain
* Picked up by a collaborator using Intel's `ifx`/`ifort`
* Rebuilt years later on a newer compiler version with stricter checks

A project that "just works" on `nvfortran` can fail outright the first time someone runs `gfortran` on it. That's not a bug in either compiler — both are allowed some latitude in how strictly they enforce certain naming rules, and the Fortran standard doesn't always spell out every case with total precision.

```{tip}
If you want code that survives a compiler switch, don't rely on what your current compiler *allows*. Rely on what good naming practice *avoids* in the first place.
```

---

## Practical Naming Guidelines (Recap)

Based on the tests above, a few concrete habits go a long way:

* Never give a variable the same name as the program or module it lives in
* Never give a module and a program the same name, even if your compiler allows it
* Avoid naming a variable the same as its own derived type
* When in doubt, test on more than one compiler before trusting that "it compiles" means "it's correct"

These are small habits, but they're the difference between code that only works for you and code that works for whoever inherits your project next.

---

## Final Takeaway

Compiler behavior is not a substitute for understanding the Fortran standard. Two compilers can look at the same file and reach different verdicts, and when that happens, the safer compiler's answer (usually the stricter one) is the one to trust. Testing across compilers early — even with something as small as the three cases above — can save a lot of debugging time later in a project's life.

For the conceptual background behind these tests — what actually distinguishes a program name, a file name, and an executable name — see the [companion article](naming_conflict.md).