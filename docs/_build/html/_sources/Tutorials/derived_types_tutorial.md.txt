---
orphan: false
myst:
  html_meta:
    "description": "Fortran derived data types tutorial from intrinsic types to type-bound procedures"
    "keywords": "Fortran, derived types, tutorial, OOP, grid example"
---

# Fortran Derived Types Tutorial

This tutorial introduces Fortran derived data types, from basic intrinsic types to object-oriented type-bound procedures.

```{admonition} Prerequisites
:class: note
- Basic Fortran knowledge (variables, subroutines)
- Fortran 90+ compiler (gfortran, ifort, nvfortran)
```

```{admonition} Learning objective
:class: tip
By the end of this tutorial, you will:
- Define and use derived types
- Apply automatic constructors
- Pass derived types to subroutines
- Implement type-bound procedures
```

## Complete Code Example

The complete tutorial code is available in `derived_types_tutorial.f90`

## Concept Breakdown

**1. Intrinsic Data Types**

Purpose: Basic variable declaration without structure.

```fortran
program grid_example
implicit none

    integer :: nx, ny, dx, dy
  
    nx = 50
    ny = 50
    dx = 2
    dy = 2

end program grid_example
```

Characteristics:

* Simple, but each variable is independent

* No logical grouping

* Limited scalability

**2. Derived Data Types**

Purpose: Group related variables into a single logical unit.

```fortran
program grid_example
implicit none
    
    type gridType
        integer :: nx, ny, dx, dy
    end type gridType

    type(gridType) :: grid

    grid%nx = 50
    grid%ny = 50
    grid%dx = 2
    grid%dy = 2
    
end program grid_example
```

Characteristics:

*   Groups multiple variables into a single composite object

*   Provides logical grouping of related data

*   Scalable and reusable across programs and subroutines

*   Access components using the `%` operator (e.g., `grid%nx`)


**3. Derived Type with Constructor**

Purpose: Group related data into a single unit.

```fortran
program test
implicit none

    type gridType
        integer :: nx, ny, dx, dy
    end type gridType

    type(gridType) :: grid

    ! Automatic structure constructor
    grid = gridType (nx=50, ny=50, dx=2, dy=2)

end program
```

Characteristics:

* `gridType` groups four integers logically

* Automatic constructor `gridType(...)` is implicitly available

*  Named arguments improve readability

**4. Subroutine with Derived Type Argument**

Purpose: Initialize derived types through procedural programming.

```fortran
program test
implicit none

    type gridType
        integer :: nx, ny, dx, dy
    end type gridType

    type(gridType) :: grid

    call createGrid (grid, nx=50, ny=50, dx=2, dy=2)

contains

    subroutine createGrid (this, nx, ny, dx, dy)
        type(gridType), intent(inout) :: this
        integer, intent(in) :: nx, ny, dx, dy

        this%nx = nx
        this%ny = ny
        this%dx = dx
        this%dy = dy

        print*, this%nx, this%ny, this%dx, this%dy

    end subroutine createGrid

end program
```

Key Points:

* Use `intent(inout)` rather than `intent(out)` (preserves existing data)

* Component access via `%` operator: `this%nx`

* Subroutine contained in program (implicit interface)

```{admonition} Tip
:class: tip
`intent(inout)` is recommended over `intent(out)` because it avoids undefined initial values if the derived type has `allocatable` components or later extensions.
```

**5. Type-Bound Procedures (Object-Oriented Fortran)**

Purpose: Attach methods directly to types (Fortran 2003+).

```fortran
module gridType_module
implicit none

    type gridType
        integer :: nx, ny, dx, dy
    contains
        procedure :: createGrid
    end type gridType

contains

    subroutine createGrid (this, nx, ny, dx, dy)
        class(gridType), intent(inout) :: this   ! Must be CLASS
        integer, intent(in) :: nx, ny, dx, dy

        this%nx = nx
        this%ny = ny
        this%dx = dx
        this%dy = dy

        print*, this%nx, this%ny, this%dx, this%dy

    end subroutine createGrid

end module gridType_module

program test
use gridType_module
implicit none

    type(gridType) :: grid

    call grid%createGrid (nx=50, ny=50, dx=2, dy=2)

end program
```

**Why `CLASS` instead of `TYPE`?**

|Type	| Use  Case|	Polymorphic|
|-----|----------|------------|
|`TYPE`	| Non-polymorphic arguments|	No|
|`CLASS`|	Polymorphic arguments (type extensions)	|Yes|

For type-bound procedures, the passed object is polymorphic, so `CLASS` is mandatory.

**Comparison Table**

|Concept	|Syntax Example|
|---------|--------------|
|Intrinsic type|	`integer :: nx`|
|Derived type definition	| `type gridType ... end type`
|Automatic constructor|	`gridType(nx=50, ny=50, dx=2, dy=2)`
|Component access	| `grid%nx`|
|Subroutine with derived type	|`call createGrid(grid, ...)`|
|Type-bound procedure|	`call grid%createGrid(...)`|
|Polymorphic argument|	`class(gridType), intent(inout) :: this`

**Compilation and Execution**

```bash
#compile
gfortran -o derived_types_tutorial derived_types_tutorial.f90

# Run
./derived_types_tutorial
```

Expected Output (Example 4):

```
50   50   2   2
```
**Key Takeaways**

```{admonition} Important
:class: important
1. **Use `intent(inout)`** for derived type arguments in subroutines
2. **Prefer `CLASS` over `TYPE`** for type-bound procedure arguments
3. **Encapsulate types in modules** for reusability and clean interfaces
4. **Use automatic constructors** with named arguments for clarity
```

**Exercises**

* Modify the grid to include a real component called `resolution`

* Add a method `printGrid()` to display formatted output

* Create an extended type `grid3DType` with an `nz` component

*  Implement validation in `createGrid` to ensure positive grid dimensions


---

## Further Reading

- [Modern Fortran Explained](https://www.elsevier.com/books/modern-fortran-explained/metcalf/978-0-08-100063-6) (Metcalf, Reid, Cohen)
- [Introduction to Programming with Fortran](https://link.springer.com/book/10.1007/978-3-319-75502-1) (Chivers , Sleightholme)
- [Intel Developer Guide - Derived Types](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2025-2/derived-data-types.html)