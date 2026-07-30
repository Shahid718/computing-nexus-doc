<style>H1{color:Black;}</style>

<style>H2{color:DarkOrange;}</style>

<style>p{color:Black;}</style>

# Why Do Type-Bound Procedures in Fortran Require CLASS Instead of TYPE? (A Simple Explanation)

`Class` allows polymorphism (working with derived types and their children), while `type` does not. 
Type-bound procedures are designed to work with inheritance.

**Example**
```fortran
! Parent type
type animal
    character(len=20) :: name
contains
    procedure :: speak => animal_speak
end type

! Child type that inherits from animal
type, extends(animal) :: dog
    character(len=20) :: breed
contains
    procedure :: speak => dog_speak  ! Override the method
end type
```

## What happens with class vs type:

**Using class ✓ (CORRECT):**

```fortran
subroutine animal_speak(this)
    class(animal), intent(in) :: this   ! Can accept ANY animal type
    print*, "Generic animal sound"
end subroutine

! This works with both:
type(animal) :: a
type(dog) :: d
call a%speak() ! ✓ Works
call d%speak() ! ✓ Works (uses dog's version)
```

**Using type ❌ (WRONG):**
```fortran
subroutine animal_speak(this)
    type(animal), intent(in) :: this   ! ONLY accepts exact animal type
    print*, "Generic animal sound"
end subroutine

! This fails:
type(dog) :: d
call d%speak()  !❌ ERROR dog can't be passed as animal
```

## Think of it like:

*	`type` = "You must bring exactly a red apple" (no substitutions)

*	`class` = "Bring any fruit" (apple, orange, banana all work)

### Technical Reason:

When you call `object%method()`, the compiler needs to know:

*   Which actual subroutine to call (especially if child types override it)

*   The object might be a child type, not just the parent type

class tells the compiler: "This could be THIS type OR ANY type that extends it"

### Bottom Line:

Type-bound procedures use class because they were designed for object-oriented programming where:

* Types can inherit from other types

*	Methods can be overridden

*	A variable might hold different type extensions at different times

type would break all these OOP features!