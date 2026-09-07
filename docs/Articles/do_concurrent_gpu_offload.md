# Does `do concurrent` Support GPU Offload for Procedures in External Modules?

`do concurrent` can be offloaded to an NVIDIA GPU automatically with `nvfortran -stdpar=gpu`. But offload of a procedure called inside `do concurrent` only works if the caller and the callee are compiled in the **same translation unit** — even when the callee is a `pure` module procedure that should, in principle, be a straightforward call.

**Example**

```fortran
do concurrent (i = 1:n_rows, j = 1:n_cols)
    call initial_micro(r(i, j), c(i, j))
    p(i, j) = c(i, j)
end do
```

Whether this offloads correctly depends entirely on *where* `initial_micro` is defined relative to the caller — not on whether it's `pure`, not on its interface, and not on whether it's a module procedure.

## What Happens in Each Case

**Case 1 — Internal procedure, single file ✓ (WORKS)**

```fortran
program random_2d_example
    implicit none
    ...
    do concurrent (i=1:n_cols, j=1:n_rows)
        call initial_micro(r(i, j), c(i, j))
        p(i, j) = c(i, j)
    end do
contains
    pure subroutine initial_micro(r_val, c_val)
        real, intent(in)  :: r_val
        real, intent(out) :: c_val
        c_val = r_val
    end subroutine initial_micro
end program random_2d_example
```

```
nvfortran -stdpar=gpu -gpu=cc61 -gpu=mem:managed -Minfo=accel test.f90 -o test
```

Compiler output confirms both the loop and the callee get GPU code:

```
random_2d_example:
     15, Generating NVIDIA GPU code
         15, Loop parallelized across CUDA thread blocks, CUDA threads(32) collapse(2)
initial_micro:
     22, Generating implicit acc routine seq
         Generating NVIDIA GPU code
```

**Case 2 — Module procedure, same file ✓ (WORKS)**

```fortran
module micro_mod
    implicit none
    private
    public :: initial_micro
contains
    pure subroutine initial_micro(r_val, c_val)
        real, intent(in)  :: r_val
        real, intent(out) :: c_val
        c_val = r_val
    end subroutine initial_micro
end module micro_mod

program random_2d_example
    use micro_mod, only: initial_micro
    implicit none
    ...
    do concurrent (i = 1:n_rows, j = 1:n_cols)
        call initial_micro(r(i, j), c(i, j))
        p(i, j) = c(i, j)
    end do
end program random_2d_example
```

Same result — still one translation unit, so it still offloads.

**Case 3 — Module in its own file, compiled separately ❌ (FAILS)**

```fortran
! micro_mod.f90
module micro_mod
    implicit none
    public :: initial_micro
contains
    pure subroutine initial_micro(r_val, c_val)
        real, intent(in)  :: r_val
        real, intent(out) :: c_val
        c_val = r_val
    end subroutine initial_micro
end module micro_mod
```

```fortran
! main.f90
program random_2d_example
    use micro_mod, only: initial_micro
    implicit none
    ...
    do concurrent (i = 1:n_rows, j = 1:n_cols)
        call initial_micro(r(i, j), c(i, j))
        p(i, j) = c(i, j)
    end do
end program random_2d_example
```

```
nvfortran -stdpar=gpu -gpu=cc61 -gpu=mem:managed -c micro_mod.f90   ! OK
nvfortran -stdpar=gpu -gpu=cc61 -gpu=mem:managed -c main.f90        ! FAILS
```

`micro_mod.f90` compiles cleanly on its own. `main.f90` does not:

```
NVFORTRAN-S-1074-Procedure call in Do Concurrent is not supported yet (main.f90: 14)
```

This is the normal way most real Fortran projects are structured — one module per file — and it's exactly the structure that breaks offload.

## The Mitigation: `include` the Module Source

Pulling the module's source into the caller's file at compile time with Fortran's `include` statement collapses Case 3 back into Case 2:

```fortran
include 'micro_mod.f90'

program random_2d_example
    use micro_mod, only: initial_micro
    implicit none
    ...
    do concurrent (i = 1:n_rows, j = 1:n_cols)
        call initial_micro(r(i, j), c(i, j))
        p(i, j) = c(i, j)
    end do
    print *, "Program executed successfully. First element p(1,1):", p(1, 1)
end program random_2d_example
```

```
nvfortran -stdpar=gpu -gpu=cc61 -gpu=mem:managed main.f90 -o main.out
./main.out
 Program executed successfully. First element p(1,1):   0.9079230
```

One command compiles and links the whole program — no separate `-c` step, because there's now only one translation unit.

### Technical Reason

`include 'file.f90'` is text substitution at compile time: the contents of `micro_mod.f90` are spliced directly into `main.f90` before parsing. The compiler never actually sees two separate files — it sees one program with the module and the caller side by side, which is the same situation as Case 2. A real separate-compilation setup (`use` across two independently compiled `.f90` files) does not give the offload analysis enough visibility into the callee, which is why Case 3 fails while Case 4 succeeds despite using the identical module source.

### Bottom Line

| Case | Procedure location | Separate translation unit? | Offloads? |
|---|---|---|---|
| 1 | Internal (`contains`), single file | No | ✓ |
| 2 | Module, same file as program | No | ✓ |
| 3 | Module, separate file, `use`d normally | Yes | ❌ `NVFORTRAN-S-1074` |
| 4 | Module, separate file, pulled in via `include` | No (merged at compile time) | ✓ |

With current `nvfortran` versions, `do concurrent` GPU offload requires the called procedure's source to be visible in the same translation unit as the caller — either because it's genuinely in the same file, or because `include` merges it in. Whether this is a fundamental restriction of the offload model or a limitation that will be relaxed with better cross-unit/LTO analysis in future releases is an open question — see the [discussion on the Fortran Discourse](https://fortran-lang.discourse.group/t/do-concurrent-nvidia-gpu-offload-does-it-support-external-procedures/11087/3) for more.
