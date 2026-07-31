<style>H1{color:Black;}</style>

<style>H2{color:DarkOrange;}</style>

<style>p{color:Black;}</style>

# Assuming Default Array Output Is Compiler-Independent

A common misconception in Fortran is that writing a two-dimensional array using the default list-directed output

```fortran
write(unit, *) array
```

will always produce a matrix-like layout in the output file. While this appears to be true with **GNU Fortran (gfortran)**, it is **not guaranteed by the Fortran standard**.

The Fortran language specifies the values that must be written but **does not define how list-directed output should be formatted**. As a result, the exact appearance of the output is **compiler dependent**.

For example:

* **gfortran** often formats a two-dimensional array as multiple rows, giving the appearance of a matrix.
* **Intel Fortran (ifx)** may instead write the same array as a continuous stream of values, resulting in an implementation-dependent layout that does not necessarily preserve the matrix structure.

## Example

Consider the following program.

```fortran
program list_directed_output
    implicit none

    integer, parameter :: nx = 4, ny = 5
    real :: A(nx, ny)

    call random_number(A)

    write(*,*) A
end program list_directed_output
```

Suppose the generated array is

```text
0.21 0.54 0.89 0.14
0.73 0.48 0.32 0.91
0.66 0.18 0.75 0.27
0.44 0.96 0.53 0.81
0.15 0.69 0.37 0.58
```

### Typical output with gfortran

```text
 0.210000 0.540000 0.890000 0.140000
 0.730000 0.480000 0.320000 0.910000
 0.660000 0.180000 0.750000 0.270000
 0.440000 0.960000 0.530000 0.810000
 0.150000 0.690000 0.370000 0.580000
```

This resembles a 5×4 matrix and is therefore easy to inspect or import into plotting software.

### Typical output with Intel ifx

```text
 0.21000000 0.73000002
 0.66000003 0.44000000
 0.15000001 0.54000002
 0.47999999 0.18000001
 0.95999998 0.68999999
 0.88999999 0.31999999
 0.75000000 0.52999997
 0.37000000 0.14000000
 0.91000003 0.27000001
 0.81000000 0.57999998
```

Notice that the values are still written in **column-major order**, but the compiler chooses where to insert line breaks. The resulting file no longer resembles the original 2D array, even though it contains exactly the same data.

This behavior is perfectly valid because **the Fortran standard does not prescribe the formatting of list-directed output**.

## Recommendation

Never rely on

```fortran
write(unit, *) array
```

when generating data files for visualization, post-processing, or portability.

Instead, explicitly write one row (or one column) at a time with explicit format appropriate for your application.

```fortran
integer :: i

do i = 1, nx
    write(unit,'(*(F10.6,1X))') A(i,:)
end do
```

Using explicit formatting makes the output deterministic, portable, and independent of the compiler.

## Understanding the Format Specifiers

Consider the following statement:

```fortran
do i = 1, nx
    write(unit,'(*(F10.6,1X))') A(i,:)
end do
```

The format string

```fortran
'(*(F10.6,1X))'
```

controls exactly how each array element is written.

* **`F10.6`** is a floating-point format descriptor.

  * **`F`** indicates fixed-point (decimal) notation.
  * **`10`** specifies the total field width of **10 characters**.
  * **`6`** specifies that **6 digits** are printed after the decimal point.

  For example,

  ```text
    3.141593
    0.250000
   12.500000
  ```

  Each value occupies exactly 10 character positions, making the columns neatly aligned.

* **`1X`** inserts **one blank space** after each number, improving readability.

* **`*()`** is an unlimited repeat specifier introduced in Fortran 2008. It tells the compiler to repeatedly apply the enclosed format (`F10.6,1X`) until every element in the array section `A(i,:)` has been written. This means the same format works for arrays of any length without needing to know the number of columns in advance.

### What Does `*` Mean in `write(unit, *)`?

The statement

```fortran
write(unit, *) A
```

uses **list-directed output**, where the asterisk (`*`) tells the compiler:

> "Choose an appropriate format automatically."

This is convenient for quick debugging because you do not have to specify a format descriptor. However, the Fortran standard **does not define the exact appearance of list-directed output**. Each compiler is free to choose how values are arranged, how many values appear on each line, and how much spacing is used.

As a result, code such as

```fortran
write(unit, *) A
```

may produce matrix-like output with one compiler (e.g., **gfortran**) and a completely different layout with another (e.g., **Intel ifx**), even though both outputs contain exactly the same numerical values.

```{tip}
For portable and reproducible output, always prefer an explicit format.
```