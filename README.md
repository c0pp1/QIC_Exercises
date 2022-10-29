# Quantum Information and Computing Exercises
This is a repository to store all the exercises done during the lessons and to have the source code easily accessible from other machines, as VMs on CloudVeneto.
The exercises are developed in `FORTAN90`

## Table of Contents
* [Exercise 01](#exercise-01)
* [Exercise 02](#exercise-02)

## Exercise 01
This simple exercise is focused on building the development environment as well as start playing with `FORTRAN90`. This GitHub repository is part of this exercise.
It is divided in three parts:  

1. **Setup**
    * Create a working directory.
    * Open a code editor - i.e. emacs, Vim, Visual Studio editor, Atom - and write your first program in `FORTRAN`.
    * Submit a test job.
    * (Optional) Connect to CloudVeneto cluster via ssh and repeat the execution.
2. **Number precision**  
    Integer and real numbers have a finite precision. Explore the limits of `INTEGER` and `REAL` in Fortran.
    * Sum $2.000.000$ and $1$ with `INTEGER*2` and `INTEGER*4`.
    * Sum $π \cdot 10^{32}$ and $\sqrt{2} \cdot 10^{21}$ with single and double precision.
3. **Performance testing**  
    Matrix-matrix multiplication is many times the bottleneck of linear algebra computations.
    * Write explicitly the matrix-matrix multiplication loop in two different orders.
    * Use the `FORTRAN` intrinsic function.
    * Increase the matrix size and track the code performance using `FORTRAN` basic date and time routines (i.e. `CPU TIME`).
    * Use different optimization flags available with the compiler and compare the performance with increasing matrix size.

## Exercise 02