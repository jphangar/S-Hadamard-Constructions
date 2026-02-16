# Algebraic Block Parametrization

This directory contains Maple code for determining block-structured parametrizations of $H_{12}$ using the algebraic method described in Section 5.2.
`algebraic_block_solver.mpl` is the script for solving the parameter block for each 4x4 permutation matrix.
The output of the solver is the following:

```m
> read "algebraic_block_solver.mpl";
Warning, (in ip) `i` is implicitly declared local |algebraic_block_solver.mpl:23|
Warning, (in sip) `i` is implicitly declared local |algebraic_block_solver.mpl:24|
Permutation block: Matrix(4, 4, [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]]), Valid solution count: 1, Trivial: false
Permutation block: Matrix(4, 4, [[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,1,0,0],[0,0,1,0],[1,0,0,0],[0,0,0,1]]), Valid solution count: 1, Trivial: false
Permutation block: Matrix(4, 4, [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,1,0,0],[0,0,0,1],[1,0,0,0],[0,0,1,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,1,0],[1,0,0,0],[0,0,0,1],[0,1,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,1,0],[0,1,0,0],[0,0,0,1],[1,0,0,0]]), Valid solution count: 1, Trivial: false
Permutation block: Matrix(4, 4, [[0,0,0,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,0,1],[1,0,0,0],[0,0,1,0],[0,1,0,0]]), Valid solution count: 1, Trivial: false
> 
```

Nontrivial solutions found by the solver are saved to `algebraic_solutions.m` in a list called `Soln_Blks` containing pairs `(P, K)` where P is the permutation block and K is a normalized parameter block.
The simplified nontrivial solutions are given in `algebraic_constructions.mpl`; these are the 3-parameter forms that we use throughout the thesis.

`compare_with_numerics.mpl` is a script for checking that each permutation and parameter block pair we found from the numerically obtained parametrizations can be constructed by one of the four algebraically obtained permutation/parameter block pairs.
This script does not print anything if the check was successful:

```m
> read "compare_with_numerics.mpl";
> 
```
