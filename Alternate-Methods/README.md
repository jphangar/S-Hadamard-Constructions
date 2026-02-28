# Alternate Approaches for Parametrizing Matrices

This directory contains code for alternate methods for parametrizing we attempted described in Section 5.4 of the thesis.

## Circulant S-Hadamard Matrices

`circulant_matrices.mpl` is the script used to generate the results given in Section 5.4.1 for checking if we can introduce one parameter to a circulant version of $H_{12}^*$ given by de Launey (1992).
The output of the script shows that only trivial parametrizations (either adding no parameters or multiplying both generating vectors by the parameter) produce S-Hadamard parametrizations of the matrix:

```m
> read "circulant_matrices.mpl";
Binary Parameter Vectors ...
Found parameterization: [1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]
Found parameterization: [a, a, a, a, a, a], [a, a, a, a, a, a]
FAIL COUNT: 4094
SUCCESS COUNT: 2
```

## Solving for Parameter Blocks

`solving_parameter_blks.mpl` is the script used to check if the $BH(12,3)$ matrices given in the Aalto repository have a block-structured parametrization in Section 5.4.2.
The script prints permutation blocks that lead to a valid parametrization and indicates if the parametrization is trivial.
All parametrizations found from this method were trivial:

```m
> read "solving_parameter_blks.mpl";
Warning, (in ip) `i` is implicitly declared local |solving_parameter_blks.mpl:35|
Warning, (in sip) `i` is implicitly declared local |solving_parameter_blks.mpl:36|
Parameter block solving for H1 ...
Permutation block: Matrix(4, 4, [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[1,0,0,0],[0,0,1,0],[0,0,0,1],[0,1,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,1,0,0],[0,0,0,1],[0,0,1,0],[1,0,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,0,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]]), Valid solution count: 1, Trivial: true
Parameter block solving for H2 ...
Permutation block: Matrix(4, 4, [[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,1,0,0],[0,0,0,1],[1,0,0,0],[0,0,1,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,1,0,0],[0,0,0,1],[0,0,1,0],[1,0,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,1,0],[0,1,0,0],[0,0,0,1],[1,0,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]), Valid solution count: 1, Trivial: true
Permutation block: Matrix(4, 4, [[0,0,1,0],[0,0,0,1],[0,1,0,0],[1,0,0,0]]), Valid solution count: 1, Trivial: true
> 
```

## SLA-Parametrization

`sla_parametrization.mpl` is the script used in Section 5.4.3 to find S-Hadamard SLA-parametrizations of the $BH(12,3)$ matrices given in the Aalto repository.
We use the basis of the kernel given in `sbh-12-3-defects.m` for both matrices as this basis will give the results shown in the thesis.
Using a different basis will lead to different results.
We also compute S-Hadamard SLA-parametrizations for normalized $H_{12}$; note that the basis for this is generated at runtime hence results will differ.
The output of the script is the following:

```m
> read "sla_parametrization.mpl";
Finding S-Hadamard SLA-parametrizations for Aalto BH(12,3) matrices ...
Idx 1, variable r[9,8] is an S-Hadamard SLA-parameterization
Idx 1, variable r[10,9] is an S-Hadamard SLA-parameterization
Idx 1, variable r[10,10] is an S-Hadamard SLA-parameterization
Idx 1, variable r[12,8] is an S-Hadamard SLA-parameterization
Idx 1, variable {r[9,8], r[10,9]} is an S-Hadamard SLA-parameterization
Idx 1, variable {r[9,8], r[12,8]} is an S-Hadamard SLA-parameterization
Idx 1, variable {r[10,9], r[12,8]} is an S-Hadamard SLA-parameterization
Idx 1, variable {r[9,8], r[10,9], r[12,8]} is an S-Hadamard SLA-parameterization
Idx 2, variable r[8,11] is an S-Hadamard SLA-parameterization
Idx 2, variable r[10,9] is an S-Hadamard SLA-parameterization
Idx 2, variable r[11,7] is an S-Hadamard SLA-parameterization
Idx 2, variable r[12,6] is an S-Hadamard SLA-parameterization
Finding S-Hadamard SLA-parametrizations for de Launey BH(12,3) matrix ...
Idx 1, variable r[10,9] is an S-Hadamard SLA-parameterization
Idx 1, variable r[11,9] is an S-Hadamard SLA-parameterization
Idx 1, variable r[11,10] is an S-Hadamard SLA-parameterization
```

Index 1 in the Aalto matrices is $H_1$ which is equivalent to $H_{12}$, the chosen basis gives a 3-parameter matrix.
Index 2 is $H_2$ which is equivalent to $H_{12}^*$, the chosen basis gives only 1-parameter matrices.
The last matrix shown in normalized $H_{12}$ which (in this run) gave only 1-parameter matrices.

There is a function in this script called `getNewBasisAndParametrize()` which we can repeatedly call to try to generate a new basis for the kernel and find all S-Hadamard SLA-parametrizations from this new basis.
By calling this function repeatedly, we observe a run where $H_1$ did not produce a 3-parameter matrix (where index 1 is $H_1$, index 2 is $H_2$, and index 3 is normalized $H_{12}$):

```m
> getNewBasisAndParametrize():
Idx 1, variable r[9,8] is an S-Hadamard SLA-parameterization
Idx 1, variable r[10,10] is an S-Hadamard SLA-parameterization
Idx 1, variable r[12,8] is an S-Hadamard SLA-parameterization
Idx 1, variable {r[9,8], r[12,8]} is an S-Hadamard SLA-parameterization
Idx 2, variable r[7,11] is an S-Hadamard SLA-parameterization
Idx 2, variable r[11,7] is an S-Hadamard SLA-parameterization
Idx 2, variable r[11,9] is an S-Hadamard SLA-parameterization
Idx 2, variable r[11,11] is an S-Hadamard SLA-parameterization
Idx 2, variable {r[11,9], r[11,11]} is an S-Hadamard SLA-parameterization
Idx 3, variable r[10,9] is an S-Hadamard SLA-parameterization
Idx 3, variable r[11,9] is an S-Hadamard SLA-parameterization
Idx 3, variable r[11,10] is an S-Hadamard SLA-parameterization
> 
```

We also observed one run where normalized $H_{12}$ did find a 3-parameter matrix:

```m
> getNewBasisAndParametrize():
Idx 1, variable r[9,8] is an S-Hadamard SLA-parameterization
Idx 1, variable r[10,9] is an S-Hadamard SLA-parameterization
Idx 1, variable r[10,10] is an S-Hadamard SLA-parameterization
Idx 1, variable r[12,8] is an S-Hadamard SLA-parameterization
Idx 1, variable {r[9,8], r[10,9]} is an S-Hadamard SLA-parameterization
Idx 1, variable {r[9,8], r[12,8]} is an S-Hadamard SLA-parameterization
Idx 1, variable {r[10,9], r[12,8]} is an S-Hadamard SLA-parameterization
Idx 1, variable {r[9,8], r[10,9], r[12,8]} is an S-Hadamard SLA-parameterization
Idx 2, variable r[8,11] is an S-Hadamard SLA-parameterization
Idx 2, variable r[10,9] is an S-Hadamard SLA-parameterization
Idx 2, variable r[11,7] is an S-Hadamard SLA-parameterization
Idx 2, variable r[12,6] is an S-Hadamard SLA-parameterization
Idx 3, variable r[8,11] is an S-Hadamard SLA-parameterization
Idx 3, variable r[10,12] is an S-Hadamard SLA-parameterization
Idx 3, variable r[11,9] is an S-Hadamard SLA-parameterization
Idx 3, variable {r[8,11], r[10,12]} is an S-Hadamard SLA-parameterization
Idx 3, variable {r[8,11], r[11,9]} is an S-Hadamard SLA-parameterization
Idx 3, variable {r[10,12], r[11,9]} is an S-Hadamard SLA-parameterization
Idx 3, variable {r[8,11], r[10,12], r[11,9]} is an S-Hadamard SLA-parameterization
> 
```

## Row-Wise Parametrization Using Backtracking

For the backtracking method described in Section 5.4.4 of the thesis, see directory `Backtracking-Method` for details.
