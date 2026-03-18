# Block-Structured Representation of Parametrizations

This directory contains Maple code and data for forming a block-structured representation of each parametrized matrix as described in Section 5.1.
De Launey's block-structured $BH(12,3)$ matrix is given as a Maple matrix in `H12.m` called `H12`.

The isomorphisms that we use to transform the base matrix of each parametrized matrix into $H_{12}$ are given in `isomorphisms.txt` in a Maple-readable format.
Note that these isomorphisms are for the specific `pslq_valid_param_mats.m` file given in this directory which is a copy of `../PSLQ/pslq_valid_param_mats.m`.
If the PSLQ parametrization script is re-run, the new parametrized matrices generated are not likely to be compatible with the given isomorphisms as parameters may have moved around.
These isomorphisms were generated using the BH equivalence algorithm (see `../BH-equivalence/pkg` for more details of usage).
The isomorphism usage is as follows in Maple:

```sh
read "isomorphisms.txt": # Three objects, CT_mats, P_mats and Q_mats

for i from 1 to nops(P_mats) do
    evalb( i in CT_mats ); # True if all omega-evaluations of i-th parametrized matrix are equivalent to H12 conjugate transpose
    # if parametrized matrix index is in CT_mats, TAKE THE CONJUGATE TRANSPOSE of the parametrized matrix before applying isomorphism

    P[i] := Matrix(P_mats[i]); # Monomial matrix acting on rows 
    Q[i] := Matrix(Q_mats[i]); # Monomial matrix acting on columns
    # isomorphism is P[i].(i-th parametrized matrix or its conjugate transpose).Q[i]
od;
```

The Maple script `block_structure.mpl` applies the isomorphisms to each parametrized matrix (or its conjugate transpose), and saves the resulting parameter matrices (i.e., matrix `M` such that `M.~H12` is the parametrized matrix) in `parameter_matrices.m` in a list called `ParMats`.
This script also determines the permutation and parameter block of each transformed parameter matrix, and saves the resulting blocks in `perm_and_parameter_blocks.m` in a list of pairs of blocks called `Blks`.
Also, the full parametrized matrices after applying the isomorphisms are saved in `block_structured_matrices.m` in a list called `TransfMats`.

Utility functions for structuring the matrices are given in `utils.mpl`.

There is no output printed by this script if each matrix could be successfully transformed into a block-structured matrix:

```m
> read "block_structure.mpl";
> 
```
