# Catalog of Butson Type S-Hadamard Matrices Referenced

This directory contains Butson type S-Hadamard matrices from the [Aalto repositiory](https://wiki.aalto.fi/spaces/Butson/pages/120482304/Matrices+up+to+monomial+equivalence) which we use as examples throughout our work.
Note that this is not an exhaustive collection, we only include matrices that are referenced in the thesis.

Each file containing S-Hadamard BH matrices is named `SBH-n-q-mpl.txt` where the file contains all $BH(n,q)$ matrices that are S-Hadamard.
The file is given in a Maple-readable format containing a list called `Evals` containing the matrices given in integer form.
For example, the `SBH-6-3-mpl.txt` file contains the following:

```{txt}
Evals := [
[[0, 0, 0, 0, 0, 0], [0, 0, 1, 1, 2, 2], [0, 1, 0, 2, 1, 2], [0, 1, 2, 0, 2, 1], [0, 2, 1, 2, 0, 1], [0, 2, 2, 1, 1, 0]]
];
```

In Maple, the integer form of the $BH(6,3)$ matrix can be formed as follows:

```{mpl}
with(LinearAlgebra):

read "SBH-6-3-mpl.txt":
M := Matrix(Evals[1]); # Integer form of matrix
```
