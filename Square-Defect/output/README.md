# Square Defect Computation Output

The provided script (`compute_defect.mpl`) generates the output files found in this directory.
Files named `sbh-n-q-defects.m` contain the square defect values for the S-Hadamard $BH(n,q)$ matrices, and `sbh-n-q-parametrizations.m` contain a subset of these matrices for which the full SLA-parametrization was S-Hadamard.
Note that all files contain Maple objects, not plain text, hence the file must be read in Maple for the contents to be visible.

The contents of the files can be read as follows in Maple:

```sh
read "sbh-12-6-defects.m":

for i from 1 to nops(sq_defects) do
    i; # Index of matrix
    sq_defects[i][1]; # Square defect value
    sq_defects[i][2]; # (n x n) parameter matrix defining the kernel of the S-Had. lin. approx.
    sq_defect[i][3]; # names of parameters
od;

read "sbh-12-6-parametrizations.m":

for i from 1 to nops(full_slas) do
    full_slas[i][1]; # Index of matrix
    full_slas[i][2]; # Integer form of BH(n,q) = H
    full_slas[i][3]; # Parameter matrix = R
    full_slas[i][4]; # Names of parameters
    # Parametrized matrix: EXP(I*2*Pi/q*H + I*R)
od;
```
