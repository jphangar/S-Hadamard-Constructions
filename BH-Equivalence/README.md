# Butson Hadamard Matrix Equivalence

This directory contains code for checking the equivalence of $BH(12,3)$ matrices used in Section 4.1 of the thesis.
All code given is in Magma.
We use a package by Egan, Flannery, and Ó Catháin, please see `pkg` directory for more details.

`omega_evaluations.txt` contains all $\omega$-evaluations of each parametrized $BH(12,3)$ matrix found via the PSLQ method from an approximation returned by the optimizer.
The file is given in a Magma friendly format.
It contains a list called `Evals` of 12x12 $BH(12,3)$ matrices given in integer form.
Each consecutive block of 27 elements in this list are all $\omega$-evaluations of one specific parametrized matrix.

`IdentifyBustonMatrices.m` goes through this list and checks if all $\omega$-evaluations of one specific parametrized matrix are equivalent to one of the two $BH(12,3)$ matrices.
If all matrices are successfully identified, there is no output from this script:

```m
> load "IdentifyButsonMatrices.m";
Loading "IdentifyButsonMatrices.m"
Loading "pkg/ButsonFunctions.m"
Loading "BH-12-3-Catalog.m"
Loading "omega_evaluations.txt"
> 
```

The script `MatchAaltoWithDeLauney.m` can be used to identify the $BH(12,3)$ matrix constructed by de Launey with the $BH(12,3)$ matrices in the Aalto repository:

```m
> load "MatchAaltoWithDeLauney.m";
Loading "MatchAaltoWithDeLauney.m"
Loading "pkg/ButsonFunctions.m"
Loading "BH-12-3-Catalog.m"
Are DL and H1 isomorphic? true
> 
```
