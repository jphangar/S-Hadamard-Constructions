# Haagerup Invariant

This directory contains Maple code for computing the Haagerup invariant of a Hadamard matrix.
`haagerup_invariant.mpl` contains functions for computing the Haagerup invariant and checking if two different invariants are equal.
The script `haagerup_inv_check.mpl` checks that the Haagerup invariant of all four algebraically obtained parametrizations is equal.
The invariant value is saved to `Habc-Haagerup-Inv.m`, and also displayed in Table 5.1 in the thesis.

This script also finds all parameter variations (of the considered set of 192 variations) that result in a parameter matrix with the same invariant value as $H(a,b,c)$.
The script prints the values of a,b,c (in this order) which lead to the same invariant value:

```m
> read "haagerup_inv_check.mpl";
--- PARAMETER VARIANTS WITH EQUAL HAAG. INV. ---- 
[a, b, c]
[1/a, 1/b, 1/c]
[a, b, a/c*b]
[1/a, 1/b, 1/a/b*c]
[a, 1/c, 1/b]
[1/a, c, b]
[1/a, c, 1/a/b*c]
[a, 1/c, a/c*b]
[a, 1/a/b*c, 1/b]
[1/a, a/c*b, b]
[1/a, a/c*b, 1/c]
[a, 1/a/b*c, c]
[b, a, c]
[1/b, 1/a, 1/c]
[b, a, a/c*b]
[1/b, 1/a, 1/a/b*c]
[b, 1/c, 1/a]
[1/b, c, a]
[1/b, c, 1/a/b*c]
[b, 1/c, a/c*b]
[b, 1/a/b*c, 1/a]
[1/b, a/c*b, a]
[1/b, a/c*b, 1/c]
[b, 1/a/b*c, c]
[1/c, a, 1/b]
[c, 1/a, b]
[c, 1/a, 1/a/b*c]
[1/c, a, a/c*b]
[1/c, b, 1/a]
[c, 1/b, a]
[c, 1/b, 1/a/b*c]
[1/c, b, a/c*b]
[c, a/c*b, a]
[1/c, 1/a/b*c, 1/a]
[c, a/c*b, b]
[1/c, 1/a/b*c, 1/b]
[1/a/b*c, a, 1/b]
[a/c*b, 1/a, b]
[a/c*b, 1/a, 1/c]
[1/a/b*c, a, c]
[1/a/b*c, b, 1/a]
[a/c*b, 1/b, a]
[a/c*b, 1/b, 1/c]
[1/a/b*c, b, c]
[a/c*b, c, a]
[1/a/b*c, 1/c, 1/a]
[a/c*b, c, b]
[1/a/b*c, 1/c, 1/b]
```
