# Equivalences via Automorphisms

The automorphism group of $H_{12}$ was computed using the BH equivalence algorithm package, `../BH-equivalence/pkg` for more details.
The automorphisms are given in a Maple-readable format in `H12_automorphisms.txt`.
The file usage is as follows in Maple:

```sh
read "H12_automorphisms.txt"; # Two lists P_mats and Q_mats
read "../Block_Structure/H12.m"; # Matrix H12
read "../Block_Structure/utils.mpl";

for i from 1 to nops(P_mats) do
    P := Matrix(P_mats[i]);
    Q := Matrix(Q_mats[i]);
    # (P,Q) is an automorphism of H12
    evalb(Equal(modOmega(P.H12.Q), H12)); # True, modOmega simplifies w powers
od;
```

The four parametrized matrices (in simplified form, over common parameters) are given in `matrix_constructions.mpl`.
To check that all four matrices are equivalent, run `equivalence_check.mpl` which finds an isomorphism between each pair of these matrices given by an automorphism of $H_{12}$.
The script prints no output when all matrices are equivalent:

```m
> read "equivalence_check.mpl";
>   
```

The representative matrix we choose for our family of matrices is denoted $H(a,b,c)$ which is saved in `Habc.m` called `Habc`.
The subgroup of $Aut(H_{12})$ formed by automorphisms of $H(a,b,c)$ is saved to `automorphism-subgroup-mats.txt` in a Magma-readable format.
The file contains a list of pairs of matrices called `Automorphisms`.
The Magma script `automorphismSubgrp.m` checks that these automorphisms indeed form a subgroup of $Aut(H_{12})$ and finds generators for the subgroup:

```m
> load "automorphismSubgrp.m";
Loading "automorphismSubgrp.m"
Loading "../BH-Equivalence/pkg/ButsonFunctions.m"
Loading "../BH-Equivalence/BH-12-3-Catalog.m"
Loading "automorphism-subgroup-mats.txt"
Order of subgroup: 9
Is subgroup normal? false
Generator 1,
 P: [w 0 0 0 0 0 0 0 0 0 0 0]
[0 w 0 0 0 0 0 0 0 0 0 0]
[0 0 w 0 0 0 0 0 0 0 0 0]
[0 0 0 w 0 0 0 0 0 0 0 0]
[0 0 0 0 w 0 0 0 0 0 0 0]
[0 0 0 0 0 w 0 0 0 0 0 0]
[0 0 0 0 0 0 w 0 0 0 0 0]
[0 0 0 0 0 0 0 w 0 0 0 0]
[0 0 0 0 0 0 0 0 w 0 0 0]
[0 0 0 0 0 0 0 0 0 w 0 0]
[0 0 0 0 0 0 0 0 0 0 w 0]
[0 0 0 0 0 0 0 0 0 0 0 w],
 Q: [w^2 0 0 0 0 0 0 0 0 0 0 0]
[0 w^2 0 0 0 0 0 0 0 0 0 0]
[0 0 w^2 0 0 0 0 0 0 0 0 0]
[0 0 0 w^2 0 0 0 0 0 0 0 0]
[0 0 0 0 w^2 0 0 0 0 0 0 0]
[0 0 0 0 0 w^2 0 0 0 0 0 0]
[0 0 0 0 0 0 w^2 0 0 0 0 0]
[0 0 0 0 0 0 0 w^2 0 0 0 0]
[0 0 0 0 0 0 0 0 w^2 0 0 0]
[0 0 0 0 0 0 0 0 0 w^2 0 0]
[0 0 0 0 0 0 0 0 0 0 w^2 0]
[0 0 0 0 0 0 0 0 0 0 0 w^2]
Generator 2,
 P: [0 0 0 0 w^2 0 0 0 0 0 0 0]
[0 0 0 0 0 0 w^2 0 0 0 0 0]
[0 0 0 0 0 0 0 w^2 0 0 0 0]
[0 0 0 0 0 w^2 0 0 0 0 0 0]
[0 0 0 0 0 0 0 0 Id(G) 0 0 0]
[0 0 0 0 0 0 0 0 0 0 Id(G) 0]
[0 0 0 0 0 0 0 0 0 0 0 Id(G)]
[0 0 0 0 0 0 0 0 0 Id(G) 0 0]
[w 0 0 0 0 0 0 0 0 0 0 0]
[0 0 w 0 0 0 0 0 0 0 0 0]
[0 0 0 w 0 0 0 0 0 0 0 0]
[0 w 0 0 0 0 0 0 0 0 0 0],
 Q: [0 0 0 0 w 0 0 0 0 0 0 0]
[0 0 0 0 0 0 0 w 0 0 0 0]
[0 0 0 0 0 w 0 0 0 0 0 0]
[0 0 0 0 0 0 w 0 0 0 0 0]
[0 0 0 0 0 0 0 0 Id(G) 0 0 0]
[0 0 0 0 0 0 0 0 0 0 0 Id(G)]
[0 0 0 0 0 0 0 0 0 Id(G) 0 0]
[0 0 0 0 0 0 0 0 0 0 Id(G) 0]
[w^2 0 0 0 0 0 0 0 0 0 0 0]
[0 0 0 w^2 0 0 0 0 0 0 0 0]
[0 w^2 0 0 0 0 0 0 0 0 0 0]
[0 0 w^2 0 0 0 0 0 0 0 0 0]
>
```

The script `coset_matrics.mpl` finds the 96 distinct matrices equivalent to $H(a,b,c)$ corresponding to each coset of the subgroup of automorphisms of $H(a,b,c)$ in $Aut(H_{12})$.
This script also checks that each of these matrices is a block-structured parametrization of $H_{12}$.
The script prints the following output, which is shown in Table 5.2 in the thesis, that shows that the matrices with the same permutation block as $H(a,b,c)$ can be constructed from $H(a,b,c)$:

```m
> read "coset_matrices.mpl";
Number of Cosets: 96
Construction form: H(a, b, c)
Construction form: H(1/a, c, b)
Construction form: H(1/c, 1/a/b*c, 1/a)
Construction form: H(c, 1/a, 1/a/b*c)
Construction form: H(b, a, b/c*a)
Construction form: H(1/a, 1/b, 1/a/b*c)
Construction form: H(1/c, b, b/c*a)
Construction form: H(1/b, b/c*a, a)
Construction form: H(1/a/b*c, 1/c, 1/b)
Construction form: H(b/c*a, c, a)
Construction form: H(1/a/b*c, b, 1/a)
Construction form: H(c, b/c*a, b)
Construction form: H(1/c, a, 1/b)
Construction form: H(b, 1/a/b*c, c)
Construction form: H(1/a, b/c*a, 1/c)
Construction form: H(a, 1/a/b*c, 1/b)
Construction form: H(b, 1/c, 1/a)
Construction form: H(b/c*a, 1/b, 1/c)
Construction form: H(c, 1/b, a)
Construction form: H(1/a/b*c, a, c)
Construction form: H(b/c*a, 1/a, b)
Construction form: H(1/b, c, 1/a/b*c)
Construction form: H(1/b, 1/a, 1/c)
Construction form: H(a, 1/c, b/c*a)
```
