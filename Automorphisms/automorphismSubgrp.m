load "../BH-Equivalence/pkg/ButsonFunctions.m";
load "../BH-Equivalence/BH-12-3-Catalog.m";

G<w> := SmallGroup(3,1);
RG := GroupAlgebra(Rationals(), G);
MAlg12 := MatrixRing(RG, 12);

load "automorphism-subgroup-mats.txt";

// de Launey's BH(12,3) matrix
H1 := IntMatrixtoGH(BH12_3_1989[1], G);

// Generate automorphism group of H1
AG := GHAutomorphismGroup(H1,G);

subgrp_P := [];
subgrp_Q := [];
subgrp_perm := [];

// Build list of P and Q matrices in subgroup
for mats in Automorphisms do
    subgrp_P := Append(subgrp_P, mats[1]);
    subgrp_Q := Append(subgrp_Q, mats[2]);
end for;


inverse := function(x);
    if x eq 0 then
        return 0;
    else
        return x^-1;
    end if;
end function;


// Identify the permutations that correspond to the subgroup
for perm in AG do
    Mats := GHPermAutToAut(perm, G);
    P := Mats[1]; Q := Mats[2]^-1;
    InvQ := Matrix([[inverse(Q[v, u]) : u in [1..Ncols(Q)]] : v in [1..Nrows(Q)]]);
    if P*H1*InvQ ne H1 then
        printf "Found permutation in automorphism group that is not an automorphism: %o", perm;
    end if;
    // Check if this is in automorphism subgroup
    idx_P := Index(subgrp_P, P);
    idx_Q := Index(subgrp_Q, InvQ);
    if idx_P ne 0 and idx_Q ne 0 then
        subgrp_perm := Append(subgrp_perm, perm);
    end if;
end for;


// Form a subgroup from the matching permutations
H := sub<AG | subgrp_perm>;
printf "Order of subgroup: %o\n", Order(H);
printf "Is subgroup normal? %o\n", IsNormal(AG, H);


// Find generators for subgroup (just one choice is fine)
generators_H := {};
elts := Generators(H); // returns the 8 non-trivial elements
S := SetToSequence(Subsets(elts, 2)); // try all pairs of generators
for gen_pair in S do
    if sub<AG | gen_pair> eq H then
        generators_H := gen_pair;
        break;
    end if;
end for;

// Convert the generators from permutations to matrices
counter := 1;
for gen_H in generators_H do
    Mats := GHPermAutToAut(gen_H, G);
    P := Mats[1]; Q := Mats[2]^-1;
    InvQ := Matrix([[inverse(Q[v, u]) : u in [1..Ncols(Q)]] : v in [1..Nrows(Q)]]);
    printf "Generator %o,\n P: %o,\n Q: %o\n", counter, P, InvQ;
    counter := counter+1;
end for;
