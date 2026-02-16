// Determining which of the two BH(12,3) matrices in the Aalto table are equivalent to de Launey's BH(12,3)

load "pkg/ButsonFunctions.m";

// Contains BH(12,3) from de Launey 1989 Example 4.6 and BH(12,3) matrices from Aalto table
load "BH-12-3-Catalog.m";

// Order 3 group
G := SmallGroup(3,1);

// de Launey's BH(12,3)
DL := IntMatrixtoGH(BH12_3_1989[1], G);

// H1 is first entry in Aalto table
H1 := IntMatrixtoGH(Alto_BH12_3[1], G);

// H2 is second entry in Aalto table
H2 := IntMatrixtoGH(Alto_BH12_3[2], G);

// Check equivalence with H1
b, sigma := GHEquivalence(DL, H1, G, G);
printf "Are DL and H1 isomorphic? %o\n", b;

// Check equivalence with H2 if H1 was not equivalent
if not b then
    b, sigma := GHEquivalence(DL, H2, G, G);
    printf "Are DL and H2 isomorphic? %o\n", b; 
end if;
