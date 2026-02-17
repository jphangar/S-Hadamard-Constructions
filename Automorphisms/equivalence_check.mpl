kernelopts(printbytes=false):

with(LinearAlgebra):
with(NumberTheory):
with(combinat):

read "H12_automorphisms.txt":
read "matrix_constructions.mpl":
read "../Block-Structure/utils.mpl":


# Function for determining if two matrices are equivalent using automorphisms of H12
isEquivalent := proc(H1::Matrix, H2::Matrix)
global P_mats, Q_mats;
local i, P, Q;
    for i from 1 to nops(P_mats) do
        P, Q := Matrix(P_mats[i]), Matrix(Q_mats[i]);

        if Equal(modOmega(P.H1.Q), H2) then
            return true;
        fi;
    od;
    return false;
end:

for H1 in Mats do
    for H2 in Mats do
        if not isEquivalent(H1, H2) then
            printf("ERROR: Failed to find isomorphism from Aut(H12)\n");
        fi;
    od;
od:
