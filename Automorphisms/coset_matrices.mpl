kernelopts(printbytes=false):

with(LinearAlgebra):
with(NumberTheory):
with(combinat):

read "H12_automorphisms.txt":
read "matrix_constructions.mpl":
read "../Block-Structure/utils.mpl":
read "Habc.m":
read "../Block-Structure/H12.m":


# Function for determining if matrix is in the list of matrices
check_member := proc(H::Matrix, mats::list)
local i;
    for i from 1 to nops(mats) do
        if Equal(H, mats[i]) then
            return true, i;
        fi;
    od;
    return false, 0;
end:


# Find 96 unique matrices equivalent to Habc using cosets, and track automorphisms belonging to each coset
eq_mats := []:
cosets_list := []:
for i from 1 to nops(P_mats) do
    P, Q := Matrix(P_mats[i]), Matrix(Q_mats[i]);
    new_mat := modOmega(P.Habc.Q);

    is_member, idx := check_member(new_mat, eq_mats);
    if not is_member then
        eq_mats := [op(eq_mats), new_mat];
        cosets_list := [op(cosets_list), [[P, Q]]];
    else
        cosets_list[idx] := [op(cosets_list[idx]), [P, Q]];
    fi;
od:
printf("Number of Cosets: %d\n", nops(cosets_list)):


# Divide matrix corresponding to cosets into quarters based on permutation block
perm_mats := [op(P1), op(P2)]:
parameter_blocks_by_perm := [[], [], [], []]:
for mat in eq_mats do
    M := getParameterMatrix(mat, 12, w);
    K := M[1..4,1..4];
    for i from 1 to 4 do
        genmat := generateParameterMatrix(perm_mats[i], K);
        if Equal(mat, H12*~genmat) then
            parameter_blocks_by_perm[i] := [op(parameter_blocks_by_perm[i]), normalizeMat(K, 4)];
            break;
        fi;
    od:
od:


# The first permutation block is Habc's permutation block, print the construction of the coset matrix using Habc
for parameter_blk in parameter_blocks_by_perm[1] do
    eval_blk := eval(N1, {a=parameter_blk[3,3], b=parameter_blk[2,4], c=parameter_blk[2,3]});
    if Equal(eval_blk, parameter_blk) then
        printf("Construction form: H(%a, %a, %a)\n", parameter_blk[3,3], parameter_blk[2,4], parameter_blk[2,3]);
    fi;
od:
