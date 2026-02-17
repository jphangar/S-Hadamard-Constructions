kernelopts(printbytes=false):

with(LinearAlgebra):
with(combinat):

read "../Square-Defect/defect_functions.mpl";


# Compute the square defects for each matrix, return a basis for the kernel
getSqDefectsForMats := proc(evals::list, q::integer)
local i, out, dfct, Rmat, freeVr;
    out := [];
    for i from 1 to nops(evals) do
        dfct, Rmat, freeVr := ComputeDefectSq( IntMatToBH(evals[i], q) );
        out := [op(out), [dfct, Rmat, freeVr]];
    od;
    out;
end:


# Return all S-Hadamard SLA-parametrizations found from given basis of kernel
findSlaParametrizations := proc(evals::list, defects::list, n::integer, q::integer)
local i, j, mats, H, p, var, var_dirs, M;
    mats := [];
    for i from 1 to nops(evals) do
    M := defects[i][2];
    H := BuildSymbolicMatrix(evals[i], M);

    if CheckSHadConditions(H, n, q) then
        printf("Idx %d is an S-Hadamard full SLA-parameterization\n", i);
        var_dirs := {op(var_dirs), p};
        mats := [op(mats), [i, H, M]];
        next;
    fi;

    # Find all one-parameter S-Hadamard matrices
    var_dirs := {};
    for p in defects[i][3] do
        M := eval(defects[i][2], { seq(var=0, var=(defects[i][3] minus {p})) });
        H := BuildSymbolicMatrix(evals[i], M);
        if CheckSHadConditions(H, n, q) then
            printf("Idx %d, variable %a is an S-Hadamard SLA-parameterization\n", i, p);
            var_dirs := {op(var_dirs), p};
            mats := [op(mats), [i, H, M]];
        fi;
    od;

    # Check combinations of one-parameter S-Hadamard matrices
    for p in choose(var_dirs) do
        if nops(p) < 2 then next fi;
        M := eval(defects[i][2], { seq(var=0, var=(defects[i][3] minus p)) });
        H := BuildSymbolicMatrix(evals[i], M); 
        if CheckSHadConditions(H, n, q) then
            printf("Idx %d, variable %a is an S-Hadamard SLA-parameterization\n", i, p);
            mats := [op(mats), [i, H, M]];
        fi;
    od;
    od;

    mats;
end:


printf("Finding S-Hadamard SLA-parametrizations for Aalto BH(12,3) matrices ...\n");
read "../SBH-Matrices/SBH-12-3-mpl.txt":
read "sbh-12-3-defects.m":
MatsAalto := findSlaParametrizations([Matrix(Evals[1]), Matrix(Evals[2])], sq_defects, 12, 3):


printf("Finding S-Hadamard SLA-parametrizations for de Launey BH(12,3) matrix ...\n");
deLH := Matrix([[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,2,1,2,1,2,2,1],[0,0,1,1,2,2,2,2,0,1,1,0],[0,0,1,1,0,1,0,1,2,2,2,2],[0,1,0,2,1,0,2,1,2,2,1,0],[0,2,1,2,0,2,1,0,2,0,1,1],[0,2,1,2,1,0,0,2,1,1,0,2],[0,1,0,2,2,1,1,0,0,1,2,2],[0,1,2,0,2,2,0,1,1,0,1,2],[0,2,2,1,2,0,1,1,0,2,0,1],[0,1,2,0,0,1,2,2,2,1,0,1],[0,2,2,1,1,1,2,0,1,0,2,0]]):
sq_defects := getSqDefectsForMats([deLH], 3):
MatsDeLauney := findSlaParametrizations([deLH], sq_defects, 12, 3):


# Redo the defect computation to get a new basis for the kernel, and try parametrizing again
evals_BH_12_3 := [Matrix(Evals[1]), Matrix(Evals[2]), deLH]:
getNewBasisAndParametrize := proc()
    local sq_defects, Mats;
    global evals_BH_12_3;
    sq_defects := getSqDefectsForMats(evals_BH_12_3, 3):
    Mats:= findSlaParametrizations(evals_BH_12_3, sq_defects, 12, 3):
end:
