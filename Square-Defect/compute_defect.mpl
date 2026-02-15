kernelopts(printbytes=false):

with(NumberTheory):

read "defect_functions.mpl":

# Compute square defect for each BH(n,q) in the given list, and print results
# The first list returned contains the results of the defect computation, 
# and the second contains only the matrices that had an S-Hadamard full SLA-parametrization
getSqDefectsForMats := proc(evals::list(list), q::integer)
local i, out_dft, out_sla, dfct, Rmat, freeVr, H, isFullSla, M;
    out_dft := []; out_sla := [];
    printf("%-10s %-15s %-10s\n", "Index", "Square Defect", "Full SLA"):
    for i from 1 to nops(evals) do
        isFullSla := false;
        H := Matrix(evals[i]);
        # Compute square defect 
        dfct, Rmat, freeVr := ComputeDefectSq( IntMatToBH(H, q) );
        # Check if we have an S-Hadamard full SLA-parametrization
        M := BuildSymbolicMatrix(H, Rmat);
        if CheckSHadConditions(M, RowDimension(M), q) then
            isFullSla := true;
            out_sla := [op(out_sla), [i, H, Rmat, freeVr]];
        fi;
        printf("%-10d %-15d %-10a\n", i, dfct, isFullSla);
        out_dft := [op(out_dft), [dfct, Rmat, freeVr]];
    od;
    out_dft, out_sla;
end:


# Computing square defect and checking full SLA-parametrization
nqPairs := [[12,3], [12,6], [12,12], [16,4], [18,3]]:
for nqPair in nqPairs do
    read sprintf("../SBH-Matrices/SBH-%d-%d-mpl.txt", nqPair[1], nqPair[2]);
    printf("Square Defect Calculations for BH(%d, %d) Matrices ...\n", nqPair[1], nqPair[2]);
    sq_defects, full_slas := getSqDefectsForMats(Evals, nqPair[2]);

    save sq_defects, sprintf("output/sbh-%d-%d-defects.m", nqPair[1], nqPair[2]);
    save full_slas, sprintf("output/sbh-%d-%d-parametrizations.m", nqPair[1], nqPair[2]);
od:


# Compute the defect for each BH(n,q) matrix in the given list, and print results
getDefectsForMats := proc(evals::list(list), q::integer)
local i, dfct, Rmat, freeVr, H;
    printf("%-10s %-15s\n", "Index", "Defect"):
    for i from 1 to nops(evals) do
        H := Matrix(evals[i]);
        # Compute defect 
        dfct, Rmat, freeVr := ComputeDefect( IntMatToBH(H, q) );
        printf("%-10d %-15d\n", i, dfct);
    od;
end:


# Get defect for BH(n,q) matrices where q is not prime
nqPairs := [[12,6], [12,12], [16,4]]:
for nqPair in nqPairs do
    read sprintf("../SBH-Matrices/SBH-%d-%d-mpl.txt", nqPair[1], nqPair[2]);    
    printf("Defect Calculations for BH(%d, %d) Matrices ...\n", nqPair[1], nqPair[2]);
    getDefectsForMats(Evals, nqPair[2]);
od:
