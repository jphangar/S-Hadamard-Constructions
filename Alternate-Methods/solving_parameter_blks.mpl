kernelopts(printbytes=false):

with(LinearAlgebra):
with(NumberTheory):
with(combinat):

read "../Block-Structure/utils.mpl":
read "../SBH-Matrices/SBH-12-3-mpl.txt":

permutation_matrices: # all 4x4 permutation matrices - from utils.mpl
n := 12:

# Solve for parameter blocks of BH(12,3) matrices in Aalto repository


IntToSymbolicMat := proc(M::list, n::integer)
local i,j,H;
    H := Matrix(n,n);
    for i from 1 to n do
    for j from 1 to n do
        H[i,j] := w^M[i,j];
    od;
    od;
    H;
end:

# BH(12,3) matrices from Aalto repository in symbolic form
AaltoMats := [IntToSymbolicMat(Evals[1], n), IntToSymbolicMat(Evals[2], n)]:

# Setup variable matrix X for parameter block
vars := {seq(x[i], i=1..16)}:
X := Matrix([[x[1], x[2], x[3], x[4]], [x[5], x[6], x[7], x[8]], [x[9], x[10], x[11], x[12]], [x[13], x[14], x[15], x[16]]]):

# Inline functions for computing inner product expressions
ip := (x,y) -> numer(normal(rem(modW(add( x[i]/y[i] , i=1..n), w), 1+w+w^2, w))):
sip := (x,y) -> numer(normal(rem(modW(add( x[i]/y[i] , i=1..n), w), 1+w+w^2, w))):


solveForParameterBlocks := proc(H::Matrix)
global permutation_matrtices, X, vars, n;
local P, Mx, Hx, ips, sips, solns, validSolns, testfxn, i, j, soln, N, trivialFlg;
    # Iterate over all permutation blocks and solve for parameter block
    for P in permutation_matrices do
        # Generate parameter matrix in x[i] variables
        Mx := generateParameterMatrix(P, X);
        Hx := H*~Mx;

        # Solve system of equations formed by inner products 
        ips := { seq( seq( ip(Hx[i],Hx[j]), j=i+1 ..n), i=1..n) } minus {0};
        sips := { seq( seq( sip(Hx[i],Hx[j]), j=i+1 ..n), i=1..n) } minus {0};
        solns := solve(ips union sips, vars);

        # Determine which solutions are valid
        validSolns:= {}:
        testfxn := mul(x[i], i=1..16):
        for soln in solns do
        if eval(testfxn, soln)<>0 then
            validSolns := validSolns union {soln};
        fi;
        od:
        
        # If there's no valid solution, move to next iteration
        if nops(validSolns) = 0 then next; fi;

        # Determine if solution is trivial by normalizing
        for soln in validSolns do
            N := normalizeMat(subs(soln, X), 4);
            trivialFlg := evalb(Equal(N, Matrix(4,4,1)));
            printf("Permutation block: %a, Valid solution count: %d, Trivial: %a\n", P, nops(validSolns), trivialFlg);
        od;
    od:
end:

printf("Parameter block solving for H1 ...\n"):
solveForParameterBlocks(AaltoMats[1]):

printf("Parameter block solving for H2 ...\n"):
solveForParameterBlocks(AaltoMats[2]):
