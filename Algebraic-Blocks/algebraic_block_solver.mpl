kernelopts(printbytes=false):

with(LinearAlgebra):
with(NumberTheory):
with(combinat):

read "../Block-Structure/utils.mpl":
read "../Block-Structure/H12.m":


# Solve for parameter block for each possible permutation block of a block-structured parametrization of H12

permutation_matrices: # all 4x4 permutation matrices - from utils.mpl
n := 12:
Soln_Blks := []:

# Setup variable matrix X for parameter block
vars := {seq(x[i], i=1..16)}:
X := Matrix([[x[1], x[2], x[3], x[4]], [x[5], x[6], x[7], x[8]], [x[9], x[10], x[11], x[12]], [x[13], x[14], x[15], x[16]]]):


# Inline functions for computing inner product expressions
ip := (x,y) -> numer(normal(rem(modW(add( x[i]/y[i] , i=1..n), w), 1+w+w^2, w))):
sip := (x,y) -> numer(normal(rem(modW(add( x[i]/y[i] , i=1..n), w), 1+w+w^2, w))):

# Iterate over all permutation blocks and solve for parameter block
for P in permutation_matrices do
    # Generate parameter matrix in x[i] variables
    Mx := generateParameterMatrix(P, X);
    H12x := H12*~Mx;

    # Solve system of equations formed by inner products 
    ips := { seq( seq( ip(H12x[i],H12x[j]), j=i+1 ..n), i=1..n) } minus {0};
    sips := { seq( seq( sip(H12x[i],H12x[j]), j=i+1 ..n), i=1..n) } minus {0};
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
    N := normalizeMat(subs(op(validSolns), X), 4);
    trivialFlg := evalb(Equal(N, Matrix(4,4,1)));
    printf("Permutation block: %a, Valid solution count: %d, Trivial: %a\n", P, nops(validSolns), trivialFlg);

    # If solution is trivial, move to next iteration
    if trivialFlg then next; fi;

    # For nontrivial solution, confirm parametrization is S-Hadamard
    if not checkSHadConditions(H12*~generateParameterMatrix(P, N)) then
        printf("ERROR: solution matrix is not S-Hadamard\n");
        next;
    fi;

    # Add permutation and parameter block for S-Hadamard parametrizations to list 
    Soln_Blks := [op(Soln_Blks), [P, N]];
od:

save Soln_Blks, "algebraic_solutions.m":
