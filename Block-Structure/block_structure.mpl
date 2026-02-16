kernelopts(printbytes=false):

with(LinearAlgebra):
with(NumberTheory):
with(combinat):

read "pslq_valid_param_mats.m":
read "isomorphisms.txt": # Three objects, CT_mats, P_mats and Q_mats
read "utils.mpl":
read "H12.m": # H12


PMats := []: QMats := []: TransfMats := []:
ParMats := []: # List of Parameter Matrices 
Blks := []: # List of pairs of permutation and parameter blocks
n := 12:

# For each matrix, map to deL mat then determine P and K for coeff matrix
for i from 1 to 79 do
    PMats := [ op(PMats), Matrix(P_mats[i]) ];
    QMats := [ op(QMats), Matrix(Q_mats[i]) ];

    Hs := validMatrices[i][2];
    if i in CT_mats then
        # Take conjugate transpose of Hs
        Hs := Transpose(Hs^~(-1));
        Hs := subs({x[1]=1/x[1], x[2]=1/x[2], x[3]=1/x[3]}, Hs);
    fi;

    # Apply monomial matrices to Hs
    TransfMats := [ op(TransfMats), modOmega(PMats[i].Hs.QMats[i]) ];

    # Get coefficient matrix
    M := getParameterMatrix(TransfMats[i], n, w);
    ParMats := [op(ParMats), M];

    # Determine 3-cycle being applied to 4x4 block of coefficients
    P := findPermutationBlock(M);
    # Set parameter block to be first block in parameter matrix
    K := M[1..4,1..4];
    Blks := [op(Blks), [P, K]];

    if not Equal(TransfMats[i], H12*~generateParameterMatrix(P, K)) then
        printf("ERROR: index %a failed check\n", i);
    fi;
od:

save TransfMats, "block_structured_matrices.m":
save ParMats, "parameter_matrices.m":
save Blks, "perm_and_parameter_blocks.m":
