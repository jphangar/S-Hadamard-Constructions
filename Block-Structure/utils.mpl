kernelopts(printbytes=false):

with(LinearAlgebra):
with(NumberTheory):
with(combinat):


modOmega := proc(H::Matrix)
    algsubs(w^3=1, algsubs(1/w=w^2, H));
end:

modW := proc(expr::algebraic)
    algsubs(w^3=1, algsubs(1/w=w^2, expr));
end:

normalizeMat := proc(H::Matrix, n::integer)
local M, i;
  M := Matrix(H);
  for i from 1 to n do
    M[i] := M[i]/M[i,1];
  od;
  for i from 1 to n do
    M[1..n,i] := M[1..n,i]/M[1,i];
  od;
  M;
end:


checkSHadConditions := proc(H::Matrix)
local i, j, ips, sips, ip, sip;
global n;
    ip := (x,y) -> numer(normal(rem(modW(add( x[i]/y[i] , i=1..n), w), 1+w+w^2, w))):
    sip := (x,y) -> numer(normal(rem(modW(add( x[i]^2/y[i]^2 , i=1..n), w), 1+w+w^2, w))):
    ips := { seq( seq( ip(H[i],H[j]), j=i+1..n), i=1..n) };
    sips := { seq( seq( sip(H[i],H[j]), j=i+1..n), i=1..n) };
    return( ips={0} and sips={0} );
end:


# Function for getting the corresponding parameter matrix from a parameterized matrix 
getParameterMatrix := proc(M::Matrix, n::integer, w::algebraic)
local i, j, T;
    T := Matrix(n):
    for i from 1 to n do
        for j from 1 to n do
            T[i,j] := lcoeff(M[i,j], w):
        od:
    od:
    T;
end:


# Generate parameter matrix using permutation matrix P and parameter block K
generateParameterMatrix := proc(P::Matrix, K::Matrix)
    return(Matrix([
        [K, K.P, K.(P^2)],
        [P.K, P.K.P, P.K.(P^2)],
        [(P^2).K, (P^2).K.P, (P^2).K.(P^2)]
    ]));
end:


# Generating all 4x4 permutation matrices
permutations_list := permute(4):
permutation_matrices := []:
for p in permutations_list do
    P := Matrix(4, 4, 0); # Initialize an n x n zero matrix
    for i from 1 to 4 do
        P[i, p[i]] := 1; # Place a 1 at (i, p[i])
    end do;
    permutation_matrices := [op(permutation_matrices), P];
end do:


# Function for determining 4x4 permutation block of a parameter matrix
findPermutationBlock := proc(T::Matrix)
local permMat;
global permutation_matrices;
    for permMat in permutation_matrices do
        if Equal(generateParameterMatrix(permMat, T[1..4,1..4]), T) then
            return permMat;
        fi;
    od;
    return FAIL;
end:
