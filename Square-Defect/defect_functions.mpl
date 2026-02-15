kernelopts(printbytes=false):

with(LinearAlgebra):
with(NumberTheory):

# Functions for computing the defect and square defect of S-Hadamard matrices,
# and util functions for building HLA/SLA-parametrizations 

# Builds a BH(n,q) matrix from its integer form
IntMatToBH := proc(H::Matrix, q::integer)
local n, i, j, M, w;
	n := RowDimension(H);
	M := Matrix(n,n);
	w := cos(2*Pi/q) + I*sin(2*Pi/q);
	for i from 1 to n do
	for j from 1 to n do
		M[i,j] := w^H[i,j];
	od;
	od;
	M;
end:


# Build symbolic matrix from integer form of a BH matrix and an affine 
# parametrization given by matrix R, the root of unity is denoted by w
BuildSymbolicMatrix := proc(H::Matrix, R::Matrix)
local n, i, k, M, j;
    n := RowDimension(H);
    M := Matrix(n,n);
    for i from 1 to n do
    for j from 1 to n do
        M[i,j] := w^H[i,j] * exp(I*R[i,j]);
    od;
    od;
    M;
end:


# Check that the given symbolic BH(n,q) matrix is S-Hadamard, where w denotes the root of unity; this may be a parametrized matrix
CheckSHadConditions := proc(H::Matrix, n::integer, q::integer)
local i,j,ips,sips,f,k;
    f := CyclotomicPolynomial(q, w);
    ips := { seq( seq( rem( numer( normal( add(expand( H[i,k]/H[j,k] ), k=1..n) ) ) , f, w) , j=i+1..n) , i=1..n) };
    sips := { seq( seq( rem( numer( normal( add(expand( H[i,k]^2/H[j,k]^2 ), k=1..n) ) ) , f, w) , j=i+1..n) , i=1..n) };
    return( ips={0} and sips={0} );
end:


# Return a normalized version of the given matrix
NormalizeMat := proc(H::Matrix, n::integer)
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


# Generate an (n x n) matrix variable R, where variable names are r[i,j]
Gen_R := proc(n::integer)
local i, j, R, vars;
	R := Matrix(n,n); vars := {};
	for i from 1 to n do
		for j from 1 to n do
		R[i,j] := r[i,j]; vars := {op(vars), r[i,j]};
		od;
	od;
	R, vars;
end:


# Generate differential of row inner product expressions
Gen_ips_eq := proc(H::Matrix, R::Matrix, n::integer)
local i, j, k, ips, ip, re_im_ips, im_ip, re_ip;
	ips := { seq( seq( simplify(numer(factor(add( H[i,k]/H[j,k] * (R[i,k] - R[j,k]) , k=1..n)))), j=i+1..n), i=1..n ) };
	# Split into real and imaginary parts for each expression
	re_im_ips := {}:
	for ip in ips do
		im_ip := coeff(ip, I);
		re_ip := expand(ip-I*im_ip);
		re_im_ips := {op(re_im_ips), im_ip=0, re_ip=0};
	od:
	re_im_ips;
end:


# Generate differential of dephasing expressions
Gen_dphs_eq := proc(R::Matrix, n::integer)
local i,j;
  	{ seq( R[1,j]=0, j=2..n), seq( R[i,1]=0, i=1..n) };
end:


# Generate differential of squared row inner product expressions
Gen_sips_eq := proc(H::Matrix, R::Matrix, n::integer)
local i, j, k, sips, sip, re_im_sips, im_sip, re_sip;
	sips := { seq( seq( simplify(numer(factor(add( H[i,k]^2/H[j,k]^2 * (R[i,k] - R[j,k]) , k=1..n)))), j=i+1..n), i=1..n ) };
	# Split into real and imaginary parts for each expression
	re_im_sips := {}:
	for sip in sips do
		im_sip := coeff(sip, I);
		re_sip := expand(sip-I*im_sip);
		re_im_sips := {op(re_im_sips), im_sip=0, re_sip=0};
	od:
	re_im_sips;
end:


# Compute the defect of the given complex Hadamard matrix
ComputeDefect := proc(M::Matrix)
local n, H, R, vars, ips, dps, soln, free_vars, var_sol;
	n := RowDimension(M);
	H := NormalizeMat(M, n);
	R, vars := Gen_R(n);
	ips := Gen_ips_eq(H, R, n);
	dps := Gen_dphs_eq(R, n);

	soln := solve(ips union dps, vars);

	# Count free vars in solution to determine dimension of solution
	free_vars := {};
	for var_sol in soln do
		if nops({op(var_sol)})=1 then
			free_vars := {op(free_vars), op(var_sol)[1]};
		fi;
	od;
	# Return defect value, matrix variable defining kernel, variables
	return nops(free_vars), subs(soln, R), free_vars;
end:


# Compute the square defect of the given S-Hadamard matrix
ComputeDefectSq := proc(M::Matrix)
local n, H, R, vars, ips, sips, dps, soln, free_vars, var_sol;
	n := RowDimension(M);
	H := NormalizeMat(M, n);
	R, vars := Gen_R(n);
	ips := Gen_ips_eq(H, R, n);
	sips := Gen_sips_eq(H, R, n);
	dps := Gen_dphs_eq(R, n);

	soln := solve(ips union dps union sips, vars);

	# Count free vars in solution to determine dimension of solution
	free_vars := {};
	for var_sol in soln do
		if nops({op(var_sol)})=1 then
			free_vars := {op(free_vars), op(var_sol)[1]};
		fi;
	od;
	# Return sq defect value, matrix variable defining kernel, variables
	return nops(free_vars), subs(soln, R), free_vars;
end:
