kernelopts(printbytes=false):

with(NumberTheory):
with(LinearAlgebra):

read "PSLQ_Try_Param_Rand.mpl";
read "PSLQ_Try_Param_Naive.mpl";

# Main function for running PSLQ-based parameterization


# Reading order 12 S-Hadamard approximations found by optimizer
read "order-12-optimizer-solns.m"; 
soln_mats: # list of 79 numeric order 12 S-Hadamard matrices to parametrize


# Helper function for checking if a parameterization found is S-Hadamard
reduceInnerProdByRtsu := proc(ips::set, rtsu::list)
local ip, cycloPoly, newIps, rt;
	newIps := ips;
	for rt in rtsu do
		cycloPoly := CyclotomicPolynomial(rt[1], rt[2]);
		newIps := { seq( rem( ip, cycloPoly, rt[2]), ip=newIps ) }:
	od:
	return(newIps);
end:


# Check if a parametrization found is S-Hadamard
checkSHadConditions := proc(Hs::Matrix, rtsu::list)
local n, i, j, k, ips, sips;
	n := ColumnDimension(Hs);
	ips := { seq( seq( numer(normal(add( Hs[i,k]/Hs[j,k], k=1..n ))), j=1..n ), i=1..n ) }:
	ips := reduceInnerProdByRtsu(ips, rtsu);
	if (ips minus {n}) <> {0} then
		return(false, "NOT_HADAMARD");
	fi;
	sips := { seq( seq( numer(normal(add( Hs[i,k]^2/Hs[j,k]^2, k=1..n ))), j=1..n ), i=1..n ) }:
	sips := reduceInnerProdByRtsu(sips, rtsu);
	if (sips minus {n}) <> {0} then
		return(false, "NOT_S_HADAMARD");
	fi;
	return(true, "S_HADAMARD");
end:


# Try to parametrize each approximation with 3 parameters using the randomized approach
validMatrices := []:
printf("%-10s %-15s %-10s %-10s\n", "Index", "Parameters", "Roots", "S-Hadamard?"):

for i from 1 to nops(soln_mats) do
	result := tryParameterizeMatrixRand(soln_mats[i], 10, 50, 500, 3, 3);
	if (result=FAIL) then
		printf("%-10d %-15d %-10a %-10a\n", i, 3, [], FAIL);
	fi;
	isValid, msg := checkSHadConditions(result[1], result[2]);
	printf("%-10d %-15d %-10a %-10a\n", i, 3, result[2], isValid);
	if isValid then
		validMatrices := [op(validMatrices), [i, result]];
	fi;
od:

# Save valid parametrized matrices
save validMatrices, "pslq_valid_param_mats.m";
