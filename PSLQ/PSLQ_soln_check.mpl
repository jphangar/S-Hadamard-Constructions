
# Main function for PSLQ based parameterization

kernelopts(printbytes=false):


with(NumberTheory):
with(LinearAlgebra):


reduceInnerProdByRtsu := proc(ips, rtsu)
local ip, cycloPoly, newIps, rt;
	newIps := ips;
	for rt in rtsu do
		cycloPoly := CyclotomicPolynomial(rt[1], rt[2]);
		newIps := { seq( rem( ip, cycloPoly, rt[2]), ip=newIps ) }:
	od:
	return(newIps);
end:


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
