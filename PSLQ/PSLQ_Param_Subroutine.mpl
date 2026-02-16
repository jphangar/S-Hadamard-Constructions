
# Subroutine for PSLQ-based parameterization of matrix

# Note that we want the matrix to be fully symbolic in our implementation, so we use labels
# for roots of unity - this differs from the algorithm description given in the thesis

with(LinearAlgebra):
with(IntegerRelations):

parameterizeMatrix := proc(H::Matrix, paramPositions::list(list), rootLabels::list, D::integer, T::integer)
local Hp, rtsAssigned, n, paramArgs, i, j, res, x, rootName;
	Digits := D;
	rtsAssigned := [];
	n := ColumnDimension(H);
	Hp := Matrix(n,n);
	paramArgs := seq( argument(H[paramPositions[i][1], paramPositions[i][2]]), i=1..nops(paramPositions) );

	for i to n do
		for j to n do

			res := PSLQ( evalf( [-argument(H[i,j]), 2*Pi, paramArgs], D ) );

			if res[1] < 0 then res := -res fi;

			if (add(x^2, x=res) > T) or (res[1] = 0) then
				return FAIL;
			else
				if res[1] = 1 then rootName := 1;
				elif res[1] = 2 then rootName := -1;
				else
					if member(res[1], rtsAssigned, 'pos') then
						rootName := rootLabels[pos];
					elif nops(rtsAssigned) < nops(rootLabels) then
						rtsAssigned := [op(rtsAssigned), res[1]];
						rootName := rootLabels[nops(rtsAssigned)];
					else
						return FAIL;
					fi;
				fi;

				Hp[i,j] := rootName^(res[2]) * mul( paramPositions[i][3]^(res[i+2]/res[1]), i=1..nops(paramPositions) );

			fi;
		od;
	od;

	return(Hp, [ seq( [rtsAssigned[i], rootLabels[i]], i=1..nops(rtsAssigned) ) ]);

end:
