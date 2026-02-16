
# Function for trying PSLQ-based parameterization iterating through all possible combinations of matrix entries

with(LinearAlgebra):
with(combinat):

read "PSLQ_Param_Subroutine.mpl";

# H - numeric matrix to parameterize
# D - precision digits
# T - threshold for PSLQ
# minParams - minimum number of parameters to use
# maxParams - maximum number of parameters to use
tryParameterizeMatrix := proc(H::Matrix, D::integer, T::integer, minParams::integer, maxParams::integer)
local n, pts, found, paramLabels, i, j, numParams, comb, Hs, rtsu, paramsUsed, out;

    n := ColumnDimension(H);
    pts := [ seq(seq( [i,j] , j=2..n) , i=2..n) ];
    found := false;
    paramsUsed := 0;
    paramLabels := [ seq('x[i]', i=1..maxParams) ];
    
    for numParams from minParams to maxParams do

        # Generate all combinations
        comb := choose(pts, numParams);
        i := 1;
        while not(found) and (i <= nops(comb)) do

            out := parameterizeMatrix(H, [seq( [ op(comb[i][j]), paramLabels[j] ], j=1..numParams )], ['w', 'z', 'v', 'u'], D, T);

            if out<>FAIL then
                Hs, rtsu := out[1], out[2];
                found := true;
                paramsUsed := numParams;
            fi;
            
            i := i + 1;
        od;
    od;

    if not(found) then
        return(FAIL);
    fi;

    return(Hs, rtsu, paramsUsed);
end:
