
# Function for trying PSLQ-based parameterization using random combinations of matrix entries

with(LinearAlgebra):

read "PSLQ_Param_Subroutine.mpl";

# H - numeric matrix to parameterize
# D - precision digits
# T - threshold for PSLQ
# trails - maximum number of parameterization attempts per number of parameters to use
# minParams - minimum number of parameters to use
# maxParams - maximum number of parameters to use
tryParameterizeMatrixRand := proc(H::Matrix, D::integer, T::integer, trials::integer, minParams::integer, maxParams::integer)
local n, found, paramLabels, i, j, numParams, comb, Hs, rtsu, paramsUsed, out, idx;

    n := ColumnDimension(H);
    found := false;
    paramsUsed := 0;
    paramLabels := [ seq('x[i]', i=1..maxParams) ];
    idx := rand(1..n);
    
    for numParams from minParams to maxParams do

        i := 1;
        while not(found) and (i <= trials) do
            # Generate combination of entries randomly
            comb := [ seq( [idx(), idx()], numParams ) ];
        
            out := parameterizeMatrix(H, [seq( [ op(comb[j]), paramLabels[j] ], j=1..numParams )], ['w', 'z', 'v', 'u'], D, T);

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
