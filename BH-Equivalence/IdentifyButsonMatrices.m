// Identify all omega-evaluations of parameterized matrices 

load "pkg/ButsonFunctions.m";

load "BH-12-3-Catalog.m";
load "omega_evaluations.txt";

G<w> := SmallGroup(3,1);

inverse := function(x)
    if x eq 0 then
        return 0;
    else
        return x^-1;
    end if;
end function;

// BH(12,3) constructed by de Launey
deL_mat := IntMatrixtoGH(BH12_3_1989[1], G);

for i in [1..79] do
    // Check all 27 omega-evaluations for this example
    deLIsoCount := 0; conjTrpIsoCount := 0;
    for j in [1..27] do
        H := IntMatrixtoGH(Evals[27*(i-1)+j], G);
        b, sigma := GHEquivalence(deL_mat, H, G, G);
        if b then
            deLIsoCount := deLIsoCount + 1;
        else
            // try conjugate transpose
            H := Transpose(Matrix([[inverse(H[i, j]) : j in [1..Ncols(H)]] : i in [1..Nrows(H)]]));
            b, sigma := GHEquivalence(deL_mat, H, G, G);
            if b then 
                conjTrpIsoCount := conjTrpIsoCount + 1;
            else
                printf "ERROR: matrix %o, evaluation %o is not isomorphic to any BH(12,3)\n", i, j;
            end if;
        end if;
    end for;

    // Check all omega-evaluations were isomorphic to the same unique BH(12,3)
    if Abs(deLIsoCount - conjTrpIsoCount) ne 27 then
        printf "ERROR: matrix %o has %o evaluations isomorphic to deLauney, and %o evaluations isomorphic to conjugate transpose\n", i, deLIsoCount, conjTrpIsoCount;
    end if;
end for;
