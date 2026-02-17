
kernelopts(printbytes=false):

with(LinearAlgebra):
with(combinat):

read "../Block-Structure/utils.mpl":

# Checking if we can introduce one parameter into circulant BH(12,3) given by de Launey, 1992

a_vec_deL := [w,1,w^2,w,w^2,1]:
b_vec_deL := [w^2,1,1,w^2,w^2,w^2]:
n := 12:


generateCircMat := proc(a_vec::list, b_vec::list)
    local A, B;
    A := Matrix(6, shape = Circulant[a_vec]):
    B := Matrix(6, shape = Circulant[b_vec]):
    return(Matrix([[A, B], [B, A]])):
end:


generateParamVec := proc(a_comb::list, b_comb::list)
    local a_vec, b_vec, i, n;
    global a_vec_deL, b_vec_deL;
    n := nops(a_comb);
    a_vec := []; b_vec := [];
    for i to n do
        a_vec := [op(a_vec), a_vec_deL[i]*a_comb[i]];
        b_vec := [op(b_vec), b_vec_deL[i]*b_comb[i]];
    od;
    return(a_vec, b_vec);
end:


generateWordsFromList := proc(L::list, n::nonnegint)
    local result, shorterWords, word, letter;
    if n = 0 then
        return [[]]; 
    else
        shorterWords := procname(L, n-1);
        result := [];
        for word in shorterWords do
            for letter in L do
                result := [op(result), [op(word), letter]];
            end do;
        end do;
        return result;
    end if;
end proc:


checkHadConditions := proc(H::Matrix)
local i, j, ips, ip;
global n;
    ip := (x,y) -> numer(normal(rem(modW(add( x[i]/y[i] , i=1..n), w), 1+w+w^2, w))):
    ips := { seq( seq( ip(H[i],H[j]), j=i+1..n), i=1..n) };
    return( ips={0} );
end:


checkAllCombinations := proc(paramCombs::list)
local failCount, print_flag, a_comb, b_comb, M_par;
    failCount := 0: 
    print_flag := 10000:
    for a_comb in paramCombs do
        for b_comb in paramCombs do
            M_par := generateCircMat(
                generateParamVec(
                    a_comb, 
                    b_comb
                )
            );

            # Just check Hadamard conditions because there are many iterations
            if checkHadConditions(M_par) then
                printf("Found parameterization: %a, %a\n", a_comb, b_comb);
            else
                failCount := failCount + 1;
                if failCount > print_flag then
                    print_flag := print_flag + 10000;
                    printf("Current fail count %d, time %a\n", failCount, Now(ProcessClock));
                fi;
            fi;
        od:
    od:
    printf("FAIL COUNT: %d\n", failCount);
    printf("SUCCESS COUNT: %d\n", nops(paramCombs)^2 - failCount);
end:


# Check all combinations of binary parameter vectors
printf("Binary Parameter Vectors ...\n"):
checkAllCombinations(generateWordsFromList([1, a], 6));


# Try ternary parameter vectors, that is 1,a,1/a as parameters
# uncomment to try, but very can be very time consuming
# printf("Ternary Parameter Vectors ...\n"):
# checkAllCombinations(generateWordsFromList([1, a, a^(-1)], 6)):
