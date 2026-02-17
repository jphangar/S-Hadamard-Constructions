kernelopts(printbytes=false):

with(LinearAlgebra):
with(combinat):

read "haagerup_invariant.mpl":
read "../Automorphisms/matrix_constructions.mpl":


# Check that Haagerup Invariant of all four algebraically obtained matrices is equal
t1 := HaagerupInvariant(Mats[1]):
for i from 2 to 4 do
    t2 := HaagerupInvariant(Mats[i]);
    if not CheckHaagerupInv(t1, t2) then
        printf("ERROR: Haagerup invariant of matrices differs\n");
    fi;
od:

# save t1 as this is invariant value of Habc
save t1, "Habc-Haagerup-Inv.m":

# uncomment to print invariant value
#for idx in {indices(t1)} do
#    printf("%a, %a\n", idx[1], t1[idx[1]]);
#od:


printf("--- PARAMETER VARIANTS WITH EQUAL HAAG. INV. ---- \n"):

# 1 = don't invert, -1 = invert (these are all binary strings of length 3)
inverse_choices := [[1,1,1], [1,1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,-1], [1,-1,1], [-1,-1,1], [-1,1,1]]:
vals := [a,b,c]:
par_choices := permute([a,b,c,a*b/c], 3):
for perm in par_choices do
    for inverse_choice in inverse_choices do
        param_choice := [ seq(perm[i]^inverse_choice[i], i=1..3) ];
        N1_perm := eval(N1, { seq(vals[i]=param_choice[i], i=1..3) });
        H_perm := H12*~generateParameterMatrix(P1[1], N1_perm);

        if CheckHaagerupInv(t1, HaagerupInvariant(H_perm)) then
            printf("%a\n", param_choice);
        fi;
    od:
od:
