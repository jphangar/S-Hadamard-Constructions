kernelopts(printbytes=false):

with(LinearAlgebra):


# Returns a table of the Haagerup multiset for the given BH(n,3) matrix
HaagerupInvariant := proc(H::Matrix)
local n, i, j, k, l, InvSet, val;
    n := RowDimension(H);
    InvSet := table();
    for i from 1 to n do
    for j from 1 to n do
    for k from 1 to n do
    for l from 1 to n do
        val := algsubs( w^3=1, algsubs( 1/w=w^2, H[i,j]*H[k,l]/(H[i,l]*H[k,j]) ) );
        if not assigned(InvSet[val]) then
            InvSet[val] := 0;
        fi;
        InvSet[val] := InvSet[val] + 1;
    od;
    od;
    od;
    od;
    return InvSet;
end:


# True if invariants are equal, otherwise false
CheckHaagerupInv := proc(t1::table, t2::table)
local idx;
    if ({indices(t1)} minus {indices(t2)}) <> {} or ({indices(t2)} minus {indices(t1)}) <> {} then
        return false;
    else # entries are same, go through and check multiplicies
        for idx in {indices(t1)} do
            if t1[idx[1]] <> t2[idx[1]] then
                return false;
            fi;
        od;
    fi;
    return true;
end:
