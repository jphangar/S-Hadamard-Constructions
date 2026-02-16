kernelopts(printbytes=false):

with(LinearAlgebra):

read "../Block-Structure/perm_and_parameter_blocks.m":
read "../Block-Structure/utils.mpl":
read "algebraic_constructions.mpl":

# Check if each parametrization found through numerical methods can be constructed by algebraically found parametrizations

for i from 1 to 79 do
    # Normalize the parameter block 
    K := normalizeMat(Blks[i][2], 4);
    P := Blks[i][1];
    
    # Check if P is in P1 or P2, this tells us which format K will match
    if Equal(P, P1[1]) or Equal(P, P1[2]) then
        # Matches N1
        K_N := subs({a=K[3,3], b=K[2,4], c=K[2,3]}, N1);
    elif Equal(P, P2[1]) or Equal(P, P2[2]) then
        # Matches N2
        K_N := subs({a=K[4,3], b=K[2,4], c=K[2,3]}, N2);
    else 
        printf("ERROR: Differs from expected permutation blocks: idx %d\n", i);
        next;
    fi;

    if not Equal(K_N, K) then
        printf("ERROR: Does not match expected parameter block form: idx %d\n", i);
    fi;
od:
