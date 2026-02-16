///////////////////////////////////////////
//
// Basic functions for manipulating Butson Hadamard (BH) 
// matrices and generalised Hadamard (GH) matrices.
//
// Section 1 - Transformation functions (GH, BH, IntMatrix)
// Section 2 - Construction of Expanded matrices
// Section 3 - Automorphism group computations
// Section 4 - Equivalence testing
//
// Important functions, throughout H is a GH over G.
// GHtoBH(H, G) - gives a BH
// BHtoGH(B, G) - B Butson, constructs H over G
// GHtoExpMatrix(H, G) - constructs the expanded matrix
// GHAutomorphismGroup - constructs the group of linear 
// automorphisms of H, requires 'nauty'.
//
////////////////////////////////////////////


//////////////////////////////////////////////
//
// Section 1
//
// Over p^th roots, BHs and GHs co-incide. For convenience, 
// we typically work with GHs, and store matrices as integer 
// matrices which encode the  exponents of a BH or GH. We provide 
// all relevant conversion functions here.
//
// (Note that the function GH->BH is non-canonical since relies 
// on a choice of faithful character of G, thus the composition 
// GH->BH->GH may not be the identity. Furthermore, the result of 
// this composition will not necessarily be Hadamard equivalent to 
// the original matrix - it will differ from the original by a 
// 'global equivalence operation', i.e. an automorphism of the cyclotomic 
// field applied identically to all entries of the matrix.
//
// The functions all take the form 
// XtoY where X and Y are one of 'BH', 'GH', 'IntMatrix'
// (and one of X or Y is 'GH' - the others are obtained as compositions)
// The functions all take two arguments: 
// H - a matrix of type X
// G - a group. If X is 'GH' this is the alphabet group of H, otherwise it is 
// the alphabet group that the returned matrix should be written over.
//
////////////////////////////////////////////
 
MatrixOfOnes := function (U, n)
   R := GroupAlgebra (Rationals (), U);
   MA := MatrixRing (R, n);
   one := Id (U);
   list := [one: i in [1..n^2]];
   return MA!list;
end function;

GenOfCyclicGroup := function(G);

	for i in G do;
		if Order(i) eq Order(G) then;
		gen := i;
		break;
		end if;
	end for;
	return gen;
end function;

GHtoIntMatrix := function(M,U);

	gen := GenOfCyclicGroup(U);

	n := NumberOfRows(M);
	Z := ZeroMatrix(Rationals(),n,n);
	index := [gen^i : i in [0..#U-1]];
	for i in [1..n] do;
		for j in [1..n] do;
			Z[i][j] +:= (Position(index,M[i][j])-1);
			end for;
		end for;
	
	return Z;

end function;



IntMatrixtoGH := function(M,U);

	gen := GenOfCyclicGroup(U);

	n := NumberOfRows(M);
	Z := MatrixOfOnes(U,n);
	for i in [1..n] do;
		for j in [1..n] do;
			Z[i][j] := Z[i][j] * gen^(M[i][j]);
			end for;
		end for;
	
	return Z;

end function;


GHtoBH := function(H, A);

	n := Exponent(A);
	F<z> := CyclotomicField(n);

	for i in CharacterTable(A) do;
		if IsFaithful(i) then;
			Chi := i;
			break;
		end if;
	end for;

	list := [];
	n := NumberOfRows(H);

	for i in [1..n] do;
		for j in [1..n] do;
			Append(~list, F!Chi(A!H[i,j]));
		end for;
	end for;

	return Matrix(n,n,list);

end function;


BHtoGH := function(H, G);

	gen := GenOfCyclicGroup(G);

	d := NumberOfRows(H);
	Entries := {H[i,j]: i,j in [1..d]};
	n := #Entries;
	l := CyclotomicOrder(Parent(H[1,1]));
	if n eq l then;
	w := Parent(H[1,1]).1;
	else 
	w := -1*Parent(H[1,1]).1;
	end if;

	Z := MatrixOfOnes(G, d);

	index := [w^i : i in [0..n-1]];

	if (Set(Entries) ne Set(index)) or (#index ne Order(gen)) then;
	print "error, matrix entries cannot be written over input group" ;
	       	return false;
	end if;
	
	for i in [1..d] do;
		for j in [1..d] do;
			Z[i][j] := Z[i][j] * gen^(Position(index,H[i][j])-1);
			end for;
		end for;
	
	return Z;
end function;


////////////////////////////////////////////////////
//
// Section 2: Expanded matrices
//
// Given H, G where H is a GH and G is the group of matrix entries,
// we construct the expanded matrix. Note that the rows are ordered 
// so that the image of the first row of H consists of rows 1, |G|+1, 
// 2*|G|+1,... of EH. There is also an inverse function.
//
/////////////////////////////////////////////////

KroneckerDelta := function(x, G);
	if (x eq G!1) eq true then;
		return 1;
	else;
		return 0;
	end if;
end function;

KroneckerMatrix := function(Mat, G);

	n := NumberOfRows(Mat);
	List := [KroneckerDelta(Mat[i,j], G): i in [1..n], j in [1..n]];

	return Matrix(Rationals(),n,List);
end function;
	
CayleyMatrix := function(G);

//	This only works for cyclic groups. 
// 	We use the (obvious) fact that the 
// 	Cayley matrix can be written in 
// 	circulant form.

	RG := GroupAlgebra(Rationals(), G);
	g := GenOfCyclicGroup(G);
	n := Order(G);
	List := [RG!g^i: i in [0..n-1]];
	// return Circulant(List); // Circulant() is not defined
    n := #List;
    return Matrix([[List[(i + j - 1) mod n + 1] : j in [1..n]] : i in [1..n]]); 
end function;


GHtoExpMatrix := function(H, G);

	Cay := CayleyMatrix(G);
	Exp := KroneckerProduct(Cay, H);

	return KroneckerMatrix(Exp, G);
end function;


ExpMatrixtoGH := function(E, G);

	nk := NumberOfRows(E);
	k := Order(G);
	n := Integers()!(nk/k);

	List := [i: i in G];

	GrpAlg := GroupAlgebra(Rationals(), G);
	Z := Zero(MatrixAlgebra(GrpAlg, n));

	for i in [0..(k-1)] do;
		T := Submatrix(E, 1, i*n+1, n, n);
		g := GrpAlg!List[i+1];
		Z := Z + g*T;
	end for;

	return Z;
end function;


//////////////////////////////////////////////////////
// 
// Section 3 - Automorphism group computations
//
// We construct a graph from the expanded matrix of H and 
// compute its automorphism group with nauty. Automorphisms 
// of the graph are not necessarily automorphisms of H in 
// the usual sense, so we compute the intersection of the 
// automorphism group of the graph with the group of equivalence 
// operations on the set of generalised Hadamard matrices isomorphic 
// to H.
//
// The default returned value of the automorphism group function 
// is a permutation group. We also provide a function for 
// transforming elements of this group into pairs of monomial matrices.
//
////////////////////////////////////////////////////////


// GHEquivalenceGroup gives the group of equivalence operations 
// on H. Reorderpoints ensures that the blocks of imprimitivity 
// of this group correspond to image sets of the rows of H in the 
// expanded matrix. (See the comment on the ordering of rows in 
// Section 2.
ReorderPoints := function(n,k);
	arith := function(i);
		return Ceiling(i/k) + n*((i-1) mod k);
	end function;
	list := [arith(i): i in [1..n*k]];
	return SymmetricGroup(n*k)!list;
end function;

GHEquivalenceGroup := function(H,G);

	n := NumberOfRows(H);
	Base := CyclicGroup(Order(G));
	Top := SymmetricGroup(n);
	AutParent := WreathProduct(Base,Top);
	return AutParent^ReorderPoints(n,Order(G));
end function;

GHAutomorphismGroup := function(H, G);

	EH := GHtoExpMatrix(H, G);
	n := NumberOfRows(H);
	k := Order(G);
	IncStr := IncidenceStructure<n*k| EH>;

// 	We compute the action of the automorphism group 
// 	on rows and columns of H as the action on points 
// 	and blocks of incidence structure of the associated 
// 	matrix.
	P := PointGroup(IncStr);
	B := BlockGroup(IncStr);
	GensP := #Generators(P);
	Gens := [];

//	For convenience we package the automorphisms
//	in a single (intransitive) group. The first 
// 	orbit gives the row action and the second the 
// 	column action.
	D, homs := DirectProduct(P, B);
	for i in [1..GensP] do;
		pointaut := homs[1](P.i);
		blockaut := homs[2](B.i);
		Append(~Gens, pointaut*blockaut);
	end for;

	BigAut := sub<D | Gens>;
	Parent := DirectProduct(GHEquivalenceGroup(H, G), 
				GHEquivalenceGroup(H, G));

	return BigAut meet Parent;
end function;

// We give a function taking a permutation of the graph of H 
// and returning a pair of monomial matrices corresponding to the 
// mapping. (This need not be an automorphism of H, it can be any 
// element of the image of the group of equivalence operations.

GHPermAutToAut := function(perm, G);

	T := Degree(Parent(perm));
	nk := Integers()!(T/2);
	Gp := PermutationGroup<T|perm>;
	Orb1 := {1..nk};
	Orb2 := {nk+1..T};

	RowAction := OrbitAction(Gp, Orb2);
	ColAction := OrbitAction(Gp, Orb1);
	Pmat := IdentityMatrix(Rationals(), nk)^RowAction(perm);
	Qmat := IdentityMatrix(Rationals(), nk)^ColAction(perm);

	P := ExpMatrixtoGH( Pmat, G);
	Q := ExpMatrixtoGH( Qmat, G);

	return [P, Q];
end function;


//////////////////////////////////////////////
//
// Section 4 - Equivalence testing
//
// Essentially we construct graphs for H_1 and H_2 
// as in the automorphism group computation and 
// use nauty to find a permutation mapping H_1 to H_2. 
// Such a permutation need not correspond to a 
// valid equivalence operation on the underlying 
// GHs, so we test this explicitly.
//
//////////////////////////////////////////////


// This function computes the intersection of the 
// coset gH with K, where H and K are permutation 
// groups of the same degree. 
CosetIntersection := function(g, H, K);

	Coset := H*g;
	Intersection := H meet K;
	Trans := Transversal(H, Intersection);
	Ans := false;
	for i in Trans do;
		if i*g in K then;
		Ans := true;
		elt := i*g;
		end if;
	end for;

	if Ans eq true then;
	return Ans, elt;
	else
	return Ans, g;
	end if;
	
end function;

GHEquivalence := function(H1,H2,G1,G2);

// First we build bipartite graphs encoding the 
// expanded matrices of each matrix.
// We colour the vertices so that rows cannot 
// be interchanged with columns.

	E1 := GHtoExpMatrix(H1, G1);
	E2 := GHtoExpMatrix(H2, G2);
	Z := Zero(Parent(E1));
	nk := NumberOfRows(E1);
	k := Order(G1);
	n := Integers()!(nk/k);

	M1 := BlockMatrix(2, 2, [Z, E1, Transpose(E1), Z]);
	M2 := BlockMatrix(2, 2, [Z, E2, Transpose(E2), Z]);
	Gph1 := Graph<2*nk|M1>;
	Gph2 := StandardGraph(Graph<2*nk|M2>);
	colours := [1: i in [1..nk]] cat [2: i in [1..nk]];
	AssignLabels(VertexSet(Gph1), colours);
	AssignLabels(VertexSet(Gph2), colours);

	bool, sigma := IsIsomorphic(Gph1, Gph2);

	if bool eq true then;


// sigma is returned as a mapping on vertices. We convert into a permutation.
	BigSym := SymmetricGroup(2*n*k);
	list := [Index(sigma(i)): i in [1..2*nk]];
	sigma := BigSym!list;

// Finally, we search for a valid equivalence operation of the underlying 
// Hadamard matrices mapping H1 to H2.

	BigAut := AutomorphismGroup(Gph1);

//	BigAut is the full automorphism group of the graph. 
//	Not every automorphism corresponds to an automorphism 
//	of the underlying Butson matrix. Those that do lie in the 
//	following group.

	RowAuts := GHEquivalenceGroup(H1, G1);
	ProperAuts := DirectProduct(RowAuts, RowAuts);

//	We attempt to find a valid automorphism. 
//	This may fail, in which case the graphs 
//	are isomorphic but the matrices are inequivalent.
//	(This may be due to a global equivalence operation.)

	properautbool,properaut:= CosetIntersection(sigma, BigAut, ProperAuts);
	if properautbool eq false then;
		return false, properaut;
	else;
		return true, properaut;
	end if;

	else;
	return false;
	end if;
	
end function;
