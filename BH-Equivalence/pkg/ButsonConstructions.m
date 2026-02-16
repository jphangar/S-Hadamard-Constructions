/////////////////////////////////////////////////////
// 
// We give code to construct the McNulty-Wiegert and 
// Butson matrices. The equivalence test in ButsonFunctions.m
// can then be used to check the equivalence of these 
// matrices at small orders. ButsonFunctions.m should 
// be loaded before this file.
//
//////////////////////////////////////////////////////


///////////////////////////////////////////////////
//
// A construction for Fourier type matrices. The input 
// is a cyclic group of order n. (Since this allows 
// construction of other matrices over the same group.)
// If one wishes the matrix to be human-readable, this should be 
// an FP-group.
//
////////////////////////////////////////////////////

TypeFMatrix := function(G);

	n := Order(G);
	y := G.1;
        RG := GroupAlgebra(Rationals(), G);
        x := RG!y;

	list := [(x^i)^j: i in [0..n-1], j in [0..n-1]];
        return Matrix(n,n,list);

end function;


//////////////////////////////////////////////////////////
//
// Butson's construction for a B(2p,p). All notation is precisely 
// as in Butson's original paper.
//
//////////////////////////////////////////////////////////

LeastNonResidue := function(p);

	for i in [1..p] do;
		if LegendreSymbol(i,p) eq -1 then;
		return Integers()!i;
		break i;
		end if;
	end for;
end function;


GHConjugateTranspose := function(H, G);
	Inv := map<G->G| g :-> g^-1>;
	n := NumberOfRows(H);
	RG := GroupAlgebra(Rationals(), G);
	list := [RG!Inv(H[j,i]): i in [1..n], j in [1..n]];
	return Matrix(n, n, list);
end function;
	

ButsonDouble := function(G);

	p := Order(G);
	n := LeastNonResidue(p);
	q := Integers()!((p-1)/2);

	V := TypeFMatrix(G);
	RG<x> := GroupAlgebra(Rationals(), G);

	Vlist := [(x^i)^j : i in [0..p-1], j in [0..p-1]];
	Wlist := [((x^i)^j)^n: i in [0..p-1], j in [0..p-1]];

	V := Matrix(p,p,Vlist);
	W := Matrix(p,p,Wlist);

	   
// We build Butson's Q-matrix

	Qlist := [(x^q)^(i^2): i in [0..p-1]];
	Q := DiagonalMatrix(Qlist);

	C := Q*V*Q;
	B := (Q^n)*W*(Q^n);

// We build the matrix P

	Plist := [];
	for i in [0..p-1] do;
	for j in [0..p-1] do;
		if ( j eq n*i mod p) then;
		Append(~Plist, 1);
		else 
		Append(~Plist,0);
		end if;
	end for;
	end for;

	P := Matrix(p, p, Plist);

	UL := Q*V;
	UR := (Q^n)*W;

	LL := GHConjugateTranspose(C, G);
	LR := Transpose(GHConjugateTranspose(B*P, G));

	r1 := HorizontalJoin([UL,UR]);
	r2 := HorizontalJoin([LL, LR]);
	K := VerticalJoin([r1,r2]);

	return K;
	
end function;

//////////////////////////////////////////////////////
//
// Normalisation for Generalised Hadamard matrices. 
// (Also works for Butson matrices.) GHReduceField 
// attempts to rewrite H over a smaller cyclotomic 
// field. These functions do not commute - it may be 
// necessary to apply a sequence to get to the smallest 
// possible field.
//
//////////////////////////////////////////////////////

GHNormalise := function(H);

	n := NumberOfRows(H);
	for i in [1..n] do;
		MultiplyRow(~H, H[i,1]^-1, i);
		MultiplyColumn(~H, H[1,i]^-1, i);
	end for;
	return H;
end function;

GHReduceField := function(H);

	n := NumberOfRows(H);
	Entries := [Minimise(H[i,j]): j in [1..n], i in [1..n]];

	return Matrix(n, Entries);
end function;

///////////////////////////////////////////////////////
//
// The McNulty-Wiegert construction of order nk
//
// Arguments are lists of length n, each list entry is a 
// BH of order k. Any matrix from the first list is unbiased 
// to any in the second. The template is a CHM of order n.
// 
////////////////////////////////////////////////////////

// We give a construction for MUBS and a product 
// operation (these are BHs)

MUBS := function(p);

	F := GF(p);
	A := AdditiveGroup(F);
	z := RootOfUnity(p);
	k := Parent(z);

	MUBS := [];

	for a in F do;
	List := [z^Integers()!(a*x^2+b*x): b in F, x in F];
	Append(~MUBS, Matrix(p,List));
	end for;

	return MUBS;

end function;

BHConjugateTranspose := function(H);

	n := NumberOfRows(H);
	list := [ComplexConjugate(H[i,j]): i in [1..n], j in [1..n]];
	return Matrix(n,n,list);

end function;

// Warning - it is not guaranteed that 
// the output of this function is over 
// the same roots of unity as the input. 
// Particularly if p is an odd p prime, it's possible 
// that -1s will be introduced. 
UnbiasedProduct := function(H1, H2);

//	k is the order of the smallest cyclotomic field containing 
//	the unbiased product. To normalise, we need to find a rotation 
//	of one of the entries of the product which is real. In all cases 
//	we attempted this was written in terms of roots of unity of order 
//	4*k. If an error occurs, the constant 4 can be increased.
	k := LCM(CyclotomicOrder(Parent(H1[1,1])), 
				   CyclotomicOrder(Parent(H2[1,1])));
	omega := RootOfUnity(4*k);
	n := NumberOfRows(H1);
	Mat := [];
	Base := BHConjugateTranspose(H1)*H2;
// We need to scale the product to have entries of unit norm. 
// So we find a real, positive phase and multiply through by this.
	phases := [omega^i: i in [1..4*k]];
	for i in phases do;
		x := Base[1,1]*i;
		if IsReal(x) and (Re(x) gt 0) then;
		scale := x^-1;
		break i;
		end if;
	end for;

	return scale*Base;
end function;


McNulty := function(List1,List2,Template);

	K := BHConjugateTranspose(DirectSum(List1));
	n := NumberOfRows(K);
	p := NumberOfRows(Template);
	List := [Template[i,j]*UnbiasedProduct(List1[i], List2[j])
	              : j in [1..p], i in [1..p]];
	H := GHNormalise(BlockMatrix(p, p, List));
// We minimise the field over which the matrix is written.
	
	return GHNormalise(GHReduceField(H));
end function; 



///////////////////////////////////////////////////////
//
// Finally we give explicit code for constructing an isomorphism 
// between the McNulty-Weigert matrix of order 10 and the 
// Butson matrix of the same order.
//
///////////////////////////////////////////////////////

//	G := SmallGroup(5, 1);
//	F5 := BHConjugateTranspose(GHtoBH(TypeFMatrix(G), G));
//	RG := Parent(F5[1,1]);
//	x := F5[2,2];
//
//	Diag := DiagonalMatrix(RG, [1, x, x^4, x^4, x]);
//	H3 := Diag^3*F5;
//	H4 := Diag^4*F5;
//
//	template := Matrix(Rationals(), 2, [1,1,1,-1]);	
//
//	B1 := GHNormalise(McNulty([F5^0, H3], [F5,H4], template));
//	H1 := BHtoGH(B1, G);
//
//	H2 := GHNormalise(ButsonDouble(G));
//
//	b, sigma := GHEquivalence(H1, H2, G, G);
//
//	Mats := GHPermAutToAut(sigma, G);
//
//	Mats[1]*H2*Mats[2]^-1 eq H1;



// Finally - the following code constructs a matrix equivalent 
// to the matrix displayed by McNulty and Weigert in their paper.
	
//	G := SmallGroup(7,1);
//	F7 := GHtoBH(TypeFMatrix(G),G);
//	x := F7[2,2];
//	D := DiagonalMatrix(Parent(x), [1, 1, x, x^3, x^6, x^3, x]);
//
//	template := Matrix(Rationals(), 2, [1,1,1,-1]);
//	B14 := McNulty([F7^0, D*F7], [F7, D^2*F7], template);
