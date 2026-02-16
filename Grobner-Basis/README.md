# Grobner Basis Approach to Parametrizing Matrices

This directory contains the code for the example used to illustrate how a numeric S-Hadamard matrix can be parametrized using Grobner bases.
This example is referred to in Section 3.2.2 and is done in Maple.

`H.m` contains a Maple matrix `H` which is a 12 x 12 S-Hadamard matrix found by our optimizer; this is the matrix we will be parametrizing.
This matrix is also the first matrix given in `../PSLQ/order-12-optimizer-solns.m`; that is, the first matrix in the list of order 12 S-Hadamard matrices returned by our optimizer.

The Maple worksheet called `GrobnerBasisMethod.mw` details the steps of the method, please refer to this worksheet to reproduce the computations.
A pdf of the worksheet is also available (`GrobnerBasisMethod.pdf`).

The resulting parametrized matrix is stored in `S.m` as a Maple matrix named `S`.
Matrix `S` is a 12 x 12 symbolic S-Hadamard matrix with three unimodular complex parameters `a`, `b`, and `c`, and `w` is used to denote the cube root of unity.
