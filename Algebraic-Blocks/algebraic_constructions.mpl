kernelopts(printbytes=false):

with(LinearAlgebra):

read "../Block-Structure/utils.mpl":
read "../Block-Structure/H12.m":

# These are the simplified 3-parameter solutions we obtain by solving algebraically for block-structured parametrizations of H12

N1 := Matrix([[1,1,1,1], [1,1,c,b], [1,1/b,a,1], [1,1/c,1,1/a]]):
N2 := Matrix([[1,1,1,1], [1,1,c,b], [1,c,c,b*c/a], [1,b,a,b]]):

P1 := [
    Matrix([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]),
    Matrix([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]])
]:
P2 := [
    Matrix([[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
    Matrix([[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]])
]:

Mats := [
    H12*~generateParameterMatrix(P1[1], N1), H12*~generateParameterMatrix(P1[2], N1),
    H12*~generateParameterMatrix(P2[1], N2), H12*~generateParameterMatrix(P2[2], N2)
]:
