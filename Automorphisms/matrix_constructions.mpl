kernelopts(printbytes=false):

with(LinearAlgebra):

read "../Block-Structure/utils.mpl":
read "../Block-Structure/H12.m":

# These are the 3-parameter parametrized matrices after applying the parameter relabeling and
# scaling by parameters. These are the forms to use for determining equivalence.

N1 := Matrix([[1,1,1,1], [1,1,c,b], [1,1/b,a,1], [1,1/c,1,1/a]]):
N2 := Matrix([[1,1,1,1], [1,1,c,b], [1/c,1,1,1/a], [1/b,1,a,1]]):

P1 := [
    Matrix(4, 4, [[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]),
    Matrix(4, 4, [[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]])
]:
P2 := [
    Matrix(4, 4, [[0, 1, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1]]),
    Matrix(4, 4, [[0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0]])
]:

Mats := [
    H12*~generateParameterMatrix(P1[1], N1), H12*~generateParameterMatrix(P1[2], N1),
    H12*~generateParameterMatrix(P2[1], N2), H12*~generateParameterMatrix(P2[2], N2)
]: