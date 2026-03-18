# Row-Wise Parametrization Using Backtracking

A python implementation of Algorithm 2 from Section 5.4.4 of the thesis is given in `backtracking.py`.
This implementation requires the input matrix be a normalized S-Hadamard $BH(n,q)$ matrices where $q$ is prime.
Note that the output parametrized matrix may differ between runs on the same input matrix as the parameter vectors are randomly ordered.

The provided script specifically works with files given in the Aalto Butson Hadamard matrix repository (see [here](https://wiki.aalto.fi/spaces/Butson/pages/120482304/Matrices+up+to+monomial+equivalence)).
We include two files taken from the Aalto repository, `BH-9-3.txt` and `BH-12-3.txt`.
Additional files can be downloaded from the linked website to try other matrices with prime order roots of unity.

The file `shad.py` is a utility for decoding matrices given in the Aalto repository.
The backtracking code uses this to read the matrix file and get the desired input matrices.

The backtracking script prints either an error message indicating that the parametrization failed (with some indication of the reason), or a parameterized matrix.
The matrix is printed with `a` denoting the parameter and `w` denoting the root of unity.

The usage of the backtracking script is shown below.
The arg `n` indicates the order of the matrix, arg `q` indicates the order of the root of unity.
Arg `all` indicates that all S-Hadamard matrices in the file named `BH-n-q.txt` should be attempted for parametrization.
Arg `i` indicates a specific index of a matrix in this file that we want to attempt.
For example, the following output is from trying all $BH(12,3)$ matrices:

```txt
Backtracking-Method % python3 backtracking.py -n 12 -q 3 -all
['w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0']
['w^0', 'a*w^0', 'w^0', 'a*w^0', 'w^1', 'a*w^1', 'w^1', 'a*w^1', 'w^2', 'a*w^2', 'w^2', 'a*w^2']
['w^0', 'a*w^0', 'w^0', 'a*w^1', 'w^0', 'w^2', 'w^2', 'a*w^2', 'w^1', 'w^1', 'w^1', 'w^2']
['w^0', 'a*w^0', 'w^1', 'a*w^2', 'w^2', 'a*w^0', 'w^1', 'a*w^2', 'w^0', 'a*w^1', 'w^2', 'a*w^1']
['w^0', 'a*w^1', 'w^0', 'a*w^2', 'w^2', 'w^1', 'w^2', 'a*w^0', 'w^2', 'w^0', 'w^1', 'w^1']
['w^0', 'a*w^1', 'w^2', 'a*w^0', 'w^1', 'w^2', 'w^0', 'a*w^2', 'w^0', 'w^2', 'w^1', 'w^1']
['w^0', 'w^1', 'w^2', 'w^1', 'w^2', 'w^0', 'w^0', 'w^1', 'w^2', 'w^1', 'w^0', 'w^2']
['w^0', 'a*w^1', 'w^2', 'a*w^2', 'w^0', 'a*w^2', 'w^1', 'a*w^1', 'w^1', 'a*w^0', 'w^2', 'a*w^0']
['w^0', 'a*w^2', 'w^1', 'a*w^0', 'w^2', 'w^0', 'w^2', 'a*w^1', 'w^1', 'w^2', 'w^1', 'w^0']
['w^0', 'a*w^2', 'w^1', 'a*w^1', 'w^0', 'w^2', 'w^1', 'a*w^0', 'w^2', 'w^2', 'w^0', 'w^1']
['w^0', 'w^2', 'w^1', 'w^2', 'w^1', 'w^1', 'w^0', 'w^2', 'w^1', 'w^0', 'w^0', 'w^2']
['w^0', 'a*w^2', 'w^2', 'a*w^1', 'w^1', 'w^1', 'w^2', 'a*w^0', 'w^0', 'w^1', 'w^2', 'w^0']

['w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0']
['a*w^0', 'a*w^0', 'w^0', 'w^0', 'a*w^1', 'w^1', 'a*w^1', 'w^1', 'a*w^2', 'a*w^2', 'w^2', 'w^2']
['w^0', 'a*w^0', 'w^0', 'w^2', 'a*w^1', 'w^2', 'a*w^2', 'w^2', 'w^0', 'w^1', 'w^1', 'w^1']
['w^0', 'a*w^0', 'w^2', 'w^1', 'a*w^2', 'w^0', 'a*w^1', 'w^2', 'w^1', 'w^0', 'w^1', 'w^2']
['a*w^0', 'a*w^1', 'w^1', 'w^2', 'a*w^0', 'w^0', 'a*w^1', 'w^2', 'a*w^2', 'a*w^2', 'w^0', 'w^1']
['w^0', 'a*w^1', 'w^2', 'w^0', 'a*w^0', 'w^2', 'a*w^2', 'w^1', 'w^1', 'w^0', 'w^2', 'w^1']
['a*w^0', 'a*w^1', 'w^2', 'w^1', 'a*w^1', 'w^2', 'a*w^0', 'w^0', 'a*w^2', 'a*w^2', 'w^1', 'w^0']
['w^0', 'a*w^1', 'w^2', 'w^2', 'a*w^2', 'w^1', 'a*w^0', 'w^1', 'w^0', 'w^1', 'w^0', 'w^2']
['w^0', 'w^2', 'w^0', 'w^1', 'w^2', 'w^1', 'w^2', 'w^0', 'w^1', 'w^2', 'w^0', 'w^1']
['w^0', 'w^2', 'w^1', 'w^0', 'w^2', 'w^0', 'w^2', 'w^1', 'w^2', 'w^1', 'w^1', 'w^0']
['w^0', 'a*w^2', 'w^1', 'w^1', 'a*w^0', 'w^2', 'a*w^1', 'w^0', 'w^0', 'w^1', 'w^2', 'w^2']
['w^0', 'a*w^2', 'w^1', 'w^2', 'a*w^1', 'w^1', 'a*w^0', 'w^2', 'w^1', 'w^0', 'w^2', 'w^0']

Backtracking-Method % 
```

The following output is from trying only the first $BH(12,3)$ matrix (which is $H_1$ which is equivalent to $H_{12}$):

```txt
Backtracking-Method % python3 backtracking.py -n 12 -q 3 -i 0
['w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0']
['a*w^0', 'a*w^0', 'w^0', 'a*w^0', 'a*w^1', 'a*w^1', 'w^1', 'a*w^1', 'a*w^2', 'a*w^2', 'a*w^2', 'w^2']
['w^0', 'a*w^0', 'w^0', 'w^1', 'a*w^0', 'a*w^2', 'w^2', 'a*w^2', 'a*w^1', 'w^1', 'a*w^1', 'w^2']
['w^0', 'w^0', 'w^1', 'w^2', 'w^2', 'w^0', 'w^1', 'w^2', 'w^0', 'w^1', 'w^2', 'w^1']
['a*w^0', 'a*w^1', 'w^0', 'a*w^2', 'a*w^2', 'a*w^1', 'w^2', 'a*w^0', 'a*w^2', 'a*w^0', 'a*w^1', 'w^1']
['a*w^0', 'a*w^1', 'w^2', 'a*w^0', 'a*w^1', 'a*w^2', 'w^0', 'a*w^2', 'a*w^0', 'a*w^2', 'a*w^1', 'w^1']
['w^0', 'a*w^1', 'w^2', 'w^1', 'a*w^2', 'a*w^0', 'w^0', 'a*w^1', 'a*w^2', 'w^1', 'a*w^0', 'w^2']
['a*w^0', 'a*w^1', 'w^2', 'a*w^2', 'a*w^0', 'a*w^2', 'w^1', 'a*w^1', 'a*w^1', 'a*w^0', 'a*w^2', 'w^0']
['a*w^0', 'a*w^2', 'w^1', 'a*w^0', 'a*w^2', 'a*w^0', 'w^2', 'a*w^1', 'a*w^1', 'a*w^2', 'a*w^1', 'w^0']
['w^0', 'w^2', 'w^1', 'w^1', 'w^0', 'w^2', 'w^1', 'w^0', 'w^2', 'w^2', 'w^0', 'w^1']
['a*w^0', 'a*w^2', 'w^1', 'a*w^2', 'a*w^1', 'a*w^1', 'w^0', 'a*w^2', 'a*w^1', 'a*w^0', 'a*w^0', 'w^2']
['w^0', 'a*w^2', 'w^2', 'w^1', 'a*w^1', 'a*w^1', 'w^2', 'a*w^0', 'a*w^0', 'w^1', 'a*w^2', 'w^0']

Backtracking-Method %
```

The following output is from trying all $BH(9,3)$ matrices:

```txt
Backtracking-Method % python3 backtracking.py -n 9 -q 3 -all 
['w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0']
['a*w^0', 'a*w^0', 'w^0', 'a*w^1', 'w^1', 'a*w^1', 'w^2', 'a*w^2', 'a*w^2']
['a*w^0', 'a*w^0', 'w^0', 'a*w^2', 'w^2', 'a*w^2', 'w^1', 'a*w^1', 'a*w^1']
['a*w^0', 'a*w^1', 'w^2', 'a*w^0', 'w^1', 'a*w^2', 'w^0', 'a*w^1', 'a*w^2']
['w^0', 'w^1', 'w^2', 'w^1', 'w^2', 'w^0', 'w^2', 'w^0', 'w^1']
['a*w^0', 'a*w^1', 'w^2', 'a*w^2', 'w^0', 'a*w^1', 'w^1', 'a*w^2', 'a*w^0']
['a*w^0', 'a*w^2', 'w^1', 'a*w^0', 'w^2', 'a*w^1', 'w^0', 'a*w^2', 'a*w^1']
['a*w^0', 'a*w^2', 'w^1', 'a*w^1', 'w^0', 'a*w^2', 'w^2', 'a*w^1', 'a*w^0']
['w^0', 'w^2', 'w^1', 'w^2', 'w^1', 'w^0', 'w^1', 'w^0', 'w^2']

['w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0']
['w^0', 'w^0', 'w^0', 'w^1', 'w^1', 'w^1', 'w^2', 'w^2', 'w^2']
['w^0', 'w^0', 'w^0', 'w^2', 'w^2', 'w^2', 'w^1', 'w^1', 'w^1']
['a*w^0', 'a*w^1', 'a*w^2', 'w^0', 'w^1', 'w^2', 'w^0', 'w^1', 'w^2']
['a*w^0', 'a*w^1', 'a*w^2', 'w^1', 'w^2', 'w^0', 'w^2', 'w^0', 'w^1']
['a*w^0', 'a*w^1', 'a*w^2', 'w^2', 'w^0', 'w^1', 'w^1', 'w^2', 'w^0']
['a*w^0', 'a*w^2', 'a*w^1', 'w^0', 'w^2', 'w^1', 'w^1', 'w^0', 'w^2']
['a*w^0', 'a*w^2', 'a*w^1', 'w^1', 'w^0', 'w^2', 'w^0', 'w^2', 'w^1']
['a*w^0', 'a*w^2', 'a*w^1', 'w^2', 'w^1', 'w^0', 'w^2', 'w^1', 'w^0']

['w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0', 'w^0']
['w^0', 'w^0', 'w^0', 'w^1', 'w^1', 'w^1', 'w^2', 'w^2', 'w^2']
['w^0', 'w^0', 'w^0', 'w^2', 'w^2', 'w^2', 'w^1', 'w^1', 'w^1']
['a*w^0', 'a*w^1', 'a*w^2', 'w^0', 'w^1', 'w^2', 'w^0', 'w^1', 'w^2']
['a*w^0', 'a*w^1', 'a*w^2', 'w^1', 'w^2', 'w^0', 'w^2', 'w^0', 'w^1']
['a*w^0', 'a*w^1', 'a*w^2', 'w^2', 'w^0', 'w^1', 'w^1', 'w^2', 'w^0']
['a*w^0', 'a*w^2', 'a*w^1', 'w^0', 'w^2', 'w^1', 'w^2', 'w^1', 'w^0']
['a*w^0', 'a*w^2', 'a*w^1', 'w^1', 'w^0', 'w^2', 'w^1', 'w^0', 'w^2']
['a*w^0', 'a*w^2', 'a*w^1', 'w^2', 'w^1', 'w^0', 'w^0', 'w^2', 'w^1']

Backtracking-Method % 
```
