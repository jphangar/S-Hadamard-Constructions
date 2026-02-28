"""
Python parser BH(n,q) matrices from Alto wiki page:
https://wiki.aalto.fi/display/Butson/Matrices+up+to+monomial+equivalence
"""

import argparse
import cmath
import sys
from pathlib import Path
from typing import List

import numpy as np


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", type=int, help="Order of matrix")
    parser.add_argument("-q", type=int, help="Order of root of unity")

    return parser.parse_args()


def decode_matrix(mat_str: str) -> np.ndarray[int]:
    """
    Decode line from file to matrix based on format given by authors on wiki page
    """

    # byte 0 -> 'B', byte 1 -> 'H', byte 2 -> version check
    version = int(mat_str[2])
    if version != 1:
        raise UnicodeDecodeError(
            reason=f"Expected version 1, but got version {version}"
        )

    # byte 3 -> 48 + rows, byte 4 -> 48 + cols
    rows = ord(mat_str[3]) - 48
    cols = ord(mat_str[4]) - 48

    # byte 5 -> 48 + number of roots of unity, q
    # roots = ord(mat_str[5]) - 48

    # print(f"Rows {rows}, Columns {cols}, Roots {roots}.")

    mat = np.zeros((rows, cols), int)

    # byte i * col + j -> 48 + [i,j]
    mat_elts_str = mat_str[6:]
    for i in range(rows):
        for j in range(cols):
            mat[i, j] = ord(mat_elts_str[i * cols + j]) - 48

    return mat


def mat_to_list(mat: np.ndarray) -> List[List[int]]:
    """Form a list of lists out of ndarray for a matrix"""
    mat_str = []
    n = mat.shape[0]
    for i in range(n):
        mat_str.append(mat[i].tolist())

    return mat_str


def herm_inner_prod_sq(r1: np.ndarray, r2: np.ndarray, q: int, n: int) -> complex:
    """Compute hermitian inner product of given n-length vectors of powers of zeta_q"""
    res = 0
    for i in range(n):
        res += (q ** (r1[i] * 2)) * (q ** (-r2[i] * 2))
    return res


def check_s_hadamard(mat: np.ndarray, q: int) -> bool:
    """True if mat is S-Hadamard, otherwise False"""
    n = mat.shape[0]
    w = cmath.exp(2 * cmath.pi * 1j / q)
    for i in range(n):
        for j in range(i + 1, n):
            if (
                herm_inner_prod_sq(mat[i, :], mat[j, :], w, n)
                > 10 ** (-9) + 10 ** (-9) * 1j
            ):
                return False
    return True


def matrix_list(n: int, q: int) -> List[List[List[int]]]:
    """Return a list of S-Hadamard matrices in the given file"""
    filepath = Path(f"BH-{n}-{q}.txt")
    mats = []

    with open(filepath, encoding="utf-8") as file:
        for line in file.readlines():
            mat = decode_matrix(line)

            if check_s_hadamard(mat, q):
                mats.append(mat_to_list(mat))

    return mats


def main() -> None:
    """Save the S-Hadamard matrices and print the number of S-Hadamard matrices found"""
    args = parse_args()

    if not args.n or not args.q:
        raise argparse.ArgumentError(
            argument=None, message="Arguments n and q are required but weren't provided"
        )

    filepath = Path(f"BH-{args.n}-{args.q}.txt")

    with open(filepath, encoding="utf-8") as file:
        # parse the file, each line is an encoded matrix
        shad_counter = 0
        mat_strs = []
        for line in file.readlines():
            mat = decode_matrix(line)
            if check_s_hadamard(mat, args.q):
                shad_counter += 1
                mat_strs.append(str(mat_to_list(mat)))

        with open(f"SBH-{args.n}-{args.q}-mpl.txt", "w", encoding="utf-8") as outfile:
            outfile.write("Evals := [\n" + ",\n".join(mat_strs) + "\n];")

        print("SHAD COUNTER: " + str(shad_counter))


if __name__ == "__main__":
    sys.exit(main())
