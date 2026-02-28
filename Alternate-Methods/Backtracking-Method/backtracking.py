"""
Backtracking algorithm for parameterizing a BH(n,q) matrix.
We only support parameterizing a normalized matrix by a single parameter
where the parameter may only take on powers 0 or 1.
Also, we assume q is prime.
"""

import argparse
import itertools
import random
import sys
from typing import List, Tuple

import shad


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", type=int, help="Order of matrix")
    parser.add_argument("-q", type=int, help="Order of root of unity")
    parser.add_argument("-i", type=int, help="Index of matrix in listing")
    parser.add_argument("-all", action="store_true", help="Parameterize all indices")

    return parser.parse_args()


def generate_words(n: int, q: int) -> List[Tuple[int]]:
    """
    Generate all binary words of length n with the number of occurrences of 1
    divisible by q but not equal to n. The equality check prevents the word of
    length n with all 1's from being generated.
    """
    words = []
    for word in itertools.product([0, 1], repeat=n):
        total_sum = sum(word)
        if total_sum != n and total_sum % q == 0:
            words.append(word)
    return words


def is_orthogonal(
    r1: List[Tuple[int, int]], r2: List[Tuple[int, int]], n: int, q: int
) -> bool:
    """
    Returns True if the given rows are orthogonal, otherwise False.
    The rows are given as a list of pairs of integers where the first element
    of the pair is the power of zeta_q and the second element is the power of
    the parameter. Note that we don't check if the rows remain orthogonal if
    each entry is squared because q is assumed to be prime so this will follow
    since f(k) = 2k mod q is a bijection for 0 <= k <= q-1.
    """

    # We compute the multiset of differences for these two rows to check orthogonality
    # If the resulting set is balanced, then the rows are orthogonal, otherwise not

    # The power of the parameter in a difference is 1, 0, or -1
    # We track for each difference of parameter powers, a multiset of the zeta_q powers
    # with a factor of the given parameter power
    diffs = {}

    # Compute the differences of the rows, update diffs as we progress
    for i in range(n):
        zeta_power = (r1[i][0] - r2[i][0]) % q  # zeta_q power of current difference
        param_power = r1[i][1] - r2[i][1]  # parameter power of current difference
        if not diffs.get(param_power):
            diffs[param_power] = {}  # init param power set

        # get the multiplicity of the power of zeta_q for parameter power
        zeta_multiplicity = diffs[param_power].get(zeta_power, 0)
        # increment the multiplicity of this power of zeta_q for this parameter power
        diffs[param_power][zeta_power] = zeta_multiplicity + 1

        # Since q is prime, the multiplicity of a single power of zeta_q cannot exceed
        # the index (n/q). If it does, then these rows are not orthogonal
        if (zeta_multiplicity + 1) > (n // q):
            return False

    # Check if differences are balanced for each parameter power. For the differences to
    # be balanced, for each parameter power that was found, each power of zeta_q must be in
    # the multiset and the multiplicity of each power of zeta_q must be the same.
    # This is equivalent to the orthogonality check because q is prime.
    for k in set(diffs.keys()):
        if len(diffs[k]) != q:
            # missing powers of zeta_q in the multiset
            return False
        if len(set(diffs[k].values())) != 1:
            # multiplicity of powers of zeta_q is not consistent
            return False

    # Differences are balanced, so rows are orthogonal
    return True


def backtracking_parameterization(
    R: List[List[Tuple[int, int]]],
    L: List[List[int]],
    n: int,
    q: int,
    words: List[List[int]],
) -> List[List[Tuple[int, int]]]:
    """
    Compute the parameterization of row r with the existing rows R.
    If r is compatible, compute the parameterization for the next row in L.
    Otherwise, the parameterization failed and we return an error.
    If there are no rows remaining in L, then we return the resulting
    matrix from appending r to R.

    Args:
        R (List[List[Tuple[int, int]]]): Rows that have already been parameterized, tuple is
            (zeta_q power, param power). The param power is binary, 0 or 1.
        L (List[List[int]]): Remaining rows to parameterize
        n (int): order of matrix
        q (int): order of roots of unity, zeta_q
        words (List[List[int]]): binary vectors of parameter choices for rows

    Returns:
        (List[List[Tuple[int, int]]]): parametrized n x n matrix formed by parametrizing
            remaining rows in L with R such that all pairs of rows are orthogonal

    Raises:
        Exception: No binary parameter vector led to a valid parametrization of the next
            remaining row with the parametrized rows in R
    """

    # Base case - no rows remaining to parametrize
    if len(L) == 0:
        return R

    # Recursive case - parametrize the first row in L such that it is orthogonal to rows in R
    new_r = L[0]  # next row to parametrize
    new_L = L[1:]  # update L by removing row we are going to parametrize

    # Iterate through parameter vectors to find a valid parametrization
    for word in words:
        # Form parametrized row with current parameter vector choice
        new_rp = [(new_r[k], word[k]) for k in range(n)]

        # Check if new_rp is orthogonal with each row in R
        orthogonal_check = True
        for row in R:
            if not is_orthogonal(new_rp, row, n, q):
                # Fail check the moment we find a row in R that is not orthogonal with new_rp
                orthogonal_check = False
                break

        # Try to parameterize next row in L with the new_rp added to R if we passed the check
        try:
            if orthogonal_check:
                new_R = R + [new_rp]
                return backtracking_parameterization(new_R, new_L, n, q, words)
        except ValueError:
            # catch exception from recursive call failing because no parameter vector worked
            # for next row in L with parametrization given be new_rp
            # This means that we should try the next parameter vector choice for new_r
            continue

    # Raise an exception if we tried all parameter vectors and did not find one that gave
    # parametrization of new_r that led to a parametrization of the remaining rows that is
    # compatible with the existing parametrized rows
    raise ValueError("No word led to valid parameterization for this row")


def check_parameterization_nontrivial(M: List[List[Tuple[int, int]]]) -> bool:
    """
    Returns True if the parametrized matrix returned by backtracking_parametrization contains
    parameters. Otherwise, returns False.
    """
    for row in M:
        for entry in row:
            if entry[1] != 0:
                return True
    return False


def pretty_print_matrix(M: List[List[Tuple[int, int]]], n: int) -> None:
    """
    Print matrix in a readable format. Use "a" to denote the parameter, and
    "w" to denote the root of unity.
    """
    for i in range(n):
        row = ["" for j in range(n)]
        for j in range(n):
            row[j] = f"a*w^{M[i][j][0]}" if M[i][j][1] == 1 else f"w^{M[i][j][0]}"
        print(row)
    print("")


def run_parameterization(
    H: List[List[int]], n: int, q: int, words: List[List[int]]
) -> None:
    """
    Run parametrization function, the result will be printed to the console. We init
    the parametrization call with the first row of the given matrix. This row will
    not be parametrized.
    """
    r = [(H[0][i], 0) for i in range(n)]
    try:
        M = backtracking_parameterization([r], H[1:], n, q, words)
    except ValueError:
        print("FAIL: Could not find parametrization for given input")
    if not check_parameterization_nontrivial(M):
        # Rerun to get different word ordering, or unsuffle and reverse words order for small orders
        print("FAIL: Trivial parametrization, try again!")
    else:
        pretty_print_matrix(M, n)  # parametrization worked and was nontrivial


def main() -> None:
    args = parse_args()

    # Using S-Hadamard matrices from specified Aalto wiki BH mats file
    mats = shad.matrix_list(args.n, args.q)

    # Generating words ...
    words = generate_words(args.n, args.q)
    triv1 = words.pop(0)  # Remove all zeros vector and save for later
    words.pop()  # Throw out vector of all ones - last element in list of words
    random.shuffle(words)  # Randomize the list of words
    words.append(triv1)  # Add the vector of zeros back to the end of the list

    # Run parametrization algorithm for each matrix selected
    try:
        if args.all:
            for H in mats:
                run_parameterization(H, args.n, args.q, words)
        else:
            H = mats[args.i]
            run_parameterization(H, args.n, args.q, words)
    except Exception as err:
        print("FAIL: Unhandled exception, see exception stacktrace")
        print(err)


if __name__ == "__main__":
    sys.exit(main())
