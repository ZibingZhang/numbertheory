from typing import List
from typing import Tuple
from typing import Set

import functools
import operator
import math
import fractions
import random


class CoprimeError(Exception):
    pass


class NotNatError(Exception):
    pass


def product(iterable: iter):
    return functools.reduce(operator.mul, iterable, 1)


def gcd(a: int, b: int) -> int:
    """
    :param a: integer
    :param b: integer
    :return: The largest non-negative integer that divides both a and b.
    """
    return math.gcd(a, b)


def lcm(a: int, b: int) -> int:
    """
    :param a: integer
    :param b: integer
    :return: The smallest non-negative integer that is divisible by both a and b.
    """
    if gcd(a, b) == 0:
        return 0

    return int(a*b/gcd(a, b))


def pyth_triple(u: int, v: int) -> Tuple[int, int, int]:
    """
    Generates a Pythagorean triple (u^2 - v^2, 2uv, u^2 + v^2).
    Conditions:
        1) u > v >= 1

    :param u: natural number
    :param v: natural number
    :return: A Pythagorean triple with a odd and b even.
    """
    if not u > v:
        raise ValueError("First argument must strictly be greater than the second.")
    elif not v >= 1:
        raise ValueError("Second argument must be greater than or equal to 1.")

    return u**2 - v**2, 2*u*v, u**2 + v**2


def prim_triple(s: int, t: int) -> Tuple[int, int, int]:
    """
    Generates a Pythagorean triple (st, (s^2 - t^2)/2, (s^2 + t^2)/2).
    Conditions:
        1) s > t >= 1
        2) gcd(s, t) == 1
        3) s and t are both odd

    :param s: natural number
    :param t: natural number
    :return: A primitive Pythagorean triple with a odd and b even.
    """
    if not s > t:
        raise ValueError("First argument must be strictly greater than the second.")
    elif not t >= 1:
        raise NotNatError("Second argument must be greater than or equal to 1.")
    elif s % 2 == 0 or t % 2 == 0:
        raise ValueError("Both arguments must be odd.")
    elif not gcd(s, t) == 1:
        raise CoprimeError("The two arguments are not coprime.")

    return s*t, int((s**2 - t**2)/2), int((s**2 + t**2)/2)


def extended_euclid(a: int, b: int) -> Tuple[int, int, int]:
    """
    Extended Euclidean Algorithm

    :param a: integer
    :param b: integer
    :return: The gcd of a and b, and x and y such that ax + by = gcd(a, b)
    """
    if a == 0 and b == 0:
        return 0, 0, 0

    equation1 = [1, 0, a]
    equation2 = [0, 1, b]

    # when we divide y by x, we say y = qx + r

    while True:
        if equation1[2] == 0:
            if equation2[2] > 0:
                return equation2[2], equation2[0], equation2[1]
            else:
                return -equation2[2], -equation2[0], -equation2[1]

        q = equation2[2] // equation1[2]
        equation2 = [equation2[index] - q * equation1[index] for index in range(3)]

        if equation2[2] == 0:
            if equation1[2] > 0:
                return equation1[2], equation1[0], equation1[1]
            else:
                return -equation1[2], -equation1[0], -equation1[1]

        q = equation1[2] // equation2[2]
        equation1 = [equation1[index] - q * equation2[index] for index in range(3)]


def prime_factors(nat: int) -> List[Tuple[int, int]]:
    """
    Every positive integer has a unique prime factorization
    (when listed in non-decreasing order).
    Conditions:
        1) nat > 0

    :param nat: natural number
    :return: An ordered list of prime factors, from least to greatest, paired with multiplicity.
    """

    def _multiplicity() -> List[Tuple[int, int]]:
        factors.sort()
        list_factors = [factors[0]]
        list_multiplicities = [1]
        for index in range(1, len(factors)):
            if factors[index] == list_factors[-1]:
                list_multiplicities[-1] += 1
            else:
                list_factors.append(factors[index])
                list_multiplicities.append(1)
        return list(zip(list_factors, list_multiplicities))

    if not nat > 0:
        raise NotNatError("Integer must be positive.")

    if nat == 1:
        return []

    factors = []
    while True:
        if nat % 2 == 0:
            nat /= 2
            factors.append(2)
        else:
            upper_bound = math.ceil(math.sqrt(nat)) + 1
            if upper_bound < 3:
                return _multiplicity()

            for factor in range(3, math.ceil(math.sqrt(nat)) + 1):
                if nat % factor == 0:
                    nat /= factor
                    factors.append(factor)
                    break
            else:
                factors.append(int(nat))
                return _multiplicity()


def lin_congruence(a: int, b: int, m: int) -> Set[int]:
    """
    Solves the linear congruence ax is congruent to b modulo m.

    ax = b (mod m)
    ax - b = my for some integer y
    ax - my = b

    Conditions:
        1) m > 0

    :param a: integer
    :param b: integer
    :param m: natural number
    :return: The solution to the linear congruence.
    """
    if not m > 0:
        raise ValueError("Modulus must be positive.")

    num_solutions = gcd(a, m)
    if b % num_solutions != 0:
        return set()
    else:
        x_naught = extended_euclid(a % m, -m)[1] * int(b / num_solutions)
        return set([(x_naught + int(k*m/num_solutions)) % m for k in range(num_solutions)])


def eulers_phi(nat: int) -> int:
    """
    The number of positive integers not exceeding nat that are relatively prime to nat.
    Conditions:
        1) nat > 0

    :param nat: natural number
    :return: phi(nat)
    """
    if not nat > 0:
        raise NotNatError("Only defined for natural numbers.")

    if nat == 1:
        return 0

    factors = [factor[0] for factor in prime_factors(nat)]
    fracs = [fractions.Fraction(p-1, p) for p in factors]
    return int(nat * product(fracs))


def inv(a: int, m: int) -> int:
    """
    Returns the inverse of a modulo m.
    Conditions
        1) gcd(a, m) == 1
        2) a != 0
        3) m > 0

    ax = 1 (mod m)
    ax - 1 = my
    ax - my = 1

    :param a: non-zero integer
    :param m: natural number
    :return: non-zero integer
    """
    if not m > 0:
        raise NotNatError("Modulus must be positive.")
    elif a == 0:
        raise ValueError("Zero has no inverse.")
    elif gcd(a, m) != 1:
        raise CoprimeError("The two arguments must be coprime.")

    return extended_euclid(a%m, -m)[1] % m


def chinese_remainder(lin_cons: List[Tuple[int, int]]) -> Tuple[int, int]:
    """
    Solves the system of linear congruences.
    The input is a list of (b_i, m_i), where x = b_i (mod m_i)
    Conditions:
        1) m's must be pairwise coprime
        2) lin_cons must be non-empty

    :param lin_cons: list of linear congruences with a = 1
    :return: Simultaneous solution to the linear congruences.
    """
    if len(lin_cons) == 0:
        raise ValueError("List of linear congruences must be non-empty.")

    for lin_con in lin_cons:
        if not lin_con[1] > 0:
            raise NotNatError("Modulus must be a natural number.")

    for index_i in range(len(lin_cons)):
        for index_j in range(index_i+1, len(lin_cons)):
            if gcd(lin_cons[index_i][1], lin_cons[index_j][1]) != 1:
                raise CoprimeError("Modulus must be pairwise coprime.")

    m_product = product([lin_con[1] for lin_con in lin_cons])
    total = 0

    for index, lin_con in enumerate(lin_cons):
        partial_m = int(m_product / lin_con[1])
        inverse = inv(partial_m, lin_con[1])
        total += lin_con[0] * partial_m * inverse

    return total % m_product, m_product


def divisors(integer: int) -> Set[int]:
    """
    Gives all positive divisors of a non-zero integer.
    https://stackoverflow.com/questions/171765/what-is-the-best-way-to-get-all-the-divisors-of-a-number

    :param integer: integer
    :return: A list of positive divisors.
    """
    if integer == 0:
        raise ValueError("Integer must be non-zero.")

    if abs(integer) == 1:
        return {1}

    list_prime_factors = prime_factors(abs(integer))
    num_unique_prime_factors = len(list_prime_factors)
    multiplicity_count = [0] * num_unique_prime_factors
    list_divisors = []
    while True:
        list_divisors.append(functools.reduce(operator.mul,
                                              [pow(list_prime_factors[x][0], multiplicity_count[x])
                                               for x in range(num_unique_prime_factors)],
                                              1))
        index = 0
        while True:
            multiplicity_count[index] += 1
            if not multiplicity_count[index] > list_prime_factors[index][1]:
                break
            multiplicity_count[index] = 0
            index += 1
            if index == num_unique_prime_factors:
                return set(list_divisors)


def sum_divisors(n: int, x: int) -> int:
    """
    Returns the sum of the positive divisors of n raised to the xth power.
    Conditions:
        1) n != 0
        2) x >= 0

    :param n: non-negative integer
    :param x: natural number
    :return: sigma_x(n)
    """
    if n == 0:
        raise ValueError("Cannot find divisors of zero.")
    elif not x >= 0:
        raise ValueError("Exponent must be non-negative.")

    if n == 1:
        return 1

    return sum([divisor**x for divisor in divisors(n)])


def mersenne_num(n: int) -> int:
    """
    Returns the nth Mersenne number.
    Conditions:
        1) n > 0

    :param n: natural number
    :return: the nth Mersenne number.
    """
    if not n > 0:
        raise NotNatError("Argument must be positive.")

    return pow(2, n) - 1


def mersenne_prime(n: int) -> int:
    """
    Returns the nth Mersenne prime.
    Conditions:
        1) n > 0

    :param n: natural number
    :return: the nth Mersenne prime
    """
    p = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
         3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243,
         110503, 132049]
    if not n > 0:
        raise NotNatError
    elif n > len(p):
        raise Exception("Larger Mersenne primes but have not been coded in...")

    return pow(2, p[n-1]) - 1


def rabin_miller(n: int, a: int = None, t: int = 1):
    """
    Rabin-Miller test for compositeness.
    Conditions:
        1) n is an odd number >= 3
        2) a is in the range [2, n-2]
        3) t >= 1

    :param n: the number being tested
    :param a: the base
    :param t: the number of times the test is being run
    :return: is the number composite? (false does not imply primality)
    """
    if not n >= 3:
        raise ValueError
    elif n % 2 == 0:
        raise ValueError
    elif not t >= 1:
        raise ValueError
    elif a is not None and not 2 <= a <= n-2:
        raise ValueError

    if n == 3:
        return False

    # n-1 = pow(2, k) * q
    k = 0
    q = n - 1
    while q % 2 == 1:
        q /= 2
        k += 1

    for trial in range(t):
        if trial != 0 or a is None:
            a = random.randint(2, n - 2)

        x = pow(a, q, n)
        if x == 1 or x == n - 1:
            continue

        for i in range(k):
            x = pow(x, 2, n)
            if x == -1:
                continue
        else:
            return True

    return False
