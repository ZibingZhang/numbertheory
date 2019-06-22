from typing import List
from typing import Tuple
from typing import Set

import functools
import operator
import math
import fractions
import random

# ======================================================================
# Errors
# ======================================================================


class CoprimeError(Exception):
    """
    Integers are not relatively prime.
    """
    pass


class NotNatError(Exception):
    """
    Integer is not a natural number.
    """
    pass


class NotPrimeError(Exception):
    """
    Integer is not prime.
    """
    pass


class UnreachableError(Exception):
    """
    Reached code that in theory is unreachable.
    """
    pass


# ======================================================================
# Utils
# ======================================================================


def product(iterable: iter):
    return functools.reduce(operator.mul, iterable, 1)


def sqrt_radical_pow(a: int, b: int, c: int, n: int) -> Tuple[int, int]:
    """
    (x + y*sqrt(c)) = (a + b*sqrt(c))^n
                    = (u + v*sqrt(c)) * (a + b*sqrt(c))

    :param a: integer
    :param b: integer
    :param c: integer
    :param n: natural number
    :return: (x, y)
    """
    if not n > 0:
        raise NotNatError

    x = a
    y = b
    for i in range(n-1):
        old_x = x
        old_y = y

        x = old_x*a + old_y*b*c
        y = old_x*b + old_y*a

    return x, y


# ======================================================================
# Functions
# ======================================================================


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
    :param m: modulus
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


def phi(nat: int) -> int:
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
    :param m: modulus
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
    :param x: exponent
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


def rabin_miller(n: int, a: int = 2, t: int = 1) -> bool:
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
    elif not 2 <= a <= n-2 and n != 3:
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


def order(a: int, m: int) -> int:
    """
    Multiplicative order of a modulo m.
    Conditions:
        1) a != 0
        2) m > 0
        3) gcd(a, m) == 1

    :param a: the base
    :param m: the modulus
    :return: Returns the order of a.
    """
    if not gcd(a, m) == 1:
        raise CoprimeError("The two parameters must be relatively prime.")
    elif not m > 0:
        raise NotNatError("Modulus must be positive.")
    elif a == 0:
        raise ValueError("a must be non-negative.")

    a %= m
    if a == 1:
        return 1

    a_k = a
    for k in range(2, phi(m) + 1):
        a_k = (a_k * a) % m
        if a_k == 1:
            return k
    else:
        raise UnreachableError


def trial_division(n: int) -> bool:
    """
    Tests the given natural number for primality through trial division.

    :param n: natural number
    :return: is the natural number a prime?
    """
    if not n > 0:
        raise NotNatError

    if n == 1:
        return False
    if n == 2:
        return True

    for i in range(2, int(math.sqrt(n) + 1)):
        if n % i == 0:
            return False
    else:
        return True


def jacobi(a: int, m: int) -> int:
    """
    The Jacobi symbol.
    Conditions:
        1) m is an odd positive integer
        2) a is a non-zero integer

    :param a: integer
    :param m: natural number
    :return: (a/m)
    """
    if a == 0:
        raise ValueError
    elif not m > 0:
        raise NotNatError
    elif m % 2 == 0:
        raise ValueError

    if gcd(a, m) != 1:
        return 0
    if a == 1 or m == 1:
        return 1

    multiplier = 1
    while True:
        a %= m

        while a % 2 == 0:
            a = a // 2
            if m % 8 == 3 or m % 8 == 5:
                multiplier *= -1

        if a == 1:
            return multiplier*a
        if a == -1:
            return multiplier*a * (1 if m % 4 == 1 else -1)

        if a % 4 == 3 and m % 4 == 3:
            multiplier *= -1
        temp = a
        a = m
        m = temp


def legendre(a: int, p: int) -> int:
    """
    The Legendre symbol.
    Conditions:
        1) p is an odd positive prime
        2) a is an odd prime

    :param a: integer
    :param p: prime
    :return: (a/p)
    """
    if not trial_division(p):
        raise NotPrimeError

    return jacobi(a, p)


def two_squares(n: int) -> Tuple[int, int]:
    """
    Returns two non-negative integers such that the sum of their squares is
    equal to the given natural number. If such a thing is not possible, then
    returns (-1, -1).
    Conditions:
        1) n > 0

    :param n: natural number
    :return: (int, int)
    """
    if not n > 0:
        raise NotNatError

    if n == 1:
        return 0, 1

    factors = prime_factors(n)

    for factor in factors:
        if factor[0] % 4 == 3 and factor[1] % 2 == 1:
            return -1, -1

    pairs = []
    m_squared = 1

    for factor in factors:
        if factor[0] % 4 == 1:
            p = factor[0]
            x = y = 1
            for i in range(int(math.sqrt(p)), p):
                if (i ** 2) % p == p - 1:
                    x = i
                    break

            while True:
                if (x**2 + y**2) % p != 0:
                    raise Exception("Shouldn't have happened.")

                m = (x**2 + y**2) // p
                if m == 1:
                    pairs.append((abs(x), abs(y)))
                    break

                # not sure if flooring is okay
                r = x % m
                if r > m/2:
                    r -= m
                s = y % m
                if s > m/2:
                    s -= m

                a = (r*x + s*y) // m
                b = (r*y - s*x) // m

                x = a
                y = b

        elif factor[0] % 4 == 3:
            m_squared *= int(math.sqrt(pow(factor[0], factor[1])))
        elif factor[0] == 2:
            for i in range(factor[1]):
                pairs.append((1, 1))
        else:
            raise UnreachableError

    current_pair = pairs[0]
    for pair in pairs[1:]:
        current_pair = (abs(current_pair[0] * pair[0] + current_pair[1] * pair[1]),
                        abs(current_pair[0] * pair[1] - current_pair[1] * pair[0]))

    current_pair = list(current_pair)
    current_pair.sort()
    return current_pair[0] * m_squared, current_pair[1] * m_squared


def square_triangular(n: int) -> Tuple[int, int, int]:
    """
    Returns the nth square triangular number, as well as the index of the
    triangular number and square number respectively.
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: a three element tuple of natural numbers
    """
    if not n >= 0:
        raise NotNatError

    if n == 0:
        return 0, 0, 0

    indices = sqrt_radical_pow(3, 2, 2, n)
    triangle_index = (indices[0]-1)//2
    triangle = (triangle_index * (triangle_index+1))//2
    square_index = indices[1]//2
    square = square_index**2
    assert triangle == square
    return triangle, triangle_index, square_index


def mobius(n: int) -> int:
    """
    Mobius Function
    Conditions:
        1) n > 0

    :param n: natural number
    :return: mu(n)
    """
    if not n > 0:
        raise NotNatError

    factors = prime_factors(n)
    for factor in factors:
        if factor[1] > 1:
            return 0
    else:
        return (-1)**len(factors)


def triangular(n: int) -> int:
    """
    Triangular Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth triangular number
    """
    if not n >= 0:
        raise ValueError

    return n*(n+1)//2


def square(n: int) -> int:
    """
    Square Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth square number
    """
    if not n >= 0:
        raise ValueError

    return n**2


def pentagonal(n: int) -> int:
    """
    Pentagonal Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth pentagonal number
    """
    if not n >= 0:
        raise ValueError

    return (3*n**2 - n)//2


def hexagonal(n: int) -> int:
    """
    Hexagonal Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth hexagonal number
    """
    if not n >= 0:
        raise ValueError

    return 2*n**2 - n


def heptagonal(n: int) -> int:
    """
    Heptagonal Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth heptagonal number
    """
    if not n >= 0:
        raise ValueError

    return n*(5*n - 3)//2


def octagonal(n: int) -> int:
    """
    Octagonal Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth octagonal number
    """
    if not n >= 0:
        raise ValueError

    return 3*n**2 - 2*n


def nonagonal(n: int) -> int:
    """
    Nonagonal Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth nonagonal number
    """
    if not n >= 0:
        raise ValueError

    return n*(7*n - 5)//2


def decagonal(n: int) -> int:
    """
    Decagonal Number
    Conditions:
        1) n >= 0

    :param n: non-negative integer
    :return: nth decagonal number
    """
    if not n >= 0:
        raise ValueError

    return 4*n**2 - 3*n
