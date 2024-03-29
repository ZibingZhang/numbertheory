import unittest

from numbertheory import *


class Tests(unittest.TestCase):
    def test_sqrt_radical_pow(self):
        self.assertEqual(sqrt_radical_pow(3, 2, 2, 1), (3, 2))
        self.assertEqual(sqrt_radical_pow(3, 2, 2, 2), (17, 12))

        self.assertRaises(NotNatError, sqrt_radical_pow, 3, 2, 2, -1)
        self.assertRaises(NotNatError, sqrt_radical_pow, 3, 2, 2, 0)

    def test_gcd(self):
        self.assertEqual(gcd(0, 0), 0)
        self.assertEqual(gcd(5, 0), 5)
        self.assertEqual(gcd(-12, 15), 3)

    def test_lcm(self):
        self.assertEqual(lcm(0, 0), 0)
        self.assertEqual(lcm(0, 5), 0)
        self.assertEqual(lcm(6, 9), 18)

    def test_pyth_triple(self):
        triple1 = pyth_triple(2, 1)
        triple2 = pyth_triple(4, 2)
        triple3 = pyth_triple(97, 5)

        self.assertEqual(triple1, (3, 4, 5))
        self.assertEqual(triple2, (12, 16, 20))
        self.assertEqual(triple3, (9384, 970, 9434))

        self.assertEqual(triple1[0] ** 2 + triple1[1] ** 2, triple1[2] ** 2)
        self.assertEqual(triple2[0] ** 2 + triple2[1] ** 2, triple2[2] ** 2)
        self.assertEqual(triple3[0] ** 2 + triple3[1] ** 2, triple3[2] ** 2)

        self.assertRaises(ValueError, pyth_triple, 5, -1)
        self.assertRaises(ValueError, pyth_triple, 5, 0)
        self.assertRaises(ValueError, pyth_triple, 3, 4)

    def test_prim_triple(self):
        triple1 = prim_triple(3, 1)
        triple2 = prim_triple(7, 5)
        triple3 = prim_triple(97, 5)

        self.assertEqual(triple1, (3, 4, 5))
        self.assertEqual(triple2, (35, 12, 37))
        self.assertEqual(triple3, (485, 4692, 4717))

        self.assertEqual(triple1[0] ** 2 + triple1[1] ** 2, triple1[2] ** 2)
        self.assertEqual(triple2[0] ** 2 + triple2[1] ** 2, triple2[2] ** 2)
        self.assertEqual(triple3[0] ** 2 + triple3[1] ** 2, triple3[2] ** 2)

        self.assertRaises(ValueError, prim_triple, 1, 5)
        self.assertRaises(NotNatError, prim_triple, 5, 0)
        self.assertRaises(NotNatError, prim_triple, 5, -1)
        self.assertRaises(ValueError, prim_triple, 6, 3)
        self.assertRaises(ValueError, prim_triple, 7, 2)
        self.assertRaises(CoprimeError, prim_triple, 15, 5)

    def test_extended_euclid(self):
        self.assertEqual(extended_euclid(0, 0), (0, 0, 0))
        self.assertEqual(extended_euclid(1, 0), (1, 1, 0))
        self.assertEqual(extended_euclid(0, 1), (1, 0, 1))
        self.assertEqual(extended_euclid(1492, 1066), (2, -5, 7))
        self.assertEqual(extended_euclid(-1492, 1066), (2, 5, 7))
        self.assertEqual(extended_euclid(1492, -1066), (2, -5, -7))
        self.assertEqual(extended_euclid(-1492, -1066), (2, 5, -7))

    def test_prime_factors(self):
        self.assertEqual(prime_factors(1), [])
        self.assertEqual(prime_factors(2), [(2, 1)])
        self.assertEqual(prime_factors(3), [(3, 1)])
        self.assertEqual(prime_factors(7), [(7, 1)])
        self.assertEqual(prime_factors(28), [(2, 2), (7, 1)])
        self.assertEqual(prime_factors(200), [(2, 3), (5, 2)])
        self.assertEqual(prime_factors(1001), [(7, 1), (11, 1), (13, 1)])
        self.assertEqual(prime_factors(1260), [(2, 2), (3, 2), (5, 1), (7, 1)])

        self.assertRaises(NotNatError, prime_factors, -10)
        self.assertRaises(NotNatError, prime_factors, 0)

    def test_lin_congruence(self):
        self.assertEqual(lin_congruence(7, 1, 31), {9})
        self.assertEqual(lin_congruence(4, 3, 7), {6})
        self.assertEqual(lin_congruence(9, 12, 15), {3, 8, 13})
        self.assertEqual(lin_congruence(42, 12, 90), {11, 26, 41, 56, 71, 86})
        self.assertEqual(lin_congruence(0, 2, 7), set())
        self.assertEqual(lin_congruence(2, 0, 7), {0})
        self.assertEqual(lin_congruence(0, 0, 7), {0, 1, 2, 3, 4, 5, 6})
        self.assertEqual(lin_congruence(2, 2, 7), {1})

        self.assertRaises(ValueError, lin_congruence, 9, 9, -1)
        self.assertRaises(ValueError, lin_congruence, 9, 9, 0)

    def test_eulers_phi(self):
        self.assertEqual(phi(1), 0)
        self.assertEqual(phi(2), 1)
        self.assertEqual(phi(3), 2)
        self.assertEqual(phi(200), 80)

        self.assertRaises(NotNatError, phi, -1)
        self.assertRaises(NotNatError, phi, 0)

    def test_inv(self):
        self.assertEqual(inv(3, 7), 5)
        self.assertEqual(inv(5, 7), 3)
        self.assertEqual(inv(1, 26), 1)
        self.assertEqual(inv(3, 26), 9)
        self.assertEqual(inv(9, 26), 3)

        self.assertRaises(NotNatError, inv, 1, -1)
        self.assertRaises(NotNatError, inv, 1, 0)
        self.assertRaises(ValueError, inv, 0, 7)
        self.assertRaises(CoprimeError, inv, 2, 6)

    def test_chinese_remainder(self):
        self.assertEqual(chinese_remainder([(1, 5), (2, 7), (3, 9), (4, 11)]), (1731, 3465))

        self.assertRaises(CoprimeError, chinese_remainder, [(1, 2), (2, 3), (3, 4)])
        self.assertRaises(CoprimeError, chinese_remainder, [(1, 3), (2, 2), (3, 4)])
        self.assertRaises(CoprimeError, chinese_remainder, [(1, 2), (2, 4), (3, 3)])
        self.assertRaises(NotNatError, chinese_remainder, [(1, 3), (2, 5), (3, -1)])
        self.assertRaises(NotNatError, chinese_remainder, [(1, 0), (2, 5), (3, 7)])
        self.assertRaises(ValueError, chinese_remainder, [])

    def test_divisors(self):
        self.assertEqual(divisors(1), {1})
        self.assertEqual(divisors(2), {1, 2})
        self.assertEqual(divisors(3), {1, 3})
        self.assertEqual(divisors(12), {1, 2, 3, 4, 6, 12})
        self.assertEqual(divisors(146), divisors(-146))

        self.assertRaises(ValueError, divisors, 0)

    def test_sum_divisors(self):
        # results from https://en.wikipedia.org/wiki/Divisor_function
        self.assertEqual(sum_divisors(1, 0), 1)
        self.assertEqual(sum_divisors(2, 0), 2)
        self.assertEqual(sum_divisors(3, 0), 2)
        self.assertEqual(sum_divisors(4, 0), 3)
        self.assertEqual(sum_divisors(5, 0), 2)
        self.assertEqual(sum_divisors(6, 0), 4)
        self.assertEqual(sum_divisors(1, 1), 1)
        self.assertEqual(sum_divisors(2, 1), 3)
        self.assertEqual(sum_divisors(3, 1), 4)
        self.assertEqual(sum_divisors(4, 1), 7)
        self.assertEqual(sum_divisors(5, 1), 6)
        self.assertEqual(sum_divisors(6, 1), 12)
        self.assertEqual(sum_divisors(1, 2), 1)
        self.assertEqual(sum_divisors(2, 2), 5)
        self.assertEqual(sum_divisors(3, 2), 10)
        self.assertEqual(sum_divisors(4, 2), 21)
        self.assertEqual(sum_divisors(5, 2), 26)
        self.assertEqual(sum_divisors(6, 2), 50)
        self.assertEqual(sum_divisors(1, 3), 1)
        self.assertEqual(sum_divisors(2, 3), 9)
        self.assertEqual(sum_divisors(3, 3), 28)
        self.assertEqual(sum_divisors(4, 3), 73)
        self.assertEqual(sum_divisors(5, 3), 126)
        self.assertEqual(sum_divisors(6, 3), 252)

        self.assertEqual(sum_divisors(125, 1), 156)
        self.assertEqual(sum_divisors(-125, 1), 156)

        self.assertRaises(ValueError, sum_divisors, 0, 5)
        self.assertRaises(ValueError, sum_divisors, 5, -5)

    def test_mersenne_num(self):
        self.assertEqual(mersenne_num(1), 1)
        self.assertEqual(mersenne_num(2), 3)
        self.assertEqual(mersenne_num(3), 7)
        self.assertEqual(mersenne_num(4), 15)
        self.assertEqual(mersenne_num(10), 1023)

        self.assertRaises(NotNatError, mersenne_num, -1)
        self.assertRaises(NotNatError, mersenne_num, 0)

    def test_mersenne_prime(self):
        self.assertEqual(mersenne_prime(9), 2305843009213693951)

        self.assertRaises(NotNatError, mersenne_prime, -1)
        self.assertRaises(NotNatError, mersenne_prime, 0)

    def test_rabin_miller(self):
        self.assertRaises(ValueError, rabin_miller, 1, 2)
        self.assertRaises(ValueError, rabin_miller, 6, 2)
        self.assertRaises(ValueError, rabin_miller, 11, 10)
        self.assertRaises(ValueError, rabin_miller, 11, 5, -1)
        self.assertRaises(ValueError, rabin_miller, 11, 5, 0)

        self.assertEqual(rabin_miller(3), False)
        self.assertEqual(rabin_miller(5), False)
        self.assertEqual(rabin_miller(7), False)
        self.assertEqual(rabin_miller(7, t=10), False)
        self.assertEqual(rabin_miller(17, t=10), False)
        self.assertEqual(rabin_miller(pow(2, 31) - 1, 10), False)
        self.assertEqual(rabin_miller(pow(2, 31) - 1, 100, 10), False)

        self.assertEqual(rabin_miller(33, 5), True)
        self.assertEqual(rabin_miller(33, 6), True)
        self.assertEqual(rabin_miller(33, 7), True)
        self.assertEqual(rabin_miller(33, 8), True)

        # even though 221 is composite, this test returns false
        self.assertEqual(rabin_miller(221, 174), False)

    def test_order(self):
        self.assertEqual(order(5, 13), 4)
        self.assertEqual(order(3, 7), 6)
        self.assertEqual(order(-4, 7), 6)
        self.assertEqual(order(10, 7), 6)
        self.assertEqual(order(-1, 7), 2)
        self.assertEqual(order(1, 7), 1)

        self.assertRaises(CoprimeError, order, 3, 6)
        self.assertRaises(NotNatError, order, 3, -1)
        self.assertRaises(CoprimeError, order, 3, 0)
        self.assertRaises(CoprimeError, order, 0, 9)
        self.assertRaises(ValueError, order, 0, 1)

    def test_trial_division(self):
        self.assertTrue(not trial_division(1))
        self.assertTrue(trial_division(2))
        self.assertTrue(trial_division(3))
        self.assertTrue(not trial_division(4))
        self.assertTrue(trial_division(5))
        self.assertTrue(not trial_division(9409))
        self.assertTrue(trial_division(997))

        self.assertRaises(NotNatError, trial_division, -1)
        self.assertRaises(NotNatError, trial_division, 0)

    def test_jacobi_sym(self):
        self.assertEqual(jacobi(1, 1), 1)
        self.assertEqual(jacobi(1, 3), 1)
        self.assertEqual(jacobi(1, 9), 1)
        self.assertEqual(jacobi(3, 1), 1)
        self.assertEqual(jacobi(9, 1), 1)
        self.assertEqual(jacobi(12, 21), 0)
        self.assertEqual(jacobi(17, 21), 1)
        self.assertEqual(jacobi(21, 17), 1)
        self.assertEqual(jacobi(21, 45), 0)
        self.assertEqual(jacobi(7, 19), 1)
        self.assertEqual(jacobi(19, 23), -1)
        self.assertEqual(jacobi(20, 23), -1)
        self.assertEqual(jacobi(21, 23), -1)
        self.assertEqual(jacobi(22, 23), -1)
        self.assertEqual(jacobi(45, 23), -1)
        self.assertEqual(jacobi(-1, 23), -1)
        self.assertEqual(jacobi(-24, 23), -1)

        self.assertRaises(ValueError, jacobi, 0, 5)
        self.assertRaises(ValueError, jacobi, 1, 6)
        self.assertRaises(NotNatError, jacobi, 3, -1)
        self.assertRaises(NotNatError, jacobi, 3, 0)

    def test_legendre_sym(self):
        self.assertEqual(legendre(9, 3), 0)
        self.assertEqual(legendre(2, 3), -1)
        self.assertEqual(legendre(1, 3), 1)
        self.assertEqual(legendre(30, 127), 1)
        self.assertEqual(legendre(157, 127), 1)
        self.assertEqual(legendre(-7, 127), 1)
        self.assertEqual(legendre(18, 83), -1)

        self.assertRaises(ValueError, legendre, 0, 5)
        self.assertRaises(NotNatError, legendre, 1, -5)

    def test_two_squares(self):
        self.assertEqual(two_squares(1), (0, 1))
        self.assertEqual(two_squares(2), (1, 1))
        self.assertEqual(two_squares(3), (-1, -1))
        self.assertEqual(two_squares(4), (0, 2))
        self.assertEqual(two_squares(5), (1, 2))
        self.assertEqual(two_squares(6), (-1, -1))
        self.assertEqual(two_squares(349), (5, 18))
        self.assertEqual(two_squares(18946512), (2004, 3864))

    def test_square_triangular(self):
        self.assertEqual(square_triangular(0), (0, 0, 0))
        self.assertEqual(square_triangular(1), (1, 1, 1))
        self.assertEqual(square_triangular(2), (36, 8, 6))
        self.assertEqual(square_triangular(3), (1225, 49, 35))

        self.assertRaises(NotNatError, square_triangular, -1)

    def test_mobius_func(self):
        self.assertEqual(mobius(1), 1)
        self.assertEqual(mobius(2), -1)
        self.assertEqual(mobius(3), -1)
        self.assertEqual(mobius(4), 0)
        self.assertEqual(mobius(5), -1)
        self.assertEqual(mobius(6), 1)
        self.assertEqual(mobius(7), -1)
        self.assertEqual(mobius(8), 0)
        self.assertEqual(mobius(9), 0)
        self.assertEqual(mobius(10), 1)

        self.assertRaises(NotNatError, mobius, -1)
        self.assertRaises(NotNatError, mobius, 0)

    def test_figurate_numbers(self):
        self.assertEqual(triangular(0), 0)
        self.assertEqual(triangular(1), 1)
        self.assertEqual(triangular(2), 3)
        self.assertEqual(triangular(9), 45)
        self.assertEqual(square(0), 0)
        self.assertEqual(square(1), 1)
        self.assertEqual(square(2), 4)
        self.assertEqual(square(9), 81)
        self.assertEqual(pentagonal(0), 0)
        self.assertEqual(pentagonal(1), 1)
        self.assertEqual(pentagonal(2), 5)
        self.assertEqual(pentagonal(9), 117)
        self.assertEqual(hexagonal(0), 0)
        self.assertEqual(hexagonal(1), 1)
        self.assertEqual(hexagonal(2), 6)
        self.assertEqual(hexagonal(9), 153)
        self.assertEqual(heptagonal(0), 0)
        self.assertEqual(heptagonal(1), 1)
        self.assertEqual(heptagonal(2), 7)
        self.assertEqual(heptagonal(9), 189)
        self.assertEqual(octagonal(0), 0)
        self.assertEqual(octagonal(1), 1)
        self.assertEqual(octagonal(2), 8)
        self.assertEqual(octagonal(9), 225)
        self.assertEqual(nonagonal(0), 0)
        self.assertEqual(nonagonal(1), 1)
        self.assertEqual(nonagonal(2), 9)
        self.assertEqual(nonagonal(9), 261)
        self.assertEqual(decagonal(0), 0)
        self.assertEqual(decagonal(1), 1)
        self.assertEqual(decagonal(2), 10)
        self.assertEqual(decagonal(9), 297)

        self.assertRaises(ValueError, triangular, -1)
        self.assertRaises(ValueError, square, -1)
        self.assertRaises(ValueError, pentagonal, -1)
        self.assertRaises(ValueError, hexagonal, -1)
        self.assertRaises(ValueError, heptagonal, -1)
        self.assertRaises(ValueError, octagonal, -1)
        self.assertRaises(ValueError, nonagonal, -1)
        self.assertRaises(ValueError, decagonal, -1)

    def test_sqrt_frac(self):
        self.assertEqual(sqrt_frac(0), [0])
        self.assertEqual(sqrt_frac(1), [1])
        self.assertEqual(sqrt_frac(2), [1, 2])
        self.assertEqual(sqrt_frac(3), [1, 1, 2])
        self.assertEqual(sqrt_frac(4), [2])
        self.assertEqual(sqrt_frac(5), [2, 4])
        self.assertEqual(sqrt_frac(6), [2, 2, 4])
        self.assertEqual(sqrt_frac(7), [2, 1, 1, 1, 4])
        self.assertEqual(sqrt_frac(8), [2, 1, 4])
        self.assertEqual(sqrt_frac(9), [3])
        self.assertEqual(sqrt_frac(10), [3, 6])
        self.assertEqual(sqrt_frac(11), [3, 3, 6])
        self.assertEqual(sqrt_frac(12), [3, 2, 6])
        self.assertEqual(sqrt_frac(13), [3, 1, 1, 1, 1, 6])
        self.assertEqual(sqrt_frac(14), [3, 1, 2, 1, 6])
        self.assertEqual(sqrt_frac(15), [3, 1, 6])
        self.assertEqual(sqrt_frac(16), [4])
        self.assertEqual(sqrt_frac(17), [4, 8])
        self.assertEqual(sqrt_frac(18), [4, 4, 8])
        self.assertEqual(sqrt_frac(19), [4, 2, 1, 3, 1, 2, 8])
        self.assertEqual(sqrt_frac(20), [4, 2, 8])
        self.assertEqual(sqrt_frac(22), [4, 1, 2, 4, 2, 1, 8])
        self.assertEqual(sqrt_frac(31), [5, 1, 1, 3, 5, 3, 1, 1, 10])

    def test_fundamental_pell(self):
        self.assertEqual(fundamental_pell(2), (3, 2))
        self.assertEqual(fundamental_pell(3), (2, 1))
        self.assertEqual(fundamental_pell(6), (5, 2))
        self.assertEqual(fundamental_pell(10), (19, 6))
        self.assertEqual(fundamental_pell(17), (33, 8))
        self.assertEqual(fundamental_pell(26), (51, 10))

        self.assertEqual(fundamental_pell(22), (197, 42))
        self.assertEqual(fundamental_pell(31), (1520, 273))
        self.assertEqual(fundamental_pell(58), (19603, 2574))
        self.assertEqual(fundamental_pell(73), (2281249, 267000))

        self.assertEqual(fundamental_pell(97), (62809633, 6377352))
        self.assertEqual(fundamental_pell(109), (158070671986249, 15140424455100))

    def test_ind(self):
        self.assertEqual(ind(3, 13, 17), 4)
        self.assertEqual(ind(3, 3, 17), 1)
        self.assertEqual(ind(20, 3, 17), 1)
        self.assertEqual(ind(3, 20, 17), 1)


if __name__ == "__main__":
    unittest.main()
