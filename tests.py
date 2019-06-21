import unittest

from numtheory import *


class Tests(unittest.TestCase):
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
        self.assertRaises(ValueError, prim_triple, 5, 0)
        self.assertRaises(ValueError, prim_triple, 5, -1)
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

        self.assertRaises(ValueError, prime_factors, -10)
        self.assertRaises(ValueError, prime_factors, 0)

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

    def test_euler_phi(self):
        self.assertEqual(euler_phi(1), 0)
        self.assertEqual(euler_phi(2), 1)
        self.assertEqual(euler_phi(3), 2)
        self.assertEqual(euler_phi(200), 80)

        self.assertRaises(ValueError, euler_phi, -1)
        self.assertRaises(ValueError, euler_phi, 0)

    def test_inv(self):
        self.assertEqual(inv(3, 7), 5)
        self.assertEqual(inv(5, 7), 3)
        self.assertEqual(inv(1, 26), 1)
        self.assertEqual(inv(3, 26), 9)
        self.assertEqual(inv(9, 26), 3)

        self.assertRaises(ValueError, inv, 1, -1)
        self.assertRaises(ValueError, inv, 1, 0)
        self.assertRaises(ValueError, inv, 0, 7)
        self.assertRaises(CoprimeError, inv, 2, 6)

    def test_chinese_remain(self):
        self.assertEqual(chinese_remain([(1, 5), (2, 7), (3, 9), (4, 11)]), (1731, 3465))

        self.assertRaises(ValueError, chinese_remain, [(1, 2), (2, 3), (3, 4)])
        self.assertRaises(ValueError, chinese_remain, [(1, 3), (2, 2), (3, 4)])
        self.assertRaises(ValueError, chinese_remain, [(1, 2), (2, 4), (3, 3)])
        self.assertRaises(ValueError, chinese_remain, [(1, 3), (2, 5), (3, -1)])
        self.assertRaises(ValueError, chinese_remain, [(1, 0), (2, 5), (3, 7)])

    def test_divisor_sum(self):
        self.assertEqual(divisor_sum(1, 0), 1)
        self.assertEqual(divisor_sum(2, 0), 2)
        self.assertEqual(divisor_sum(3, 0), 2)
        self.assertEqual(divisor_sum(4, 0), 3)


if __name__ == "__main__":
    unittest.main()
