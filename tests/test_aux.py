# stdlib
import unittest
# third party
import numpy as np
# brainnets
from brainnets import aux


class AuxTest(unittest.TestCase):

    def setUp(self):
        return

    def test_get_n_largest_indices(self):
        arr = np.array([1, 2, 3, 1, 2.5])
        res = aux.get_n_largest_indices(arr, 3)
        assert list(res) == [2, 4, 1]
        self.assertRaises(AssertionError, aux.get_n_largest_indices, arr, 6)

