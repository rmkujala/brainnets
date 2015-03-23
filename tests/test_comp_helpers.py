# stdlib
import unittest
# brainnets
from brainnets import comp_helpers


def dummy_func(x):
    return x*x

class CompHelpersTest(unittest.TestCase):

    # def setUp(self):
    #     pass

    def test_run_in_parallel(self):
        # test that the mapping is correct
        args = [1, 2, 3, 4, 5, 6, 7, 8]
        res = comp_helpers.run_in_parallel(dummy_func, args, 1)
        assert res == [1, 4, 9, 16, 25, 36, 49, 64]
        comp_helpers.run_in_parallel(dummy_func, args, 3)
        assert res == [1, 4, 9, 16, 25, 36, 49, 64]



