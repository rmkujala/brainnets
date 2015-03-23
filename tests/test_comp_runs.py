import unittest

from brainnets import compprops
from brainnets import complinkdistances as cld
from brainnets import compcoms
from brainnets.testscripts import mycfg
import numpy as np

class RunTest(unittest.TestCase):

    def setUp(self):
        self.np_set_err = np.seterr(all="raise")
        self.cfg = mycfg.get_config()

    def test_compprops_runs(self):
        cld.comp_link_distances(self.cfg)
        cld.comp_paired_common_and_differing_link_distances(
            self.cfg
        )
        cld.comp_consistent_link_distances(self.cfg)
        compprops.comp_link_sim_matrices(self.cfg)

    def tearDown(self):
        # set back to default np errors:
        np.seterr(**self.np_set_err)

