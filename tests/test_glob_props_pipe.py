# stdlib
import unittest
# brainnets
from brainnets import compprops, compstats, plots
from brainnets.testscripts import mycfg


class GlobPropPipeTest(unittest.TestCase):

    def setUp(self):
        self.cfg = mycfg.get_config()
        return

    def test(self):
        compprops.compute_global_properties(self.cfg, weighted=False)
        compstats.comp_glob_prop_stats(self.cfg, False)

        compprops.compute_global_properties(self.cfg, weighted=True)
        compstats.comp_glob_prop_stats(self.cfg, True)

        compprops.comp_link_sim_matrices(self.cfg)

        plots.plot_global_uw_props(self.cfg)
        plots.plot_global_w_props(self.cfg)
