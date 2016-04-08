import unittest

from brainnets import plots
from brainnets.testscripts import mycfg
import matplotlib


"""
This module runs tests for plots using testdata
and saves the plots to brainnets/testdata/output
as *.pdf files.
"""


class Test(unittest.TestCase):

    def setUp(self):
        # to be able to run the script also without a window:
        # matplotlib.use('Agg')
        self.cfg = mycfg.get_config()

    def test_corr_dist_plots(self):
        assert plots.plot_pooled_correlation_dists_by_condition(self.cfg)
        assert plots.plot_corr_dists_subjectwise(self.cfg)
        assert plots.plot_link_dist_probs_by_condition(self.cfg)
        assert plots.plot_pooled_corr_t_val_dists(self.cfg)
        assert plots.plot_global_w_props(self.cfg)
        assert plots.plot_global_uw_props(self.cfg)
