# stdlib
import os
import unittest
from copy import copy
# third party
import numpy as np
# brainnets
from brainnets.testscripts import mycfg
from brainnets import dataio


class DataioTest(unittest.TestCase):

    def setUp(self):
        self.cfg = mycfg.get_config()

    def test_blacklist_index_fname_to_blacklist_bool_list(self):
        ind_bl_fname = copy(self.cfg['blacklist_fname']).replace("bool", "ind")
        out_bl_fname = self.cfg['blacklist_fname'].replace("bool", "bool_test")
        assert "all_fnames" in self.cfg
        assert self.cfg["all_fnames"]
        dataio.blacklist_index_fname_to_blacklist_bool_list(
            ind_bl_fname,
            out_bl_fname,
            self.cfg["all_fnames"][0]
        )
        bl = dataio.get_ok_nodes(out_bl_fname)
        assert type(bl) is np.ndarray
        assert bl.dtype == bool
        # specific for the blacklist provided in the testdata:
        assert not bl[0]
        assert not bl[-1]
        os.remove(out_bl_fname)
        assert not os.path.exists(out_bl_fname)

    def test_index_blacklist_to_bool_blacklist(self):
        ind_bl = [1, 3]
        tot_n_nodes = 5
        bl = dataio.index_blacklist_to_bool_blacklist(ind_bl, tot_n_nodes)
        assert len(np.zeros(len(bl))[bl]) == np.sum(bl)
        assert (bl == np.array([True, False, True, False, True])).all()

    def test_expand_1D_node_vals_to_non_blacklisted_array(self):
        vals = np.array([-1, -2, -3, -4])
        ok_nodes = np.array([True, False, True, True, False, True], dtype=bool)
        def_val = 0
        should_be = np.array([-1, 0, -2, -3, 0, -4])
        expanded = dataio.expand_1D_node_vals_to_non_blacklisted_array(
            vals, ok_nodes, default_value=def_val)
        assert (should_be == expanded).all()
