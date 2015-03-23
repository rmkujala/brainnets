# stdlib
import unittest
# third party
# import numpy as np
# brainnets
from brainnets import settings
from brainnets import fname_conventions as fnc
from brainnets import config


class Test(unittest.TestCase):

    def setUp(self):
        return

    def test_get_fname(self):
        cfg = {"outdata_dir": "/tmp/",
               "include_mst": True, "density": 0.1}
        settings.common_links_tag
        config.check_property(cfg, "include_mst")
        config.check_property(cfg, "outdata_dir")
        config.check_property(cfg, "density")

        fname = fnc.get_fname(cfg, settings.common_links_tag)
        self.assertEqual(
            fname, "/tmp/" + settings.common_links_tag + "_with_mst.pkl")

        fname = fnc.get_fname(cfg, settings.louvain_cluster_tag)
        self.assertEqual(
            fname,
            "/tmp/" + settings.louvain_cluster_tag + "_" +
            str(cfg['density']) + "_with_mst.pkl"
        )
