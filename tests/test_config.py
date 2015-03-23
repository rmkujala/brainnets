# stdlib
import unittest
import warnings
from copy import copy
# brainnets
from brainnets import config
from brainnets.testscripts import mycfg


class ConfigTest(unittest.TestCase):

    def test_get_default_config(self):
        cfg = config.get_default_config()
        assert type(cfg) is dict
        for key in config.CONFIG_KEYS:
            assert key in cfg
        for key in ["all_fnames",
                    "blacklist_fname",
                    "group_1_mat_fnames",
                    "group_2_mat_fnames",
                    "paired",
                    "outdata_dir",
                    "node_info_fname",
                    "paired"]:
            assert key in cfg
            self.assertRaises(AssertionError,
                              config.check_property, cfg, key)
        keys = ["bootstrap_samples",
                "bootstrap_coverage",
                "density",
                "density_range",
                "group_1_label",
                "group_2_label",
                "group_1_color",
                "group_2_color",
                "include_mst",
                "n_it_permutation"]
        config.require(cfg, keys)

    def test_check_full_config_with_mycfg(self):
        #testscripts/mycfg.py should always have a full copy of config
        cfg = mycfg.get_config()
        config.check_full_config(cfg)

    def test_default_config_consistency(self):
        cfg = config.get_default_config()
        cfg_copy = copy(cfg)
        cfg = config.fill_defaults(cfg)
        for key in cfg:
            assert cfg is not cfg_copy
            assert cfg[key] is cfg_copy[key]

    def test_some_defaults(self):
        cfg = {"paired": True}
        cfg = config.set_default(cfg, "n_it_permutation")
        assert cfg["n_it_permutation"] is "all"

    def test_check_property(self):
        cfg = {
            "group_1_color": "asdf",
            "group_2_color": (1.0, 1.0)
        }
        for key in cfg:
            self.assertRaises(AssertionError, config.check_property, cfg, key)

    def test_check_paired(self):
        cfg = {
            "paired": True,
            "n_it_permutation": "all"
        }
        for key in cfg:
            assert config.check_property(cfg, key)

    def test_warn_for_low_n_it_permutation(self):
        cfg = {
            "n_it_permutation": 50
        }
        for key in cfg:
            with warnings.catch_warnings(record=True) as w:
                # Cause all warnings to always be triggered.
                warnings.simplefilter("always")
                # Should trigger a warning.
                config.check_property(cfg, key)
                # Verify some things
                assert len(w) == 1
                assert "Low number of permutations" in str(w[-1].message)

    def test_warn_for_no_validity_check(self):
        cfg = {
            "this_should_not_be_a_config_param": ("plaa", 1)
        }
        for key in cfg:
            with warnings.catch_warnings(record=True) as w:
                # Cause all warnings to always be triggered.
                warnings.simplefilter("always")
                # Should trigger a warning.
                config.check_property(cfg, key)
                # Verify some things
                assert len(w) == 1
                assert "Omitting the check" in str(w[-1].message)
