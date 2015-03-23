# stdlib
import unittest
# third party
import matplotlib
matplotlib.use('Agg')
# brainnets
from brainnets import compprops, compstats, compcoms
from brainnets import plots, settings, dataio
from brainnets import visualizations as viz
from brainnets.testscripts import mycfg


class Test(unittest.TestCase):

    def setUp(self):
        self.cfg = mycfg.get_config()
        return

    def test_glob_props(self):
        compprops.compute_global_properties(self.cfg, weighted=False)
        compstats.comp_glob_prop_stats(self.cfg, False)

        compprops.compute_global_properties(self.cfg, weighted=True)
        compstats.comp_glob_prop_stats(self.cfg, True)

        plots.plot_global_uw_props(self.cfg)
        plots.plot_global_w_props(self.cfg)

    def test_node_props(self):
        compprops.comp_node_props(self.cfg)
        compstats.comp_node_pvals(self.cfg, False)

        all_fnames = self.cfg['all_fnames']
        for fname in [all_fnames[0], all_fnames[-1]]:
            fig = viz.viz_node_props_for_ind(
                fname, settings.degree_tag, self.cfg)
            assert fig

    def test_modules(self):
        compcoms.comp_louvain_communities(self.cfg)
        # viz individual
        individual = 0
        viz.viz_com_stru_for_ind(
            self.cfg['all_fnames'][individual],
            self.cfg
        )

        consensus_out_fname = (self.cfg['outdata_dir'] +
                               settings.louvain_consensus_tag +
                               "_" + "all_fnames" + "_" +
                               str(self.cfg['density']) +
                               ".pkl")
        compcoms.comp_consensus_partition(
            self.cfg,
            'all_fnames',
            consensus_out_fname)

        data = dataio.load(consensus_out_fname)
        partition = data[settings.louvain_cluster_tag]
        filtered_partition = partition[
            dataio.get_ok_nodes(self.cfg['blacklist_fname'])
        ]
        fig = viz.viz_com_structure_using_slices(filtered_partition, self.cfg)
        fig.savefig(consensus_out_fname.replace(".pkl", ".pdf"), format="pdf")

        fig = viz.comp_and_viz_cluster_diff_matrices(
            partition,
            self.cfg,
            self.cfg['all_fnames'],
            len(self.cfg['group_1_mat_fnames']),
            vmin=-5,
            vmax=5,
            suptitle="",
            recompute=True
        )
        fig.savefig(consensus_out_fname.replace(".pkl", "_coarse-grained.pdf"))
