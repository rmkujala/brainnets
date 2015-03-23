# stdlib
import unittest
# third party
import numpy as np
import igraph
# brainnets
# from brainnets.testscripts import mycfg
from brainnets import gencomps, settings


class Test(unittest.TestCase):

    def setUp(self):
        n_nodes = 5
        self.g1 = igraph.Graph(n_nodes)
        edges1 = [(0, 1), (0, 2), (1, 2), (2, 3)]
        self.g1.add_edges(edges1)
        weights1 = [1, 1, 2, 2]
        self.g1.es['weight'] = weights1

        self.g2 = igraph.Graph(n_nodes)
        edges2 = [(0, 1), (0, 2), (1, 2), (3, 4)]
        self.g2.add_edges(edges2)
        # weights2 = [1, 1, 2, 2]
        # self.g2.es['weight'] = weights2

        self.g3 = igraph.Graph(n_nodes)
        edges3 = [(0, 1), (0, 2), (0, 3), (0, 4)]
        self.g3.add_edges(edges3)

        self.g4 = igraph.Graph(7)
        edges4 = [(1, 2), (1, 3), (1, 4), (1, 5)]
        self.g4.add_edges(edges4)

        self.g5 = igraph.Graph(7)
        edges5 = [(1, 2), (4, 5)]
        self.g5.add_edges(edges5)

        self.membership_lists = np.array([
            [1, 1, 3, 4],
            [1, 1, 3, 3],
            [2, 1, 1, 1]]
        )

        # create different kinds of networks
        self.aaaeq = np.testing.assert_almost_equal

    def test_get_global_uw_props(self):
        props = settings.global_uw_props
        results = gencomps.get_global_uw_props(
            self.g1, props=props)
        rdict = dict(zip(props, results))
        assert len(self.g1.es) == 4
        assert len(self.g1.vs) == 5
        self.assertAlmostEqual(
            rdict[settings.avg_clustering_tag], (1 + 1 + 1 / 3.) / 5.0)
        self.assertAlmostEqual(rdict[settings.global_clustering_tag], 3 / 5.)
        # average over all distances over connected components!
        self.assertAlmostEqual(
            rdict[settings.avg_path_length_tag], (1 + 1 + 1 + 1 + 2 + 2) / 6.)
        self.assertAlmostEqual(rdict[settings.max_kshell_tag], 2)
        # compute assortativity manually
        degrees_linkwise = np.array([[2, 2], [2, 3], [2, 3], [3, 1]])
        degrees_linkwise_other_dir = degrees_linkwise[:, [1, 0]]
        degrees_linkwise_all = np.vstack((
            degrees_linkwise, degrees_linkwise_other_dir))
        assortativity = np.corrcoef(
            degrees_linkwise_all[:, 0], degrees_linkwise_all[:, 1])[0, 1]
        self.assertAlmostEqual(
            rdict[settings.assortativity_tag], assortativity)
        self.assertAlmostEqual(rdict[settings.max_degree_tag], 3)

    def test_get_global_w_props(self):
        props = settings.global_w_props
        results = gencomps.get_global_w_props(
            self.g1, props=props)

        rdict = dict(zip(props, results))
        assert len(self.g1.es) == 4
        assert len(self.g1.vs) == 5

        self.assertAlmostEqual(
            rdict[settings.weighted_average_path_length_tag], float("inf"))
        self.assertAlmostEqual(
            rdict[settings.max_strength_tag], 5)

        # self.assertAlmostEqual(
        #     rdict[settings.weighted_clustering_tag], -1)

        strengths_linkwise = np.array([[2, 3], [2, 5], [3, 5], [5, 2]])
        strengths_linkwise_other_dir = strengths_linkwise[:, [1, 0]]
        strengths_linkwise_all = np.vstack((
            strengths_linkwise, strengths_linkwise_other_dir))
        assortativity = np.corrcoef(
            strengths_linkwise_all[:, 0], strengths_linkwise_all[:, 1])[0, 1]
        self.assertAlmostEqual(
            rdict[settings.weighted_assortativity_tag], assortativity)
        self.assertAlmostEqual(
            rdict[settings.avg_weight_tag], 1.5)

    def test_get_node_props(self):
        props = settings.node_props
        results = gencomps.get_node_props(self.g1, props=props)
        rdict = dict(zip(props, results))
        self.aaaeq(rdict[settings.degree_tag], [2, 2, 3, 1, 0])
        self.aaaeq(rdict[settings.strength_tag], [2, 3, 5, 2, 0])
        self.aaaeq(rdict[settings.k_shell_tag], [2, 2, 2, 1, 0])
        self.aaaeq(rdict[settings.node_clustering_tag], [1, 1, 1 / 3., 0, 0])

    def test_get_link_props(self):
        pass

    def test__unwrap(self):
        d = {"a": 1, "b": 2, "c": 3}
        val_list = gencomps._unwrap(d, ["b", "c"])
        self.aaaeq(val_list, [2, 3])

    def test_get_node_props_from_mat(self):
        pass

    def test_get_link_props_from_mat(self):
        pass

    def test_comp_link_sim_mat(self):
        mat = gencomps.comp_link_sim_mat([self.g1, self.g2, self.g3])
        self.aaaeq(mat, [[4, 3, 2], [3, 4, 2], [2, 2, 4]])

    #
    # Module stuff:
    #

    def test_get_best_louvain_partition(self):
        # wigh g1
        rdict = gencomps.get_best_louvain_partition(self.g1, False, 10)
        partition = rdict[settings.louvain_cluster_tag]
        assert len(np.nonzero(partition == partition[4])) == 1
        # with g2:
        rdict = gencomps.get_best_louvain_partition(self.g2, False, 10)
        partition = rdict[settings.louvain_cluster_tag]
        assert (partition[0] == partition[1] == partition[2]
                and partition[3] == partition[4])
        # wight g4
        rdict = gencomps.get_best_louvain_partition(self.g4, False, 10)
        partition = rdict[settings.louvain_cluster_tag]
        assert len(np.nonzero(partition == partition[0])) == 1
        assert len(np.nonzero(partition == partition[-1])) == 1

        rdict = gencomps.get_best_louvain_partition(self.g5, False, 10)
        partition = rdict[settings.louvain_cluster_tag]
        assert len(np.nonzero(partition == partition[0])) == 1
        assert len(np.nonzero(partition == partition[3])) == 1
        assert len(np.nonzero(partition == partition[6])) == 1
        assert partition[1] == partition[2]
        assert partition[4] == partition[5]


    def test_get_louvain_partitions(self):
        # big random graphs should not output clear communities:
        assert_true=False
        n = 100
        d = 0.5
        for i in range(10):
            self.random_graph = igraph.Graph.Erdos_Renyi(n=n, m=int(n*(n-1)/2*d))
            rdict = gencomps.get_louvain_partitions(self.random_graph, False, 10)
            partitions = rdict[settings.louvain_cluster_tag]
            for i in range(1, len(partitions)):
                same = (np.array(partitions[0]) == np.array(partitions[i])).all()
                if not same:
                    assert_true = True
            if assert_true:
                break
        assert assert_true, "large random networks yield same output each time!"


    def test_comp_consensus_partition(self):
        """ Just simple sanity checks """
        clu1 = [1, 2, 3, 4, 5, 6]
        clu2 = [2, 3, 4, 5, 6, 7]
        clus = [clu1, clu2]
        conclu = gencomps.comp_consensus_partition(clus)
        assert len(np.unique(conclu)) == len(clu1)
        clu1 = [3, 3, 3, 4, 4, 4]
        clu2 = [1, 1, 1, 5, 5, 5]
        clus = [clu1, clu2]
        conclu = gencomps.comp_consensus_partition(clus)
        assert ((np.array(clu1) == 3) == (conclu == conclu[0])).all()
        """ Something more tricky """
        clu1 = [1, 2, 2, 3, 3, 3]
        clu2 = [1, 1, 1, 2, 2, 2]
        clus = [clu1, clu2]
        conclu = gencomps.comp_consensus_partition(clus)
        assert conclu[1] == conclu[2]
        assert conclu[3] == conclu[4] == conclu[5]
        assert conclu[1] != conclu[3]

    def comp_partition_sim_mats(self):
        partition_sim_mats = gencomps.comp_partition_sim_mats(
            self.membership_lists)
        vi = partition_sim_mats['vi']
        assert vi[0, 0] == vi[1, 1] == vi[2, 2] == 0
        assert vi[0, 1] > vi[0, 0]
        assert vi[0, 2] > vi[1, 2]
        nmi = partition_sim_mats['nmi']
        assert nmi[0, 1] < nmi[0, 0]
        assert nmi[0, 2] < nmi[1, 2]
        assert nmi[0, 0] == nmi[1, 1] == nmi[2, 2] == 1

        arand = partition_sim_mats['adjusted_rand']
        assert arand[0, 1] < arand[0, 0]
        assert arand[0, 2] < arand[1, 2]

    def test_comp_scaled_inclusivity(self):
        #  [
        #     [1, 1, 3, 4],
        #     [1, 1, 3, 3],
        #     [2, 1, 1, 1]
        # ]
        sis = gencomps.comp_scaled_inclusivity(
            self.membership_lists, normalize=False)
        self.assertAlmostEqual(sis[0], 1 + 1 / 2. + 1 / 2.)
        self.assertAlmostEqual(sis[1], 1 + 1 / 6. + 1 / 6.)
        self.assertAlmostEqual(sis[2], 1 / 2. + 1 / 3. + (2 ** 2) / (2. * 3))
        self.assertAlmostEqual(sis[2], sis[3])
        sis_norm = gencomps.comp_scaled_inclusivity(
            self.membership_lists, True)
        self.assertTrue(((0 <= sis_norm) * (sis_norm <= 1)).all())
        self.aaaeq(sis / 3, sis_norm)

    def test_comp_scaled_inclusivity_for_ref_partition(self):
        ref_partition = np.array(self.membership_lists[0])
        other_partitions = self.membership_lists[1:3]
        other_partitions = np.array(other_partitions)
        sis = gencomps.comp_scaled_inclusivity_for_ref_partition(
            ref_partition, other_partitions)
        self.assertAlmostEqual(sis[0], 1 + 1 / 2.)
        self.assertAlmostEqual(sis[1], 1 + 1. / (2 * 3))
        self.assertAlmostEqual(sis[2], 1 / 2. + 1 / 3.)
        self.assertAlmostEqual(sis[2], sis[3])
        sis_norm = gencomps.comp_scaled_inclusivity_for_ref_partition(
            ref_partition, other_partitions, True)
        self.assertTrue(((0 <= sis_norm) * (sis_norm <= 1)).all())
        self.aaaeq(sis / 2, sis_norm)

    # def test_matchClustersHungarianAlgo(self):
    #     assert False
