# -*- coding: utf-8 -*-
# stdlib
import unittest
# third party
import igraph
import numpy as np
# brainnets
from brainnets import netgen
from brainnets.testscripts import mycfg


class Test(unittest.TestCase):

    def setUp(self):
        wMatrix = np.arange(25).reshape((5, 5)) + 1
        wMatrix = (wMatrix + wMatrix.T) + np.arange(5)
        self.wMatrix = wMatrix
        self.ok_nodes = np.array([True, True, False, True, True])
        #self.ok_nodes = np.array([0, 1, 3, 4])
        self.initMat = wMatrix.copy()
        wMatrix = wMatrix[self.ok_nodes, :]
        wMatrix = wMatrix[:, self.ok_nodes]
        self.triuMat = np.triu(wMatrix, 1)
        self.zero_density = 0.0
        self.half_density = 0.50
        self.high_density = 0.90
        self.maxNLinks = 6

        self.diffmat = np.array([[0, 11, 10, 1],
                                 [0, 0, 10, 0],
                                 [0, 0, 0, -1],
                                 [0, 0, 0, 0]])
        self.diffmat_raw = np.array(
            [[0, 11, 200, 10, 600, 1],
             [0, 0, -100, 10, -987, 0],
             [0, 0, -100, 10, 1234, 0],
             [0, 0, -100, 0, 58, -1],
             [1000, -1000, -1000, -1000, -1000, -1000],
             [0, 0, -100, 0, 1254, 0]]
        )
        #self.diffmat_oknodes = np.array([0, 1, 3, 5])
        self.diffmat_oknodes = np.array([True, True, False, True, False, True])

    # TOPOLOGICAL TESTS
    def assertHalfPercentageStuff(self, g):
        self.assertEqual(len(g.es), 3)
        links = [(g.es[0].source, g.es[0].target),
                 (g.es[1].source, g.es[1].target),
                 (g.es[2].source, g.es[2].target)]
        self.assertTrue(((2, 3) in links) != ((3, 2) in links))
        self.assertTrue(((1, 3) in links) != ((3, 1) in links))
        self.assertTrue(((0, 3) in links) != ((3, 0) in links))

    def test_makeBinNetWithoutMST(self):
        g = netgen.make_net_from_unfiltered_data(self.wMatrix, self.ok_nodes, self.zero_density, False, False)
        self.assertEqual(len(g.es), 0)

        g = netgen.make_net_from_unfiltered_data(self.wMatrix, self.ok_nodes, self.half_density, False, False)
        self.assertHalfPercentageStuff(g)

        g = netgen.make_net_from_unfiltered_data(self.wMatrix, self.ok_nodes, self.high_density, False, False)
        self.assertEqual(len(g.es), 5)
        self.assertEqual((g.es[4].source, g.es[4].target), (0, 2))

    def test_MakeBinNetWithMST(self):

        g = netgen.make_net_from_unfiltered_data(self.wMatrix, self.ok_nodes, self.half_density, True, False)
        self.assertEqual(len(g.es), 3)

    def test_MakeNetWithNegativeValues(self):
        mat = self.triuMat.copy()
        mat = np.triu(mat - 10, 1)
        g = netgen.make_net_from_unfiltered_data(self.wMatrix, self.ok_nodes, self.half_density, False, False)
        self.assertHalfPercentageStuff(g)

    def test_make_net_from_unfiltered_data(self):
        g = netgen.make_net_from_unfiltered_data(
            self.initMat - 10, self.ok_nodes,
            self.half_density, False)
        self.assertHalfPercentageStuff(g)

    def test_makeNeWithDifficultMST(self):
        g = netgen.make_net_from_unfiltered_data(self.diffmat, np.array([True, ]*4) , self.half_density, True, False)
        self.assertTrue((g.es[0].source, g.es[0].target), (0, 1))
        self.assertEqual((g.es[2].source, g.es[2].target), (0, 3))

    def test_difficultMSTWeights(self):
        g = netgen.make_net_from_unfiltered_data(self.diffmat, np.array([True, ]*4), self.half_density, True, True)
        self.assertEqual(g.es[2]["weight"], 1)
        self.assertEqual(g.es[1]["weight"], 10)
        self.assertEqual(g.es[0]["weight"], 11)

    def test_difficultMSTWeighsFromUnfilteredData(self):
        g = netgen.make_net_from_unfiltered_data(
            self.diffmat_raw,
            self.diffmat_oknodes,
            self.half_density,
            True, True
        )

        self.assertEqual(len(g.es), 3)
        self.assertEqual(g.es[2]["weight"], 1)
        self.assertEqual(g.es[1]["weight"], 10)
        self.assertEqual(g.es[0]["weight"], 11)

    def test_difficultMSTWeighsFromUnfilteredData2(self):
        print self.diffmat_raw
        g = netgen.make_net_from_unfiltered_data(
            self.diffmat_raw, self.diffmat_oknodes, 1, True, True
        )
        links = [(g.es[0].source, g.es[0].target),
                 (g.es[1].source, g.es[1].target),
                 (g.es[2].source, g.es[2].target),
                 (g.es[3].source, g.es[3].target),
                 (g.es[4].source, g.es[4].target),
                 (g.es[5].source, g.es[5].target)]
        print "links:", links

    def test_sort_links_by_weight(self):
        triumat = np.array([[0, 2, 0, 1],
                           [0, 0, 5, 0],
                           [0, 0, 0, 1],
                           [0, 0, 0, 0]])
        links = netgen.sort_links_by_weight(triumat, np.array([True, ]*4), False)
        self.assertEqual(links[0][0], 1)
        self.assertEqual(links[0][1], 2)
        self.assertEqual(len(links), 6)

    def test_make_full_weighted_net_from_weight_mat(self):
        triumat = np.array([[0, 2, 0, 1],
                           [0, 0, 5, 0],
                           [0, 0, 0, 1],
                           [0, 0, 0, 0]])
        g = netgen.make_full_weighted_net_from_weight_mat(triumat, np.array([True, ]*4))
        assert g[0, 1] == 2
        assert g[1, 0] == 2
        assert g[0, 3] == 1
        assert g[0, 2] == 0

    def test_get_graph_from_bare_data(self):
        cfg = mycfg.get_config()
        g = netgen.get_graph_from_bare_data(cfg['group_1_mat_fnames'][0],
                                            cfg['blacklist_fname'],
                                            cfg['density'],
                                            cfg['include_mst'],
                                            weighted=False)
        assert type(g) is igraph.Graph


if __name__ == '__main__':
    unittest.main()
