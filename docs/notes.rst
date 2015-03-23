====================================
Important notes on using `brainnets`
====================================

* Avg. clustering is computed using c_i=0, when k_i<2
* The computation of *weighted* network properties may not
  be currently fully reasonable, especially for connectivity
  metrics different from correlation.
* Also the avg. path length is currently not fully working with multiple components.
* Some node and link properties, like **betweenness centrality**
  etc. are not yet tested.
  (Although they should be tested in igraph)

