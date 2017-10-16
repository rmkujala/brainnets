In this directory you can find an example set-up how to reproduce the network coarse-graining pipeline
used in:

Kujala, R., Glerean, E., Pan, R. K., Jääskeläinen, I. P., Sams, M. and Saramäki, J. (2016), Graph coarse-graining reveals differences in the module-level structure of functional brain networks. Eur J Neurosci, 44: 2673–2684. doi:10.1111/ejn.13392

If everything is set-up correcly, simply running

	python pipeline.py

should  

1. compute louvain communities
2. compute their consensus communities
3. match the communities (used for visualization)
4. plot the communities on brain
5. plot the alluvial
6. plot the coarse_grained matrices





