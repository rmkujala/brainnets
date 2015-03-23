brainnets
---------

This is currently a semi-organized collection of functions and scripts for analyzing fMRI data using networks using Python.
Currently the framework is capable of dealing with networks up to thousands of nodes.

Some of the code is currently 'tailored' for use within the Dept. of Neuroscience and Biomed. Engineering. 
However, with a little effort it can be made to meet the standards of a wider audience.

So please have a look what we've got, and let us know if you're interested in using this code.
Then we can make some efforts to make the code and its functions accessible for a wider audience.

(Main) author of the code:
Rainer.Kujala@aalto.fi



Documentation
-------------
Compile documentation with 

make docs # requires rst2html command line tool

Then open docs/build/html/index.html with your favourite browser and enjoy the documentation :)


Licencing
---------

In addition to the provided Python code, this software package distributes also open source code from other projects:

1. The ClusterPack library by A. Strehl for computing 'consensus communities' 

	http://strehl.com/download/ClusterPackV10.zip


2. The generic Louvain algorithm (by E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto) for providing high-quality partitions of networks.

	https://sites.google.com/site/findcommunities/
	http://sourceforge.net/projects/louvain/


3. R-code from A. Alexander-Bloch.:

	http://sourceforge.net/projects/brainnetworks/
	

For more information on these pieces of software and their licensing, see the subdirectories of **external_code**.

For general info on the code licensing, see LICENSE.txt





