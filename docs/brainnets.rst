
=====================================================
Brainnets - network analysis on fMRI data with Python
=====================================================

.. contents::

What can you do with Brainnets?
===============================

Perform network analysis of fMRI data given in the format outputted by
the bramila_ preprocessing pipeline.


Background
==========
This package contains research code for analyzing fMRI networks.
The code-basee is a collection of functions that enable somewhat pipelined network analysis of fMRI data for certain use cases.
In future, this code-base can hopefully be useful for a wider audience.

Some features of the library:
=============================
- Computation of various network metrics.
- Tools for community detection and network coarse-graining
- Plotting and visualization.
- Some documentation
- Some testing



Installation & Dependencies
===========================

This code has been tested and used only on **Linux** computers.
Here are some installation instructions for the brainnets library and dependencies.
There are (quite) many dependencies at the moment -- sorry for that!

1. Go to the directory you want to install the library in::

    cd **CODEDIR**

2. Clone the repository with::

    git clone git@github.com:rmkujala/brainnets.git

    And add these to your ~/.bashrc (or equivalent)::

        # add the main brainnets library into your PYTHONPATH:
        export PYTHONPATH=${PYTHONPATH}:**CODEDIR**/
        # add the Louvain code to your PATH, so that it can be executed from anywhere
        export PATH=${PATH}:**CODEDIR**/brainnets/Louvain/

3. Install igraph_ (both the C library as well as the Python bindings)::

    # At Aalto, we need to specify the following environment variables to get things up and running.
    export C_INCLUDE_PATH="/path/to/igraph-0.6/include/igraph:"
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/igraph-0.6/lib/
    export PYTHONPATH=${PYTHONPATH}:/path/to/python-igraph-0.6/lib/python2.7/site-packages/

4. Compile the Louvain code used by community detection::

    cd **CODEDIR**/external_code/gen-louvain/
    make


5. Get verkko_ code library (used for code parallelization & statistics).::

    cd **CODEDIR**
    git clone git@git.becs.aalto.fi:complex-networks/verkko.git
    (Required for computing some statistics + visualization/plotting tools)

6. Some parts of the code also use Matlab_, and require a working installation of mlabwrap_ for pipelining the processes.
    Now working on getting matlab.engine to work, see its installation instructions
    `here <http://se.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html>`_

    Also, the NIFTI toolbox is used for some visualizations.
    Please adjust the corresponding row in settings.py if needed.

Getting started
===============

See tutorial_, for writing scripts and commands.

Module structure
================

**Core modules:**

* `aux.py`
    many small auxillary functions
* `communities.py`
    various functions related to communities
* `compcoms.py`
    automated community computations
* `comp_helpers.py`
    generic helper functions (e.g. parallellism)
    the user should not need to use these
* `compprops.py`
    automated computation of basic network properties
* `compstats.py`
    automated statistics computations
* `complinkdistances.py`
    automatic computations of different aspects of link distances
* `config.py`
    input parameter validation and default parameters
* `dataio.py`
    file input and output
* `fname_conventions.py`
    conventions for naming files
* `gencomps.py`
    internal computations from igraph.Graph instances
* `genplots.py`
    some generic plotting functions
* `netgen.py`
    generating networks from matrices
* `plots.py`
    automated plotting scripts
* `settings.py`
    various global settings for the computations/analysis
    of the brainnets package.
    Normal user should not need to touch these.
* `visualizations.py`
    visualizations on brain slices + alluvial diagrams

**Helper script for using** triton_:

* `slurm_submit.py`
    simple script for submitting jobs to triton_

**Modulest to be removed / deprecated**

* `statistics.py`
    old statistics module, *NOT IN USE, WILL BE REMOVED!*
    Use the permtests module of verkko_ instead!
* `exports.py`
    some functions to export network stuff in other formats
    currently *outdated* = may or may not work (properly)
* `playground.py`
    Some miscellaneous stuff, not in use.


Testing
========

- The test suite for the code is located in directory ``tests``
- To test the whole suite use nosetests_::

    nosetests **CODEDIR**/tests/


Working with Triton @ Aalto University
======================================

For submitting jobs to triton_ (at Aalto University School of Science) put this to your ~/.bashrc (at triton)::

    alias trisub='python **TRITON_PATH_TO_CODEDIR**/brainnets/slurm_submit.py'

The file **TRITON_PATH_TO_CODEDIR**/brainnets/slurm_submit.py contains some code which helps you to get started using triton.

Go to the directory where you have your script and submit by typing::

    trisub script_to_submit_to_triton.py  # submits to play queue (for testing etc.)
    trisub script_to_submit_to_triton.py batch # submits to batch queue (use for real jobs)

Remember to check the contents of the slurm_submit.py, so that you understand what's going on!
Hopefully this thing helped you going with using Triton.


---------------------------------------------------------------------------

.. _python igraph: http://igraph.org/python/
.. _igraph: _python igraph
.. _bramila: https://git.becs.aalto.fi/bml/bramila
.. _verkko: https://git.becs.aalto.fi/complex-networks/verkko
.. _tutorial: tutorial.html
.. _triton: https://wiki.aalto.fi/display/Triton/Triton%20User%20Guide
.. _BECS: http://www.becs.aalto.fi
.. _nosetests: https://nose.readthedocs.org/en/latest/
.. _Matlab: http://se.mathworks.com/products/matlab/
.. _mlabwrap: http://mlabwrap.sourceforge.net/
