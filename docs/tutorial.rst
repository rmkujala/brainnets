==================
Brainnets tutorial
==================

(Part of the brainnets documentation. See readme_)

.. contents::

Required data (outputted from bramila_)
=======================================

Example data can be found from **CODEDIR**/testdata/

node_info_fname:

- a matlab file containing information about the network nodes
- example: roi_info.mat

corr_mats:

- correlation matrices
- examples: movie_*X*.mat and rest_*X*.mat

blacklist:

- node blacklist, in format::

	[True, False, True, True, True, False, True, ... ]

  invalid nodes should be masked false
  **this is actually more of a whitelist**

- If the original blacklist is in the form of indices::

	[1, 2, 5, ..] # matlab indexing

  Use :py:func:`dataio.blacklist_index_fname_to_blacklist_bool_list`
  to first convert it to the required format.

**brainnets config dictionary**:

- Brainnets functionality is mostly steered with a config object
  (which is just a python dictionary)
- It is often enough for the user to give only the config object as a
  parameter to the functions
- Use :py:func:`config.check_full_config` which makes a *lightweight*
  sanity check of the input parameters (i.e. the existence of paths/files
  are checked, but in general the contents of any files are not validated)
- The meaning of the input parameters is shortly presented in the code,
  where all the different :py:attr:`config.CONFIG_KEYS` are listed
- For an example config dictionary, look at testscritpts/myconfig.py
- For a new project, the best practice is probably to store the config
  object to a separate python module (something like in `myconfig.py`).
- After you think you've got everything set up, check the validity of
  your config dicti using :py:func:`config.check_full_config` to avoid
  simple, easily avoidable problems.


Basic pipelines
===============
See example data processing pipelines from
``tests/test_run_pipes.py``
(There they are in form of a test, more explicit examples/demos to become.)


Notes on data input and output
==============================
The file format used internally by ``brainnets`` is Python's pickle_
format.
This choise has been made to enable easy input and output by
avoiding unnecessary file parsing.
("pickling" is actually the same as dumping the objects from memory.)
Most of the time this does not require anything extra from the user.

To read a pickle-file, do the following (e.g. in ``ipython`` console)::

	import pickle
	data = pickle.load(open("my_pickle_file.pkl", "r"))

or using brainnets, simply::

	from brainnnets import dataio
	data = dataio.load_pickle("my_pickle_file.pkl")

	# This works as well:
	data = dataio.load("my_pickle_file.pkl")


----------------------------------------

.. _readme: readme.html
.. _pickle: https://docs.python.org/2/library/pickle.html
.. _bramila: https://git.becs.aalto.fi/bml/bramila