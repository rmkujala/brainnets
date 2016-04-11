
# third party
import os
import scipy.io
import numpy as np
import cPickle as pickle
# import os
# import collections
# brainnets imports
import settings
import fname_conventions as fnc


class DataIOError(Exception):

    """
    Exception raised for errors in the dataio module.

    Parameters
    ----------
    expr : expr
        input expression in which the error occurred
    """

    def __init__(self, msg):
        self.msg = msg


def blacklist_index_fname_to_blacklist_bool_list(fname, outfname,
                                                 sample_adj_mat_fname):
    """
    Converts an (matlab based) index-based blacklist::

        [3, 10, 12, ...]

    to a list of booleans::

        [True, True, False, True, ...]

    Parameters
    ----------
    fname : str
        path to the (matlab-indexed) blacklist
    outfname : str
        path to the to-be-created bool blacklist
    sample_adj_mat_fname : str
        path to one instance of a adjacency matrix
    """
    data = load_mat(fname, squeeze_me=True)
    mat = load_adj_matrix_from_mat(sample_adj_mat_fname)
    tot_n_nodes = len(mat)
    ind_bl = data['blacklist'] - 1  # matlab indexing
    if len(ind_bl) > 0:
        assert isinstance(ind_bl[0], (int, long, np.integer))
        #assert type(ind_bl[0]) is int
    bl = index_blacklist_to_bool_blacklist(ind_bl, tot_n_nodes)
    # because matlab seems to be unable to read
    # boolean arrays created by scipy.savemat
    bl = np.array(bl, dtype=float)
    out_dict = {'blacklist': bl}
    save(outfname, out_dict)


def index_blacklist_to_bool_blacklist(ind_bl, tot_n_nodes):
    """
    Converts an python-indexed blacklist to a bool blacklist.

    Parameters
    ----------
    ind_bl : list of ints
        the index blacklist containing the indices
    tot_n_nodes : int
        the total number of nodes in the non-blacklisted data

    Examples
    --------
    >>> ind_bl = [1, 2, 4]
    >>> tot_n_nodes = 6
    >>> index_blacklist_to_bool_blacklist(ind_bl, tot_n_nodes)
    array([True, False, False, True, False, True])

    """
    blacklist = np.ones(tot_n_nodes, dtype=bool)
    blacklist[ind_bl] = False
    return blacklist


def load_pickle(fname):
    """
    Simple wrapper for loading pickled data.

    Parameters
    ----------
    fname : str
        path to the pickle filename.

    Returns
    -------
    data : the data in the file
    """
    f = open(fname, "rb")
    data = pickle.load(f)
    f.close()
    return data


def load(fname, squeeze_me=True):
    """
    Load data from pickle of mat files.

    Parameters
    ----------
    fname : str
        path to the data
    squeeze_me : bool
        whether or not to squeeze data, if .mat file is loaded
    """
    assert type(fname) is str
    format = fname[-4:]
    if format == ".pkl":
        return load_pickle(fname)
    elif format == ".mat":
        return load_mat(fname, squeeze_me=squeeze_me)
    else:
        raise DataIOError('Unknown file format: ' + format +
                          "\nSupported formats: '.pkl' and '.mat' ")


def save(fname, data_dict):
    """
    Save data in pickle (.pkl) or matlab (.mat) file.
    The format is deciphered from the fname ending (.pkl or .mat supported).

    Parameters
    ----------
    fname : str
    data_dict : dict
        a dictionary containing the data to be outputted
    """
    assert type(data_dict) == dict
    format = fname[-4:]
    if format == ".pkl":
        save_pickle(fname, data_dict)
    elif format == ".mat":
        scipy.io.savemat(fname, data_dict)
    else:
        raise DataIOError('Unknown file format: ' + format +
                          "\nSupported formats: '.pkl' and '.mat' ")


# def update_mapping(old, new):
#     """
#     Updates the old dictionary with the new.

#     Source:

#         http://stackoverflow.com/questions/3232943/
#         update-value-of-a-nested-dictionary-of-varying-depth

#     Tested it a couple of times and it seems to work -- Rainer
#     """
#     for k, v in new.iteritems():
#         if isinstance(v, collections.Mapping):
#             r = update_mapping(old.get(k, {}), v)
#             old[k] = r
#         else:
#             old[k] = new[k]
#     return old


def save_pickle(fname, data_dict):  # , overwrite=True):
    """
    Saves the data_dict to a pickle file.

    Parameters
    ----------
    fname : str
        the path where the data is stored
    data_dict : dict
        dictionary containing the data to be outputted
    """
    with open(fname, 'wb') as f:
        pickle.dump(data_dict, f, -1)


def load_adj_matrix_from_mat(fname, var_name='adj'):
    """
    Load a correlation/adjacency matrix using .mat file format

    Parameters
    ----------
    fname : str
        path to the .mat file containing the  matrix
    var_name : str
        the variable name of the matrix
    """
    assert fname[-4:] == ".mat", "Trying to load incorrect file format"
    adjMatrix = load_mat(fname, squeeze_me=False)[var_name]
    adjMatrix = np.triu(adjMatrix, 1)
    return adjMatrix


def out_put_node_prop_to_nii(nodevals, out_fname,
                             node_info_fname, blacklist_fname):
    """
    Could do this using nibabel (python) but so far just porting to the
    matlab scripts by Enrico

    Parameters
    ----------
    nodevals : numpy array
        the list of (unblacklisted) node values
    out_fname : str
        the name of the nii file
    node_info_fname : str
        the name of filename containing the information of
        the ROIS (from Enrico)
    blacklist_fname : str
        the path to .mat file listing the blacklist

    Returns
    -------
    success : bool
        whether the operation was successful or not
    """
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath(settings.package_dir + "mfiles")
    nodevals = matlab.double(list(nodevals))
    to_return = True
    if os.path.exists(settings.path_to_NIFTI_toolbox):
        eng.addpath(settings.path_to_NIFTI_toolbox)
        eng.data2niivol(nodevals, blacklist_fname, node_info_fname, out_fname, settings.ext_data_path)
    else:
        print "Path to NIFTI toolbox is not set correctly at settings.path_to_NIFTI_toolbox"
        print "Now: " + settings.path_to_NIFTI_toolbox
        print "No nii file has been produced"
        to_return = False
    eng.quit()
    return to_return


def get_ok_nodes(blacklist_fname):
    """
    Loads list of ok/bad nodes as bool array from the file defined by the
    blacklist_fname

    Parameters
    ----------
    blacklist_fname : str
        the path to the blacklist .mat file
    """
    data = load(blacklist_fname, squeeze_me=True)['blacklist']
    return np.array(data, dtype=bool)


def expand_link_val_array_to_non_bl_mat(values, blacklist_fname):
    """
    Expands a 1D array of linkwise val to the (original) 2D non-blacklisted
    matrix.
    This function is used sometimes to simplify visulization etc.
    """
    ok_nodes = get_ok_nodes(blacklist_fname).astype(np.int32)
    new_to_old_mat_indices = np.nonzero(
        np.triu(np.outer(ok_nodes, ok_nodes), 1)
    )
    valArray = np.zeros((len(ok_nodes), len(ok_nodes)))
    valArray[new_to_old_mat_indices] = values
    return valArray


def expand_1D_node_vals_to_non_blacklisted_array(
        values, ok_nodes, default_value=float("nan")):
    """
    Expands a 1d node val array to indices used before removing the
    blacklisted nodes.
    OkNodes as returned by the get_ok_nodes
    """
    assert np.sum(ok_nodes) == len(values)
    newVals = np.ones(
        len(ok_nodes), dtype=np.array(values).dtype) * default_value
    newVals[ok_nodes] = values
    assert len(newVals) == len(ok_nodes)
    return newVals


def load_mat(fname, squeeze_me=False):
    """
    Loads an arbitrary .mat file.

    Parameters
    ----------
    fname : str
        path to the .mat file
    squeeze_me : bool
        whether or not to squeeze additional
        matrix dimensions used by matlab.

    Returns
    -------
    data_dict : dict
        a dictionary with the 'mat' variable names as keys,
        and the actual matlab matrices/structs etc. as values

    Limited support is also available for HDF-matlab files.
    """
    try:
        data_dict = scipy.io.loadmat(fname, squeeze_me=squeeze_me)
    except NotImplementedError as e:
        if e.message == "Please use HDF reader for matlab v7.3 files":
            import h5py
            data = h5py.File(fname, "r")
            data_dict = {}
            for key in data.keys():
                if squeeze_me:
                    try:
                        data_dict[key] = np.squeeze(data[key])
                    except:
                        data_dict[key] = data[key]
                else:
                    data_dict[key] = data[key]
        else:
            raise e
    return data_dict


def get_blacklist_filtered_adj_mat(fname, blacklist_fname):
    """
    Get a adj_mat from which the blacklist nodes have been
    filtered out.

    Parameters
    ----------
    fname : str
        path to the .mat file containing the adj_mat
    blacklist_fname : str
        path to the blacklist (.mat file)
    """
    adj_mat = load_adj_matrix_from_mat(fname)
    ok_nodes = get_ok_nodes(blacklist_fname)
    adj_mat = adj_mat[ok_nodes, :][:, ok_nodes]
    return adj_mat


def get_blacklist_filtered_and_flattened_adj_mat(fname, blacklist_fname):
    """
    Get the blacklist filtered adj_mat in flattened (1D) form.
    See :py:func:get_blacklist_filtered_adj_mat

    Parameters
    ----------
    fname : str
        path to the .mat file containing the adj_mat
    blacklist_fname : str
        path to the blacklist (.mat file)

    Returns
    -------
    flat_adj_mat : numpy array
    """
    adj_mat = get_blacklist_filtered_adj_mat(fname, blacklist_fname)
    two_d_indices = np.triu_indices_from(adj_mat, 1)
    return adj_mat[two_d_indices].flatten()


def get_blacklist_filtered_and_flattened_adj_mats(fnames, blacklist_fname):
    """
    Get all flattened+filtered corr/adj matrices corresponding to the fnames.

    Parameters
    ----------
    fnames : list of strings
        a list of (.mat) filenames corresponding to the matrices
    """
    # Some old notes:
    # computing these beforehand + pickling -> error with unpickling
    # (a bug in the pickle library it seems, memory stuff..)
    assert type(fnames) is list
    corr_mats = []
    for filename in fnames:
        corr_mats.append(
            get_blacklist_filtered_and_flattened_adj_mat(
                filename, blacklist_fname
            )
        )
    return np.array(corr_mats)


def load_filtered_node_data(node_info_fname, blacklist_fname):
    """
    Get the blacklist filtered info about the network nodes.

    Parameters
    ----------
    node_info_fname : str
        the path to the node-info filename outputted
        from `bramila <https://git.becs.aalto.fi/bml/bramila>`_
    """
    node_info = load_mat(node_info_fname, squeeze_me=True)['rois']
    ok_nodes = get_ok_nodes(blacklist_fname)
    return node_info[ok_nodes]


def merge_and_load_props_data(fnames, props_tag,
                              props, cfg):
    """
    Load the individual properties (glob UW, glob W, node, link)
    and combine them to a joint file with structure:
    data[prop][subjectN(movie...rest)][nodeIndex] = value

    Parameters
    ----------
    fnames : list of strings
        paths to the original filenames
    props_tag : str
        the tag of the property, as specified in :py:mod:`settings`
        e.g. "global_uw_props" or "node_props" etc.
    props : list of property tags
    cfg : dict
        a brainnets config dict
    """
    data_dict = {}
    counter = 0
    n = len(fnames)
    # loop over modes
    for fname in fnames:
        ind_data = load_pickle(fnc.get_ind_fname(fname, cfg, props_tag))
        # loop over properties
        for prop in set(ind_data.keys()).intersection(props):
            # initialize new properties:
            if counter == 0:
                # l = number of nodes / percetages / or zero (modularity stuff)
                try:
                    l = len(np.array(ind_data[prop]))
                except:
                    l = 1
                # modularities are a bit falsely outputted:
                if isinstance(ind_data[prop], basestring):
                    l = 1
                data_dict[prop] = np.zeros((n, l))
            # add the data to the dict:
            data_dict[prop][counter, :] = ind_data[prop]
        counter += 1
        try:
            data_dict[settings.densities_tag] = \
                ind_data[settings.densities_tag]
        except:
            pass
    return data_dict
