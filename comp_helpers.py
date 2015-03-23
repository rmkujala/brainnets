import sys
import multiprocessing
# brainnets:
import settings
import dataio


def do_start(fname, blacklist_fname):
    """ Prints which work has started and returns the adj_mat """
    if settings.be_verbose:
        print "started " + fname
        sys.stdout.flush()
    adj_mat = dataio.load_adj_matrix_from_mat(fname)
    ok_nodes = dataio.get_ok_nodes(blacklist_fname)
    return adj_mat, ok_nodes


def do_end(fname, out_fname, results):
    """ Prints which work ended and saves the results"""
    dataio.save_pickle(out_fname, results)
    assert(settings.config_tag in results)
    if settings.be_verbose:
        print "finished " + fname
        sys.stdout.flush()


def run_in_parallel(work_func, arg_list, n_cpus, chunksize=1):
    """
    Run ``work_func(args)`` with n_cpus number of processors in parallel

    Returns:
        [work_func(args) for args in arg_list]
    """
    # mainly for debugging purposes and generality
    if n_cpus == 1:
        result_list = []
        for args in arg_list:
            result_list.append(work_func(args))
    else:
        pool = multiprocessing.Pool(processes=n_cpus)
        result_list = \
            pool.map_async(work_func, arg_list, chunksize=1).get(31536000)
            # To enable keyboard interuption
    return result_list
