import numpy as np
from scipy.spatial.distance import squareform
import warnings


def get_n_largest_indices(array, N):
    """
    Put in a numpy array and get the *indices* of the N largest elements in the
    array. ``nan``'s are assumed as smallest values

    Parameters
    ----------
    array: a 1d numpy array
        the array whose largest indices are sought for
    N : int
        the number of indices corresponding to N largest elements
        to be returned (largest first)

    Returns
    -------
    indices : a numpy array
        array containing the **indices** of the ``N`` largest elements
        of the ``array``
    """
    dtype = array.dtype
    assert len(array.shape) is 1
    assert len(array) >= N
    if issubclass(dtype.type, np.integer):
        minval = np.iinfo(dtype).min
    else:
        if issubclass(dtype.type, np.float):
            minval = np.finfo(dtype).min
        else:
            raise
    a = array.copy()
    a[np.isnan(array)] = minval
    return np.argsort(a)[-N:][::-1]


def sliding_window_average(array, window_width):
    """
    Get the sliding window average (i.e. running mean).

    :param array: a numpy array whose first axis will
    :param window_width: the length of the window

    >>> l = [1,2,3,4,5,4,3,2,1]
    >>> ww = 2
    >>> sliding_window_average(l, ww)
    array([2, 3, 4, 4.33333, 4, 3, 2]
    >>> len(sliding_window_average) == len(l) - ww + 1
    True

    """
    window = np.ones(int(window_width)) / float(window_width)
    return np.convolve(array, window, 'valid')


def expand_up_diag_vec_to_sym_mat(vec, k=0):
    """
    Expands a vector corresponding the entries in a upper diagonal
    (numpy triu k=0) matrix to the matrix form.

    :param vec: the elements
    :param int k: either 0 (upper diagonal) or 1
        (strictly upper diagonal)

    >>> l = [1, 2, 3]
    >>> expand_up_diag_vec_to_sym_mat(l)
    array([[1, 2],
          [2, 3]])
    >>> expand_up_diag_vec_to_sym_mat(l, k=1)
    array([[ 0.,  1.,  2.],
           [ 1.,  0.,  3.],
           [ 2.,  3.,  0.]])

    """
    warnings.warn('This function is not properly tested yet')
    if k == 0:
        vec = np.array(vec)
        dtype = vec.dtype
        n = int(np.sqrt(2 * len(vec) + 1 / 4.))
        resultMat = np.zeros((n, n), dtype=dtype)
        triuIndices = np.triu_indices_from(resultMat)
        resultMat[triuIndices] = vec
        resultMat += resultMat.T
        diagIndices = np.diag_indices_from(resultMat)
        resultMat[diagIndices] = resultMat[
            diagIndices] / np.array(2, dtype=dtype)
        return resultMat
    if k == 1:
        return squareform(vec)


def extract_semi_diagonal(square_array):
    """
    From a 2D matrix/array (with size 6), this function extracts the
    numbers corresponding to the indices outputted by
    :py:func:`get_semi_diag_indices`

    :param square_array: a np.array instance (or similar), with a (n,n) shape
    """
    warnings.warn('This function is not properly tested yet')
    shape = square_array.shape
    _assertSquareArray(shape)
    return np.diag(square_array, shape[0] / 2)


def get_semi_diag_indices(square_array):
    """
    From a 2D matrix/array (with size 6), this function would now
    extract the indices of the
    matrix marked by "1"s::

        0 0 0 1 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0

    Parameters
    ----------
    square_array: a np.array
        array with a (n,n) shape
    """
    warnings.warn('This function is not properly tested yet')
    shape = square_array.shape
    _assertSquareArray(shape)
    firstAxisIndices = np.arange(0, shape[0] / 2, dtype=np.int32)
    secondAxisIndices = np.arange(shape[0] / 2, shape[0], dtype=np.int32)
    assert len(secondAxisIndices) == len(firstAxisIndices)
    return (firstAxisIndices, secondAxisIndices)


def get_up_right_indices(square_array):
    """
    From a 2D matrix/array (with size 6), this function would now
    extract the indices of the
    matrix marked by "1"s::

        0 0 0 1 1 1
        0 0 0 1 1 1
        0 0 0 1 1 1
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0

    Parameters
    ----------
    square_array: a np.array
        array with a (n,n) shape
    """
    warnings.warn('This function is not properly tested yet')
    shape = square_array.shape
    _assertSquareArray(shape)
    firstAxisIndices = np.repeat(
        np.arange(0, shape[0] / 2, dtype=np.int32), shape[0] / 2)
    secondAxisIndices = np.tile(
        np.arange(shape[0] / 2, shape[0], dtype=np.int32), shape[0] / 2)
    assert len(secondAxisIndices) == len(firstAxisIndices)
    return (firstAxisIndices, secondAxisIndices)


def get_up_left_indices_without_diag(square_array):
    """
    From a 2D matrix/array (with size 6), this function would now
    extract the indices of the
    matrix marked by "1"s::

        0 1 1 0 0 0
        1 0 1 0 0 0
        1 1 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0

    Parameters
    ----------
    square_array: a np.array
        array with a (n,n) shape
    """
    warnings.warn('This function is not properly tested yet')
    shape = square_array.shape
    _assertSquareArray(shape)
    firstAxisIndices = np.repeat(
        np.arange(0, shape[0] / 2, dtype=np.int32), shape[0] / 2 - 1)
    secondAxisIndices = []
    for diagIndex in range(0, shape[0] / 2):
        secondAxisIndices.extend(range(0, diagIndex))
        secondAxisIndices.extend(range(diagIndex + 1, shape[0] / 2))
    secondAxisIndices = np.array(secondAxisIndices)
    assert len(secondAxisIndices) == len(firstAxisIndices)
    return (firstAxisIndices, secondAxisIndices)


def get_low_right_indices_without_diag(square_array):
    """
    From a 2D matrix/array (with size 6), this function would now
    extract the indices of the
    matrix marked by "1"s::

        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 1 1
        0 0 0 1 0 1
        0 0 0 1 1 0

    Parameters
    ----------
    square_array: a np.array
        array with a (n,n) shape
    """
    warnings.warn('This function is not properly tested yet')
    shape = square_array.shape
    _assertSquareArray(shape)
    firstAxisIndices = np.repeat(
        np.arange(shape[0] / 2, shape[0], dtype=np.int32), shape[0] / 2 - 1)
    secondAxisIndices = []
    for diagIndex in range(0, shape[0] / 2):
        secondAxisIndices.extend(range(0, diagIndex))
        secondAxisIndices.extend(range(diagIndex + 1, shape[0] / 2))
    secondAxisIndices = np.array(secondAxisIndices) + shape[0] / 2
    assert len(secondAxisIndices) == len(firstAxisIndices)
    return (firstAxisIndices, secondAxisIndices)


def get_up_right_indices_without_semi_diag(square_array):
    """
    Get upper right indices without the 'semi-diagonal'.
    From a 2D matrix/array (with size 6), this function would now
    output the **indices** of the
    matrix marked by "1"s::

        0 0 0 0 1 1
        0 0 0 1 0 1
        0 0 0 1 1 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0

    Parameters
    ----------
    square_array: a np.array
        array with a (n,n) shape
    """
    warnings.warn('This function is not properly tested yet')
    shape = square_array.shape
    _assertSquareArray(shape)
    firstAxisIndices = np.repeat(
        np.arange(0, shape[0] / 2, dtype=np.int32), shape[0] / 2 - 1)
    secondAxisIndices = []
    for diagIndex in range(0, shape[0] / 2):
        secondAxisIndices.extend(range(0, diagIndex))
        secondAxisIndices.extend(range(diagIndex + 1, shape[0] / 2))
    secondAxisIndices = np.array(secondAxisIndices) + shape[0] / 2
    assert len(secondAxisIndices) == len(firstAxisIndices)
    return (firstAxisIndices, secondAxisIndices)


def extract_diagonal(square_array):
    """
    Extracts the diagonal of a matrix.

    Parameters
    ----------
    square_array: a np.array
        array with a (n,n) shape
    """
    warnings.warn('This function is not properly tested yet')
    shape = square_array.shape
    _assertSquareArray(shape)
    return np.diag(square_array)


def _assertSquareArray(shape):
    warnings.warn('This function is not properly tested yet')
    assert len(shape) == 2
    assert shape[0] == shape[1]
    assert shape[0] % 2 == 0


def expand_first_axis(array):
    """
    Expands the first axis of a numpy array.

    Parameters
    ----------
    array: a np.array
        a numpy array with at least 2 dimensions

    Returns
    -------
    out_array : a np.array
        the expanded numpy array

    Examples
    --------

    >>> a1 = [[1, 1], [2, 2]]
    >>> a2 = [[3, 3], [4, 4]]
    >>> a = np.array([a1, a2])
    >>> expand_first_axis(a)
    array([[1, 1], [2, 2], [3, 3], [4, 4]])
    >>> len(expand_first_axis(a)) == 4
    True
    """
    s = array.shape
    assert len(s) > 1, "Not enough axes to expand."
    return array.reshape((s[0] * s[1], -1))


def tex_output_array(array, out_fname, prec=2):
    """
    Outputs an numpy array in late format, which can be
    then imported to a .tex document.

    :param array: a numpy (2d) array to be outputted
    :param str out_fname: name of the created file,
        any contents will be overwritten
    :param int prec: the precision = number of decimals of the outputted
        numbers
    """
    string = " \\\\\n".join([" & ".join(map(('{0:.' + str(prec) + 'f}').format,
                                            line)) for line in array])
    f = open(out_fname, "w")
    f.write(string)
    f.close()


def get_lin_bins(n_bins, bin_min, bin_max):
    bins = np.linspace(bin_min, bin_max, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    return bins, bin_centers


def get_bin_counts(values, bins):
    return np.bincount(np.digitize(values, bins), minlength=len(bins))[1:]


def space2MNI(xyz):
    """
    Def. takes the indices of the 3d volume  (2mm space with 91x109x91)
    and returns the MNI coordinates.
    """
    xMNI = 2 * xyz[0] - 92
    yMNI = 2 * xyz[1] - 126
    zMNI = 2 * xyz[2] - 72
    return np.array([xMNI, yMNI, zMNI])
