"""
This module contains some low-level functions.
Some functions are extracted from the Futile Utils of PyBigDFT.
"""

from mppi.Utilities import Constants as C
import numpy as np
import matplotlib.pyplot as plt
import os

def damp_ft(ft, time, t_initial, damp_type="LORENTZIAN", eta=0.1,time_units='fs'):
    """
    Apply a damping to a function in the time domain to avoid spurious oscillations in the spectrum.
    The damping function is applied as a multiplicative factor to the function in the time domain and it is defined as

    - LORENTZIAN : exp(-|t-t_initial|*eta/hbar) 
    - GAUSSIAN : exp(-(t-t_initial)**2*(eta/hbar)**2
    
    where eta has the dimension of energy and is expressed in eV, and hbar is the reduced Planck constant
 
    Args:
        ft (:py:class:`numpy.ndarray`): The function in the time domain to be damped
        time (:py:class:`numpy.ndarray`): array with the time values 
        t_initial (:py:class:`float`): The switch on time of the external field, used as reference for the damping.      
        damp_type (:py:class:`str`): The type of damping function to apply. Can be "LORENTZIAN" or "GAUSSIAN". Default is "LORENTZIAN". 
        eta (:py:class:`float`): The damping factor to apply expressed in eV. Default is 0.1.
        time_units (:py:class:`str`): set the units to time sampling. Default is 'fs' and the other possible choice is 'au'.
    
    Return:
        :py:class:`numpy.ndarray`: The damped function in the time domain
    """

    if eta == 0.0:
        print("No damping applied to the F_t.")
        return ft
    if time_units == 'fs':
        damp_factor=eta/C.Planck_reduced_ev_fs
    elif time_units == 'au':
        damp_factor=eta/C.HaToeV # (hbar = 1 in au)
    else:
        print("Unknown time units. No damping applied.")
        return 0

    ft_damped=np.empty_like(ft)
    if damp_type.upper() == "LORENTZIAN":
        ft_damped[:]=ft[:]*np.exp(-abs(time[:]-t_initial)*damp_factor)
    elif damp_type.upper() == "GAUSSIAN":
        ft_damped[:]=ft[:]*np.exp(-(time[:]-t_initial)**2*damp_factor**2)
    else:
        print("Unknown damping type. No damping applied.")
        return 0

    return ft_damped

def Plot_Array(xvalues, data, data2=None, xlim=None, label=None, label2=None, figsize=(8,6)):
    """
    Plot real array data.
    
    Supports:
    - data: 1D or 2D
    - data2: optional, same first dimension as data
    
    Args:
        xvalues (:py:class:`array`): x values
        data (:py:class:`array`): real data (1D or 2D)
        data2 (:py:class:`array`, optional): second dataset
        xlim (:py:class:`tuple`, optional): x-axis limits
        label (:py:class:`str`, optional): label for data
        label2 (:py:class:`str`, optional): label for data2
        figsize (:py:class:`tuple`, optional): figure size for the plot (width, height. Default is (8, 6)
    """
    
    char_size = 12
    fig, axes = plt.subplots(figsize=figsize)

    if label is None:
        label = 'F'
    if data2 is not None and label2 is None:
        label2 = 'F2'

    title = f'Plot of the array data {label}'
    if data2 is not None:
        title += f' and {label2}'
    axes.set_title(title,fontsize=char_size
    )

    data = np.asarray(data)
    def plot_dataset(data, base_label):
        if data.ndim == 1:
            axes.plot(xvalues, data, label=base_label)

        elif data.ndim == 2:
            n_cols = data.shape[1]
            for i in range(n_cols):
                lbl = f"{base_label}[{i}]"
                axes.plot(xvalues, np.real(data[:, i]), label=lbl)
        else:
            raise ValueError("data must be 1D or 2D array")

    plot_dataset(data, label)
    if data2 is not None:
        data2 = np.asarray(data2)
        if data2.shape[0] != data.shape[0]:
            raise ValueError("data and data2 must have the same first dimension")
        plot_dataset(data2, label2)

    
    axes.legend(fontsize=char_size - 2)
    if xlim is not None:
        axes.set_xlim(xlim)

    plt.show()

def Plot_3dArray(xvalues, data, data2=None, xlim=None, label=None, label2=None):
    """
    Plot array data in the three Cartesian directions.
    
    Args:
        xvalues (:py:class:`array`): x values
        data (:py:class:`array`): 2D array with shape (3, N)
        data2 (:py:class:`array`, optional): second dataset, same shape as data
        xlim (:py:class:`tuple`, optional): x-axis limits
        label (:py:class:`str`, optional): label for data
        label2 (:py:class:`str`, optional): label for data2
    """
    
    char_size = 12
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 10))

    if label is None:
        label = 'F'
    if data2 is not None and label2 is None:
        label2 = 'F2'

    axes[0].set_title(
        'Plot of the array in the three cartesian directions',
        fontsize=char_size
    )

    data = np.asarray(data)
    if data.shape[0] != 3:
        raise ValueError("data must have shape (3, N)")

    directions = ['x', 'y', 'z']

    def plot_dataset(data, base_label):
        for i in range(3):
            axes[i].plot(
                xvalues,
                data[i],
                label=f"{base_label}_{directions[i]}"
            )

    plot_dataset(data, label)
    if data2 is not None:
        data2 = np.asarray(data2)

        if data2.shape != data.shape:
            raise ValueError("data2 must have the same shape as data")

        plot_dataset(data2, label2)

    for ind in range(3):
        axes[ind].legend(fontsize=char_size - 2)

    if xlim is not None:
        for ind in range(3):
            axes[ind].set_xlim(xlim)

    plt.tight_layout()
    plt.show()

def Plot_ComplexArray(xvalues, data, data2=None, xlim=None, label=None, label2=None):
    """
    Plot real and imaginary parts of complex data.
    
    Supports:
    - data: 1D or 2D
    - data2: optional, same first dimension as data
    
    Args:
        xvalues (:py:class:`array`): x values
        data (:py:class:`array`): complex data (1D or 2D)
        data2 (:py:class:`array`, optional): second dataset
        xlim (:py:class:`tuple`, optional): x-axis limits
        label (:py:class:`str`, optional): label for data
        label2 (:py:class:`str`, optional): label for data2
    """
    
    char_size = 12
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))

    if label is None:
        label = 'F'
    if data2 is not None and label2 is None:
        label2 = 'F2'

    axes[0].set_title(
        'Plot of the real and imaginary parts of complex data',
        fontsize=char_size
    )

    data = np.asarray(data)
    def plot_dataset(data, base_label):
        if data.ndim == 1:
            axes[0].plot(xvalues, np.real(data), label='Re ' + base_label)
            axes[1].plot(xvalues, np.imag(data), label='Im ' + base_label)

        elif data.ndim == 2:
            n_cols = data.shape[1]
            for i in range(n_cols):
                lbl = f"{base_label}[{i}]"
                axes[0].plot(xvalues, np.real(data[:, i]), label='Re ' + lbl)
                axes[1].plot(xvalues, np.imag(data[:, i]), label='Im ' + lbl)
        else:
            raise ValueError("data must be 1D or 2D array")

    plot_dataset(data, label)
    if data2 is not None:
        data2 = np.asarray(data2)
        if data2.shape[0] != data.shape[0]:
            raise ValueError("data and data2 must have the same first dimension")
        plot_dataset(data2, label2)

    for ind in range(2):
        axes[ind].legend(fontsize=char_size - 2)

    if xlim is not None:
        for ind in range(2):
            axes[ind].set_xlim(xlim)

    plt.tight_layout()
    plt.show()

def file_to_list(filename,skip='#'):
    """
    Read the filename and append all the lines that do not start
    with the skip string, to a list. Empty lines are skipped

    Args:
        filename (str): name of the file
        skip (str): first elements of the skipped lines
    """
    lines = []
    with open(filename) as f:
        for l in f:
            if not l.startswith(skip) and len(l.strip()) > 0:
                lines.append(l)
    return lines

def floats_from_string(line,sep=None):
  """
  Split a string using blank spaces and convert the elements to float. If an element
  cannot be converted it is skipped.

  Args:
        line (:py:class:`string`): string that contains a line of the file
        sep (:py:class:`string`) : Delimiter at which splits occur. If `None` the string is splitted at whitespaces

  Return:
        :py:class:`list` : list with the floats extracted from the line

  """
  line_float = []
  for value in line.split(sep=sep):
      try: line_float.append(float(value))
      except ValueError: pass
  return line_float

def file_parser(filename,skip='#',sep=None):
    """
    Parse a file. All the lines the start with the skip string are skipped.
    The lines of the file are converted in floats (elements are separated using the
    sep parameter). Elements that cannot be converted to float are ignored.

    Args:
        filename (:py:class:`string`): name of the file
        skip (:py:class:`string`): first elements of the skipped lines
        sep (:py:class:`string`) : Delimiter at which splits occur. If `None` the string is splitted at whitespaces

    Return:
        :py:class:`array' : array with the floats extracted from the file sorted in columns. So columns[0] contains
            the first column of the file and so on


    """
    lines = file_to_list(filename)
    splitted = []
    for line in lines:
         splitted.append(floats_from_string(line,sep=sep))

    columns = np.array(splitted).transpose()
    return columns

def array_dump_to_file(data, filename, columns=None, sep="       ", decimals=6):
    """
    Writes a 2D Python array to a text file.
    
    Args:
        data (list of lists): 2D array to write to file
        filename (str): name of the output text file
        columns (list of str, optional): list of column names. If None, default names "col1", "col2", ... are used.
        sep (str, optional): separator string between columns. Default is 7 spaces. 
        decimals (int, optional): number of decimal places for numeric values. Default is 6.    

    """
    ncols = len(data[0])

    # Generate column names automatically if needed
    if columns is None:
        columns = [f"col{i+1}" for i in range(ncols)]
    else:
        if len(columns) != ncols:
            raise ValueError("Number of column names does not match data width.")

    # Validate row lengths
    for row in data:
        if len(row) != ncols:
            raise ValueError("All rows must have the same number of columns.")

    # Detect numeric or string
    def is_numeric(x):
        return isinstance(x, (int, float))

    # Format function for numbers
    def fmt(x):
        if is_numeric(x):
            return f"{x:.{decimals}f}"
        return str(x)

    # Transform everything into strings
    str_data = [[fmt(x) for x in row] for row in data]

    # Compute column widths
    col_widths = [
        max(len(columns[i]), max(len(str_data[r][i]) for r in range(len(str_data))))
        for i in range(ncols)
    ]

    # Format header (right aligned)
    def format_header():
        return sep.join(columns[i].rjust(col_widths[i]) for i in range(ncols))

    # Format a row (right aligned for all values)
    def format_row(row):
        return sep.join(row[i].rjust(col_widths[i]) for i in range(ncols))

    # Write the file
    with open(filename, "w") as f:
        f.write("# " + format_header() + "\n")   # header line
        for row in str_data:
            f.write(format_row(row) + "\n")


def convertTonumber(x):
    """
    Check if the input string can be converted to an
    integer or to a float variable

    Args:
        x (:py:class:`string`) : string to be converted

    """
    try :
        if x.lstrip('-').isnumeric() : return int(x)
        else : return float(x)
    except (TypeError, ValueError):
        return x

def is_point_near_line(point, line_start, line_end, tolerance):
    """
    Check if a point is near a line segment defined by two endpoints.

    Parameters:
        point (tuple): The point to check (x0, y0)
        line_start (tuple): Starting point of the line (x1, y1)
        line_end (tuple): Ending point of the line (x2, y2)
        tolerance (float): Maximum allowed distance to consider the point "near"

    Returns:
        bool: True if the point is within the tolerance from the line segment, False otherwise
    """
    # Convert input points to NumPy arrays
    P = np.array(point)
    A = np.array(line_start)
    B = np.array(line_end)

    # Vector from A to B and from A to P
    AB = B - A
    AP = P - A

    # Length of the segment
    norm_AB = np.linalg.norm(AB)
    if norm_AB == 0:
        # Degenerate case: A and B are the same point
        return np.linalg.norm(P - A) <= tolerance

    # Projection scalar of AP onto AB
    t = np.dot(AP, AB) / np.dot(AB, AB)

    # Clamp t to [0, 1] to stay within the segment
    t = max(0, min(1, t))
    closest_point = A + t * AB

    # Distance from P to the closest point on the segment
    distance = np.linalg.norm(P - closest_point)

    return distance <= tolerance

# Can be removed
# def dict_set(inp,*subfields):
#     """Ensure the provided fields and set the value

#     Provide a entry point to the dictionary.
#     Useful to define a key in a dictionary that may not have the
#     previous keys already defined.

#     Arguments:
#        inp (dict): the top-level dictionary
#        subfields (str,object): keys, ordered by level, that have to be retrieved from topmost level of ``inp``.
#           The last item correspond to the value to be set .

#     Example:

#        >>> inp={}
#        >>> dict_set(inp,'dft','nspin','mpol',2)
#        >>> print (inp)
#        {'dft': {'nspin': {'mpol': 2}}}

#     """
#     if len(subfields) <= 1:
#         raise ValueError('invalid subfields, the sequence should be longer than one item as the last one is the value to be given')
#     keys=subfields[:-1]
#     tmp,key=push_path(inp,*keys)
#     tmp[key]=subfields[-1]

def dict_get(inp,*subfields):
    """Find the value of the provided sequence of keys in the dictionary, if available.

    Retrieve the value of the dictionary in a sequence of keys if it is available.
    Otherwise it provides as default value the last item of the sequence ``subfields``.

    Args:
       inp (dict): the top-level dictionary. Unchanged on exit.
       subfields (str,object): keys, ordered by level, that have to be retrieved from topmost level of ``inp``.
              The last item correspond to the value to be set.

    Returns:
       The value provided by the sequence of subfields if available, otherwise the default value given as the last item of the ``subfields`` sequence.

    """
    if len(subfields) <= 1:
        raise ValueError('invalid subfields, the sequence should be longer than one item as the last one is the value to be given')
    tmp=inp
    keys=subfields[:-1]
    val=subfields[-1]
    for key in keys:
        tmp=tmp.get(key)
        if tmp is None: return val
    return tmp

def sort_lists(sort_by,ascending,*lists):
    """
    Sort lists altogether following the lists indicated by the ``sort_by`` index.

    Args:

       sort_by (int):  the index of the list which has to be taken as reference for sorting
       ascending (bool): Sort is performed in ascending order if True

       *lists: sequence of lists to be mutually sorted. They have to be of the same length.

    Returns:
       tuple of sorted lists

    Example:

        >>> l1=[5,3,4]
        >>> l2=['c','t','q']
        >>> l3=[6,3,7]
        >>> print (sort_lists(0,True,l1,l2,l3))
        >>> print (sort_lists(2,True,l1,l2,l3))
        [(3, 4, 5), ('t', 'q', 'c'), (3, 7, 6)]
        [(3, 5, 4), ('t', 'c', 'q'), (3, 6, 7)]
    """
    import operator
    return zip(*sorted(zip(*lists), reverse=not ascending, key=operator.itemgetter(sort_by)))


def dict_merge(dest, src):
    """
    Recursive dict merge. Inspired by :meth:`dict.update`, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``src`` is merged into ``dest``.

    Arguments:
       dest (dict): dict onto which the merge is executed
       src (dict): dict merged into dest

    """
    import collections
    for k, v in src.items():
        if (k in dest and isinstance(dest[k], dict)
                and isinstance(src[k], collections.Mapping)):
            dict_merge(dest[k], src[k])
        else:
            dest[k] = src[k]

def ensure_dir(file_path):
    """
    Guarantees the existance on the directory given by the (relative) file_path

    Args:
       file_path (str): path of the directory to be created

    Returns:
       bool: True if the directory needed to be created, False if it existed already
    """
    directory = file_path
    created=False
    if not os.path.exists(directory):
        os.makedirs(directory)
        created=True
    return created
