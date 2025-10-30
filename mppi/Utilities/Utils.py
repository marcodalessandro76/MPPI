"""
This module contains some low-level functions.
Some functions are extracted from the Futile Utils of PyBigDFT.
"""

import numpy as np
import os

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

def dict_set(inp,*subfields):
    """Ensure the provided fields and set the value

    Provide a entry point to the dictionary.
    Useful to define a key in a dictionary that may not have the
    previous keys already defined.

    Arguments:
       inp (dict): the top-level dictionary
       subfields (str,object): keys, ordered by level, that have to be retrieved from topmost level of ``inp``.
          The last item correspond to the value to be set .

    Example:

       >>> inp={}
       >>> dict_set(inp,'dft','nspin','mpol',2)
       >>> print (inp)
       {'dft': {'nspin': {'mpol': 2}}}

    """
    if len(subfields) <= 1:
        raise ValueError('invalid subfields, the sequence should be longer than one item as the last one is the value to be given')
    keys=subfields[:-1]
    tmp,key=push_path(inp,*keys)
    tmp[key]=subfields[-1]

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
    to an arbitrary depth, updating keys. The ``src`` is merged into
    ``dest``.  From :ref:`angstwad/dict-merge.py`

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
