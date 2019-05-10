# -*- coding: utf-8 -*-

"""

The below functions can be used to import delimited data files into Numpy or
Matlab database format.

"""

import argparse
import copy
import glob
import math
import os
import re
from enum import Enum

import numpy as np
import pkg_resources
# pylint: disable=no-member
import scipy.io


class _Colors:
    """

    A collection of colors that can be used to highlight terminal outputs.

    """

    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class _TextSnippets(Enum):
    """

    Text snippets to be used when merging delimited files.

    """
    header = "This file was automatically generated using the merge_del\n" \
             "function of the Python tribology package, version {}.\n" \
             "\n" \
             "See here for more information:\n" \
             "https://pypi.org/project/tribology/\n"\
             "\n"\
             "The file contains data from the following source files " \
             "(in  order):\n"

    seperator = "\n" \
                "Beginning of file:\n" \
                "{}\n"


def __make_dir(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)
    return dirpath


def __get_outpath(outdir):
    if outdir:
        outpath = __make_dir(outdir)
    else:
        outpath = os.getcwd()
    return outpath


def __get_outfile(in_file, idx, out_ext):
    fname = ''.join(in_file.split('.')[:-1]).split(os.sep)[-1]
    return '{}-{}.{}'.format(fname, str(idx), out_ext)


def __num_char(char):
    return bool(char.isdigit() or char == '-')


def split_del(file, deli='\t', ext='txt', cmin=3, hspan=1, outdir=None,
              force=False):
    """

    Split a delimited data file into several separate data files, if the file
    contains more than one block of data. Blocks of data are typically
    separated by at least one line of column headers. The first data column
    of each data block has to be numeric.

    This function is meant to be used on data files where different blocks of
    data have different numbers of columns or different column headers. After
    splitting the data file into individual data files, import methods like
    :code:`import_del` can be used on the individual files. If all data should
    be merged into a single database afterwards, the :code:`merge_npz` function
    can be used.

    Parameters
    ----------
    file: str
        Path to the data file.
    deli: str, optional
        Delimiter used to separate data columns in :code:`file`
    ext: str, optional
        File extension of output files. Default is :code:`txt`
    cmin: int, optional
        Minimum number of columns that a line of data needs to have in order to
        be classified as data.
    hspan: int, optional
        Maximum number of non-data lines above each data block that should be
        written to individual data files (usually equal to number of lines
        spanned by the column headers).
    outdir: str, optional
        Path to output directory. Default is current working directory.
    force: bool
        If True, existing output files will be overwritten. Will raise an
        exception if file exists and force is False.

    Returns
    -------
    outfiles: list
        Paths to output files.

    """
    outpath = __get_outpath(outdir)

    outfiles = []
    idx = 0
    f_out = None
    write = False
    to_write = []

    with open(file) as infile:
        for line in infile:

            # if first character of line is not numeric
            if not __num_char(line[0]):
                write = False
                to_write.append(line)

                while len(to_write) > hspan:
                    del to_write[0]

            else:
                # if numeric line has at least 'cmin' columns
                if len(line.split(deli)) >= cmin and not write:
                    write = True

                    idx += 1
                    f_out = os.sep.join([outpath,
                                         __get_outfile(file, idx, ext)])
                    if f_out not in outfiles:
                        outfiles.append(f_out)

                    if os.path.isfile(f_out):
                        if force:
                            os.remove(f_out)
                        else:
                            raise OSError("output file exists. "
                                          "use argument 'force' to overwrite.")

            if write and f_out:
                with open(f_out, "a") as out:
                    for element in to_write:
                        out.write(element)
                    to_write = []
                    out.write(line)

    return outfiles


def __verify_merge(in_files, accum):
    """

    Check if all npz files have same set of keys and contain all keys in accum.
    Throw exception if not.

    Parameters
    ----------
    in_files: list
        Paths to database files to merge. Files are merged in order.
    accum: list
        Database keys for which values should be accumulated. Values must be
        numeric.

    """
    ref_keys = []
    for idx, file in enumerate(in_files):
        keys = sorted(np.load(file).keys())
        if idx == 0:
            ref_keys = copy.deepcopy(keys)
        if keys != ref_keys:
            raise KeyError('keys in npz databases 0 and {} differ'.format(idx))
        if accum and not all(key in keys for key in accum):
            raise KeyError('key(s) defined in accum not in npz database {}'
                           .format(file))


def merge_npz(in_files, accum=None, safe=True):
    """

    Merge npz databases by concatenating all databases in :code:`in_files`.
    Databases are concatenated in the order given in :code:`in_files`.

    Database keys for which values are to be accumulated can be given as a list
    using the :code:`accum` argument. For examples, if all databases have the
    key :code:`time`, then :code:`accum=['time']` will produce a continuous
    time axis, adding the last time value of the first database to all time
    values of the second database (and so on).

    Parameters
    ----------
    in_files: list
        Paths to database files to merge. Files are merged in order.
    accum: list
        Database keys for which values should be accumulated. Values must be
        numeric.
    safe: bool
        If True, checks will be performed to ensure that all databases share the
        exact same set of keys and that all keys in :code:`accum` are in all
        databases. An exception (type KeyError) will be raised if not.

    Returns
    -------
    merged: dict
        Merged data.

    """

    if safe:
        __verify_merge(in_files, accum)

    merged = {}
    for file in in_files:
        in_dat = np.load(file)
        for key in in_dat.keys():
            if key in merged:
                if accum and key in accum:
                    merged[key] = np.append(merged[key],
                                            in_dat[key] + merged[key][-1])
                else:
                    merged[key] = np.append(merged[key], in_dat[key])
            else:
                merged[key] = in_dat[key]

    return merged


def __get_version(package):
    """

    Get the version of a Python package.

    Parameters
    ----------
    package: str
        The name of the package

    Returns
    -------
    Version number as string.

    """
    return pkg_resources.get_distribution(package).version


def __long_substr(strings):
    """

    Returns longest common substring of list of strings. taken from:
    # https://stackoverflow.com/questions/2892931/longest-common-substring-
    from-more-than-two-strings-python

    Parameters
    ----------
    strings: list
        A list of strings.

    Returns
    -------
    substr: str
        The longest common substring of all list elements. For a list with only
        one element, the list element is returned; for an empty list, and empty
        string is returned.

    """
    substr = ''
    if len(strings) > 1 and len(strings[0]) > 0:
        for i in range(len(strings[0])):
            for j in range(len(strings[0]) - i + 1):
                if j > len(substr) and all(strings[0][i:i + j] in x for x in
                                           strings):
                    substr = strings[0][i:i + j]
    return substr


def merge_del(in_files, out_file=None):
    """

    Merge several delimited data files into a single file. The merged
    file contains all data from the data files, in the order given in the
    :code:`in_files` argument.

    No checks are performed to ensure that the data files
    have a compatible format, for example the same number of data columns.

    Parameters
    ----------
    in_files: list
        File paths to the files to be merged. Files will be merged in order.
    out_file: str, optional
        Path to output file, including file extension. If no path is provided,
        a file name is generated based on the input file names.

    Returns
    -------
    out_file_abs: str
        Absolute path to the merged file.

    """
    if len(in_files) == 0:
        raise ValueError('need at least one file to merge')

    in_files_abs = [os.path.abspath(file) for file in in_files]

    if out_file:
        out_file_abs = os.path.abspath(out_file)
    else:
        out_file = __long_substr(in_files_abs).split('.')[0]
        out_file_abs = out_file + 'xxx-merged.txt'

    max_len_path = max(len(file) for file in in_files_abs)

    with open(out_file_abs, "w") as txt_file:

        # write header
        txt_file.write(str(_TextSnippets.header.value).format(
            __get_version("tribology")))
        for in_file in in_files_abs:
            txt_file.write(in_file + "\n")

        # write files
        for in_file in in_files_abs:
            txt_file.write('\n' + '#' * max_len_path)
            txt_file.write(str(_TextSnippets.seperator.value).format(in_file))
            txt_file.write('#' * max_len_path + '\n')
            with open(in_file) as file:
                for line in file:
                    txt_file.write(line)

    return out_file_abs


def __print_status(message, status_color=_Colors.ENDC):
    """

    Print a color-coded message to the terminal.

    Parameters
    ----------
    message: str
        The message to print to the terminal.
    status_color:
        The color in which to print the message.

    """
    print(status_color + message + _Colors.ENDC)


def __is_floatable(num):
    """

    Check if 'num' can be converted to float. If yes, return :code:`True`, else
    return :code:`False`.

    """
    try:
        float(num)
        return True
    except ValueError:
        return False


def __to_float(num):
    """

    Try to convert 'num' to float, return 'num' if it's not possible, else
    return converted :code:`num`.

    """
    try:
        float(num)
        return float(num)
    except ValueError:
        return num


def __assemble_data_table(num_data_tables, max_num_data_length):
    """

    Assemble the complete data table from a list of data tables.

    """
    num_data = np.zeros((
        (len(num_data_tables) - 1) * max_num_data_length +
        num_data_tables[-1].shape[0],
        num_data_tables[-1].shape[1]), dtype=object)
    for idx, data_table in enumerate(num_data_tables):
        # do this for all but the last data table
        if idx + 1 < len(num_data_tables):
            num_data[idx * max_num_data_length:
                     (idx + 1) * max_num_data_length, :] = data_table
        # do this for the last data table
        else:
            num_data[idx * max_num_data_length:, :] = data_table
    return num_data


def __write_to_out_dict(num_data, column_headers, pcs=False):
    """

    Extract the data columns from the num_data array and write them to a
    dictionary.

    Parameters
    ----------
    num_data: ndarray
        The data extracted from the delimited file, stored in a single table.
    column_headers: ndarray
        The column headers corresponding to the columns in :code:`num_data`

    Returns
    -------
    out_dict: dict
        A dictionary containing all data that is to be saved to the output
        database. Keys are based on column headers, values are data columns of
        num_data.

    """
    out_dict = {'column_headers': column_headers}
    for idx, column in enumerate(column_headers):
        # explicitly take care of the fact that PCS forgot a '\tab' character in
        # their data export implementation
        if column == 'image_file_name' and \
                math.isnan(float(num_data[0, idx])) and not \
                column_headers[column_headers.tolist().index(column) - 1] and \
                pcs is True:
            out_dict[column] = num_data[:, idx - 1].astype(object)[:, None]
        # take care of all other columns
        # if empty data columns are not padded with tabs
        elif column:
            if idx >= num_data.shape[1]:
                out_dict[column] = np.zeros(num_data.shape[1]) * float('nan')
            else:
                # if data is of numeric type
                if __is_floatable(num_data[0, idx]):
                    out_dict[column] = num_data[:, idx].astype(float)[:, None]
                # if data is of other type (string)
                else:
                    out_dict[column] = num_data[:, idx].astype(object)[:, None]
    return out_dict


def __process_header(heads):
    """

    Process the column headers by removing special characters and converting to
    Matlab-optimized data type.

    Parameters
    ----------
    prev_line: list of strings
        The column headers of the delimited file.

    Returns
    -------
    col_heads: ndarray (dtype = object)
        The re-formated column headers.

    """

    merge = []
    # merge colum headers if they span several lines
    for i in range(len(heads[0])):
        merge.extend([' '.join([heads[row][i] for row in range(len(heads))])])

    # replace non-alphanumeric characters and trailing underscores
    col_heads = [re.sub(r"\W+", '_', item.lower()).strip('_') for item in merge]

    # convert data type for easy matlab export
    col_heads = np.asarray(col_heads, dtype='object')

    return col_heads


def __process_data(split_line, num_dat, max_len, num_data_tables):
    """

    Append a data line to the current data table. If the length of the current
    data table exceeds the maximum permitted data table length, save the current
    data table to a list of data tables and initialise a new one.

    Parameters
    ----------
    split_line: ls
        The data that is to be appended to the table.
    num_dat: ndarray
        The current data table to which the last line of data was appended.
    max_len: positive int
        The maximum length of a data table.
    num_data_tables: ls
        The complete list of data tables.

    Returns
    -------
    num_dat: ndarray
        The data table to which the current line of data was appended.

    """
    # if data table becomes large, make new data table and add old
    # table to table list (for speed)
    if num_dat.shape[0] == max_len:
        num_data_tables.append(num_dat)
        num_dat = np.asarray(
            [__to_float(item.rstrip('\n')) for item in
             split_line]).reshape((1, len(split_line)))
    # else simply append to data table
    else:
        num_dat = np.append(num_dat, np.asarray(
            [__to_float(item.rstrip('\n')) for item in split_line])
                            .reshape((1, len(split_line))), axis=0)
    return num_dat


def __process_file(in_file, dec_mark, deli, pad=0, colheadlines=1):
    """

    Extract data from a delimited text file and return a dictionary containing
    all data.

    Parameters
    ----------
    in_file: str
        The file handle of the delimited file that is to be imported.
    dec_mark: str
        The decimal mark of the data file.
    deli: str
        The delimiter used to separate data columns in the delimited file.
    pad: positive int
        Ignore the first :code:`n` leading columns in the delimited file, where
        :code:`n = pad`. For example, if pad = 8, the first 8 columns
        are ignored.

    Returns
    -------
    out_dict: dict
        A dictionary containing all data that is to be saved to the output
        database. Keys are based on column headers, values are data columns of
        num_data.

    """
    max_len = 1000
    num_dat = []
    col_heads = []
    num_data_tables = []
    prev_lines = []

    with open(in_file) as dat_file:
        for line in dat_file:
            split_line = line.replace(dec_mark, '.').split(deli)

            if len(split_line) > pad:
                split_line = split_line[pad:]

            # get rid of trailing newline characters
            if split_line[-1] == '\n':
                split_line[-1] = ''

            # check if first character is not (digit or minus symbol (hyphen))
            # to identify non-data lines. skip non-data lines.
            if not (line[0].isdigit() or line[0] == '-') or \
                    len(split_line) <= 1:
                if split_line != ['']:
                    prev_lines.append(split_line)
                    if len(prev_lines) > colheadlines:
                        del prev_lines[0]
                continue

            # if line contains data, split line into data fields, fill empty
            # fields with 'nan'
            split_line[:] = (item or 'nan' for item in split_line)
            # if this is the first data-containing line...
            if not len(col_heads):
                # get the column headers
                col_heads = __process_header(prev_lines)
                # write the first line to the data table
                num_dat = np.asarray(
                    [__to_float(item.rstrip('\n'))
                     for item in split_line]).reshape((1, len(split_line)))
            else:
                num_dat = __process_data(split_line, num_dat, max_len,
                                         num_data_tables)

    # assemble the complete data table and create output dictionary
    num_data_tables.append(num_dat)
    num_dat = __assemble_data_table(num_data_tables, max_len)

    return num_dat, col_heads


def __get_file_handles(in_dir, ext, recursive=False):
    """

    Get file handles for all delimited files that are to be imported.

    Parameters
    ----------
    in_dir: str
        The directory in which the delimited files are stored.
    ext: str
        The file extension of the delimited files.
    recursive: bool, optional
        If :code:`True`, delimited files are imported for all child directories
        of :code:`directory` (including :code:`directory`). If :code:`False`,
        only files in :code:`directory` are imported. Default is :code:`False`.

    Returns
    -------
    in_files: ls of strings
        The file handles to all delimited files that are to be imported.

    """
    if not recursive:
        in_files = sorted(glob.glob('{}{}*.{}'.format(in_dir, os.sep, ext)))
    else:
        in_files = []
        dir_list = [x[0] + os.sep for x in os.walk(in_dir)]
        for directory in dir_list:
            in_files.extend(sorted(glob.glob('{}*.{}'.format(directory, ext))))
        # in_files = [f.replace(in_dir, '').lstrip(os.sep) for f in in_files]
    return in_files


def __save_out_file(out_file, out_dict, out_ext):
    """

    Save the imported data to an output database, either in Numpy or Matlab
    format.

    Parameters
    ----------
    out_file: str
        A handle to the output file that was generated during import.
    out_dict: dict
        The output data stored in a dictionary where keys correspond to column
        headers, values correspond to data.
    out_ext: str
        The file extension (format) of the output file. Options are :code:`npz`
        for Numpy format and :code:`mat` for Matlab database format.

    Returns
    -------
    out_file: str
        A handle to the output file that was generated after import.

    """
    if out_ext == 'mat':
        out_file = '{}.mat'.format(out_file)
        scipy.io.savemat(out_file, out_dict)
    elif out_ext == 'npz':
        out_file = '{}.npz'.format(out_file)
        np.savez(out_file, **out_dict)
    return out_file


def __get_out_file(in_file, out_dir):
    """

    Get the path of the output file.

    Parameters
    ----------
    in_file: str
        Path to input file.
    out_dir: str
        Path to output directory.

    Returns
    -------
    file_no_ext: str
        The file name without extension.
    out_dir: str
        The path to the output directory.
    out_file: str
        The path of the output file.

    """
    if out_dir == '':
        out_dir = os.path.dirname(in_file)
    file_no_ext = os.path.splitext(in_file)[0].split(os.sep)[-1]
    if out_dir == '':
        out_dir = '.'
    out_file = '/'.join([out_dir, file_no_ext])
    return file_no_ext, out_dir, out_file


def __import_file(in_file, out_file, out_ext, force=False, deli='\t',
                  dec_mark='.', pad=0, colheadlines=1):

    import_status = None
    num_dat = None
    col_heads = None
    out_file_exists = os.path.isfile('{}.{}'.format(out_file, out_ext))

    if (not out_file_exists) or (force is True):
        try:
            num_dat, col_heads = __process_file(in_file, dec_mark, deli,
                                                pad=pad,
                                                colheadlines=colheadlines)
            import_status = True
        except (ValueError, AttributeError):
            import_status = False

    return num_dat, col_heads, import_status


def import_del(in_file, force=False, deli='\t', dec_mark='.', out_ext='npz',
               out_dir='', pad=0, colheadlines=1):
    """

    Import a delimited data file into Numpy or Matlab database format. The file
    must have at least two data columns that are separated by :code:`deli`.

    Parameters
    ----------
    in_file: str
        The file handle of the delimited file that is to be imported.
    force: bool, optional
        If :code:`True`, existing output files will be overwritten during
        import. Default is :code:`False`.
    deli: str, optional
        The delimiter used to separate data columns in the delimited file.
        Default is tab.
    dec_mark: str, optional
        The decimal mark of the data file. Default is dot.
    out_ext: str, optional
        The file extension (format) of the output file. Default is :code:`npz`
        for Numpy database format. Alternative is :code:`mat` for Matlab
        database format.
    out_dir: str, optional
        The absolute or relative path to the output directory. Default is the
        current working directory.
    pad: positive int
        The numbers of data columns to skip. For :code:`pad = n`, the first
        :code:`n` data columns will not be imported.
    colheadlines: int, optional
        The number of lines spanned by the column headers. If several lines are
        spanned, the lines will be merged to generate the column keys in the
        output dictionary.

    Returns
    -------
    out_file: str
        A handle to the output file that was generated during import.
    import_status: str
        The import status of :code:`in_file`. If :code:`True`, the file was
        successfully imported. If :code:`False`, file import was attempted and
        failed. If :code:`None`, file import was not attempted (most likely
        because an output file with the same name already exists).
    out_dict: dict
        The data that was imported from :code:`in_file`.

    """
    _, out_dir, out_file_no_ext = __get_out_file(in_file, out_dir)
    out_dict = None

    num_dat, col_heads, import_status = \
        __import_file(in_file, out_file_no_ext, out_ext, force=force, deli=deli,
                      dec_mark=dec_mark, pad=pad, colheadlines=colheadlines)

    if import_status is True:
        out_dict = __write_to_out_dict(num_dat, col_heads)
        out_file = __save_out_file(out_file_no_ext, out_dict, out_ext)
    else:
        out_file = None

    return out_file, import_status, out_dict


def __gen_acc_time(step_time, steps, outformat='npz'):
    """

    For files produced by PCS Instrument test rigs, generate a continuous time
    axis by combining all step times from all steps.

    """
    # get index of last data point of each step
    current_step_end = np.where(np.subtract(step_time[1:], step_time[0:-1]) < 0)
    step_end = np.append(current_step_end[0], [step_time.shape[0] - 1])

    # get index of first data point of each step
    step_start = np.append([0], [step_end[0:-1] + 1])

    # add empty steps for mapper steps
    step_start_with_other = []
    step_end_with_other = []
    idx = 0
    for step_type in steps:
        if step_type == 'data':
            step_start_with_other.append(step_start[idx])
            step_end_with_other.append(step_end[idx])
            idx += 1
        elif step_type == 'other':
            if step_start_with_other:
                step_start_with_other.append(step_end_with_other[-1])
                step_end_with_other.append(step_end_with_other[-1])
            else:
                step_start_with_other.append(0)
                step_end_with_other.append(0)

    # loop over steps and create continuous time axis
    time_accumulated_s = copy.copy(step_time)
    offset = 0

    for step in range(1, len(step_end)):
        offset += step_time[step_end[step - 1]]
        time_accumulated_s[step_start[step]:step_end[step] + 1] += offset

    # save data to dictionary
    if outformat == 'mat':
        sub_dict = {'time_accumulated_s': time_accumulated_s,
                    'step_start': [s + 1 for s in step_start_with_other],
                    'step_end': [s + 1 for s in step_end_with_other]}
    else:
        sub_dict = {'time_accumulated_s': time_accumulated_s,
                    'step_start': step_start_with_other,
                    'step_end': step_end_with_other}
    return sub_dict


def __post_process_image_data(out_dict):
    """

    For SLIM Mapper Analysis files produced by PCS Instrument test rigs,
    extract the (x, y) coordinate system, generate an (x, y) grid and map the
    film thickness data to the grid.

    """
    img_dat = {}

    # get (unique) x and y axis values and allocate film thickness matrix
    x_ax = out_dict['x']
    y_ax = out_dict['y']
    x_uniq = np.unique(x_ax)
    y_uniq = np.unique(y_ax)
    x_index = np.zeros(len(x_ax))
    y_index = np.zeros(len(y_ax))
    film = np.zeros((len(x_uniq), len(y_uniq))) * float('nan')

    # get unique rank index for each element in x and y
    for idx, rank_value in enumerate(sorted(x_uniq)):
        x_index[np.where(x_ax == rank_value)[0]] = idx
    for idx, rank_value in enumerate(sorted(y_uniq)):
        y_index[np.where(y_ax == rank_value)[0]] = idx

    # combine x and y indices in a list that can be used to index the film array
    arr_idx = [x_index.astype(int), y_index.astype(int)]

    # assign all measured film thickness values to film thickness matrix
    film[arr_idx] = out_dict['film'][:, 0]

    # create variables that simplify plotting of film thickness data
    img_dat['film_surf'] = film
    img_dat['x_set'] = np.asarray(list(x_uniq))[:, None]
    img_dat['y_set'] = np.asarray(list(y_uniq))[:, None]
    img_dat['x_grid'], img_dat['y_grid'] = \
        np.meshgrid(img_dat['x_set'], img_dat['y_set'], indexing='ij')

    return img_dat


def __get_pcs_steps(in_file):
    """

    Get a list indicating the type of step for each step in a PCS data file.

    Parameters
    ----------
    in_file: str
        Path to PCS file

    Returns
    -------
    steps: list
        A list of step types. for numeric data, the step type is 'data', for
        other step types 'other'

    """
    steps = []
    with open(in_file) as dat_file:
        for line in dat_file:
            if line.startswith('Step ') and ' started at ' in line:
                steps.append('data')
            if line.lower().startswith('step type	mapper	') or \
                line.lower().startswith('step type	zero_check	') or \
                line.lower().startswith('step type	film_zero	') or \
                line.lower().startswith('step type	heating	'):
                steps[-1] = 'other'
    return steps


def import_pcs(in_file, force=False, out_ext='npz', out_dir=''):
    """

    Import a delimited data file that was produced by an MTM or EHD2 test rig
    manufactured by PCS Instruments. The method calls the :code:`import_del`
    method to perform a basic import of a delimited text file, and generates
    additional output variables that simplify data analysis.

    Parameters
    ----------
    in_file: str
        The file handle of the delimited file that is to be imported.
    force: bool, optional
        If :code:`True`, existing output files will be overwritten during
        import. Default is :code:`False`.
    out_ext: str, optional
        The file extension (format) of the output file. Default is :code:`npz`
        for Numpy database format. Alternative is :code:`mat` for Matlab
        database format.
    out_dir: str, optional
        The absolute or relative path to the output directory. Default is the
        current working directory.

    Returns
    -------
    out_file: str
        A handle to the output file that was generated during import.
    import_status: str
        The import status of :code:`in_file`. If :code:`True`, the file was
        successfully imported. If :code:`False`, file import was attempted and
        failed. If :code:`None`, file import was not attempted (most likely
        because an output file with the same name already exists).
    out_dict: dict
        The data that was imported from :code:`in_file`.

    """
    _, out_dir, out_file_no_ext = __get_out_file(in_file, out_dir)
    out_dict = None
    out_file = None

    num_dat, col_heads, import_status = \
        __import_file(in_file, out_file_no_ext, out_ext, force=force, deli='\t',
                      dec_mark='.', pad=8)

    steps = __get_pcs_steps(in_file)

    if import_status is True:
        out_dict = __write_to_out_dict(num_dat, col_heads, pcs=True)

        try:
            if 'step_time_s' in out_dict:
                t_dict = \
                    __gen_acc_time(out_dict['step_time_s'].astype(float), steps,
                                   out_ext)
                out_dict = {**out_dict, **t_dict}
            out_dict = {**out_dict, **__post_process_image_data(out_dict)}
        except KeyError:
            pass
        except IndexError:
            out_dict = None
            import_status = False

        if import_status:
            out_file = __save_out_file(out_file_no_ext, out_dict, out_ext)

    return out_file, import_status, out_dict


def __print_import_stats(in_file, status):
    """

    Print the import status to the console.

    Parameters
    ----------
    in_file: str
        The file name of the file for which to print the status.
    status: bool or None
        The import status of :code:`in_file`.

    """
    if status is False:
        out_col = _Colors.FAIL
    elif status is True:
        out_col = _Colors.OKGREEN
    else:
        out_col = _Colors.WARNING

    out_str = '\t'.join([str(status), str(in_file)])
    __print_status(out_str, out_col)


def __parse_args():
    """

    Parse all parser arguments that are provided when the script is running in
    a terminal.

    Returns
    -------
    args: Namespace
        The parsed parser arguments.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--force', action="store_true", default=False,
                        help='overwrite existing database files during import')
    parser.add_argument('-e', '--extension', action="store", default='txt',
                        help='specify file extension. default is "txt"')
    parser.add_argument('-d', '--delimiter', action="store", default='\t',
                        help='specify column delimiter. default is tab (\\t)')
    parser.add_argument('-m', '--mark', action="store", default='.',
                        help='specify decimal mark for numeric data. default is'
                             ' dot (.)')
    parser.add_argument('-o', '--outformat', action="store", default='npz',
                        help='specify output database format. default is "npz"'
                             ' for numpy database. use "mat" for matlab '
                             ' database format.')
    parser.add_argument('-r', '--recursive', action="store_true", default=False,
                        help='recursively walk through all sub-directories of'
                             ' current working directory')
    parser.add_argument('-p', '--pcs', action="store_true", default=False,
                        help='indicate if files are pcs files.')
    parser.add_argument('-c', '--colheadlines', action="store", default='1',
                        help='number of lines spanned by the column headers')
    args = parser.parse_args()
    return args


def import_dir(in_dir, in_ext='txt', recursive=False, force=False, deli='\t',
               dec_mark='.', out_ext='npz', out_dir='', print_stat=False,
               pcs=False, colheadlines=1):
    """

    Import all delimited data files in a directory into Numpy or Matlab
    database format. Optionally, all data files in a directory and all its
    child directories can be imported. The method can be applied to regular
    delimited files as well as files generated by test rigs made by PCS
    Instruments. All files must have at least two data columns that are
    separated by :code:`deli`.

    Parameters
    ----------
    in_dir: str
        Path to directory for which to import all files with extension
        :code:`in_ext`. If :code:`recursive=True`, imports are performed for all
        files with extension :code:`in_ext` in the directory tree with parent
        :code:`in_dir`.
    in_ext: str, optional
        File extension of files to import (without dot). Default is :code:`txt`.
    recursive: bool, optional
        If :code:`True`, all files in :code:`in_dir` and all its child
        directories are imported. Default is :code:`False`.
    force: bool, optional
        If :code:`True`, existing output files will be overwritten during
        import. Default is :code:`False`.
    deli: str, optional
        The delimiter used to separate data columns in the delimited file.
        Default is tab.
    dec_mark: str, optional
        The decimal mark of the data file. Default is dot.
    out_ext: str, optional
        The file extension (format) of the output file. Default is :code:`npz`
        for Numpy database format. Alternative is :code:`mat` for Matlab
        database format.
    out_dir: str, optional
        The path to the output directory where output databases are stored after
        import. By default, files are stored in :code:`in_dir` if
        :code:`recursive=False`. If :code:`recursive=True`, files are stored in
        the respective child directories of :code:`in_dir` if :code:`out_dir`
        is not specified.
    print_stat: bool, optional
        If :code:`True`, the current import status is printed to the console.
        Default is :code:`False`.
    pcs: bool, optional
        If :code:`True`, the delimited files are treated like files that were
        generated using an MTM or EHD2 test rig manufactured by PCS Instruments.
    colheadlines: int, optional
        The number of lines spanned by the column headers. If several lines are
        spanned, the lines will be merged to generate the column keys in the
        output dictionary.

    Returns
    -------
    in_files: ls of strings
        The file handles of all files for which import was attempted.
    out_files: ls of strings
        The file handles of all output files that were generated during the
        import process.
    import_status: ls of bools
        The import status of each file in :code:`in_files`. If :code:`True`,
        the file was successfully imported. If :code:`False`, file import was
        attempted and failed. If :code:`None`, file import was not attempted
        (most likely because an output file with the same name already exists).

    """
    in_files = __get_file_handles(in_dir, in_ext, recursive)  # type: ls
    out_files = []
    import_status = []

    if print_stat:
        print('importing {} files'.format(len(in_files)))
        print('status\tfilename\n'
              '======\t========')

    for in_file in in_files:

        if pcs is False:
            out_file, status, _ = import_del(in_file, force=force, deli=deli,
                                             dec_mark=dec_mark, out_ext=out_ext,
                                             out_dir=out_dir,
                                             colheadlines=colheadlines)
        else:
            out_file, status, _ = import_pcs(in_file, force=force,
                                             out_ext=out_ext,
                                             out_dir=out_dir)
        out_files.append(out_file)
        import_status.append(status)

        if print_stat:
            __print_import_stats(in_file, status)

    return in_files, out_files, import_status


if __name__ == "__main__":
    # if the file is executed as a script, import all data files in the
    # current working directory based on the parser arguments provided.
    ARGS = __parse_args()
    import_dir(os.getcwd(), in_ext=ARGS.extension, recursive=ARGS.recursive,
               force=ARGS.force, deli=ARGS.delimiter, dec_mark=ARGS.mark,
               out_ext=ARGS.outformat, out_dir=os.getcwd(), print_stat=True,
               pcs=ARGS.pcs, colheadlines=int(ARGS.colheadlines))
