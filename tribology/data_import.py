"""

The below methods can be used to import delimited data files into Numpy or
Matlab database format.

"""

import argparse
import glob
import os
import re

import numpy as np
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


def __print_status(message, status_color=_Colors.ENDC):
    """

    Print a color-coded message to the terminal.

    Parameters
    ----------
    message: str
        The message to print to the terminal.
    status_color:
        The color in which to print the message.

    Returns
    -------
        None

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


def __write_to_out_dict(num_data, column_headers):
    """

    Extract the data columns from the num_data array and write them to a
    dictionary.

    Parameters
    ----------
    num_data: ndarray
        The data extracted from the delimited file, stored in a single table.
    column_headers: list of strings
        The column headers corresponding to the columns in :code:`num_data`

    Returns
    -------
    output_dict: dict
        A dictionary containing all data that is to be saved to the output
        database. Keys are based on column headers, values are data columns of
        num_data.

    """
    output_dict = {'column_headers': column_headers}
    for idx, column in enumerate(column_headers):
        # if empty data columns are not padded with tabs
        if idx >= num_data.shape[1]:
            output_dict[column] = np.zeros(num_data.shape[1]) * float('nan')
        else:
            # if data is of numeric type
            if __is_floatable(num_data[0, idx]):
                output_dict[column] = num_data[:, idx].astype(float)[:, None]
            # if data is of other type (string)
            else:
                output_dict[column] = num_data[:, idx].astype(object)[:, None]
    return output_dict


def __process_header(prev_line):
    """

    Process the column headers by removing special characters and converting to
    Matlab-optimized data type.

    Parameters
    ----------
    prev_line: list of strings
        The column headers of the delimited file.

    Returns
    -------
    col_headers: ndarray (dtype = object)
        The re-formated column headers.

    """

    # replace non-alphanumeric characters and trailing underscores
    col_headers = [re.sub("\W+", '_', item.lower()).strip('_')
                   for item in prev_line]
    # convert data type for easy matlab export
    col_headers = np.asarray(col_headers, dtype='object')
    return col_headers


def __process_data(split_line, num_dat, max_len, num_data_tables):
    """

    Append a data line to the current data table. If the length of the current
    data table exceeds the maximum permitted data table length, save the current
    data table to a list of data tables and initialise a new one.

    Parameters
    ----------
    split_line: list
        The data that is to be appended to the table.
    num_dat: ndarray
        The current data table to which the last line of data was appended.
    max_len: positive int
        The maximum length of a data table.
    num_data_tables: list
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


def __process_file(in_file, dec_mark, deli, padding=0):
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
    padding: positive int
        Ignore the first :code:`n` leading columns in the delimited file, where
        :code:`n = padding`. For example, if padding = 8, the first 8 columns
        are ignored.

    Returns
    -------
    output_dict: dict
        A dictionary containing all data that is to be saved to the output
        database. Keys are based on column headers, values are data columns of
        num_data.

    """
    max_len = 1000
    num_dat = []
    col_headers = []
    num_data_tables = []
    prev_line = ''

    with open(in_file) as dat_file:
        for line in dat_file:
            split_line = line.replace(dec_mark, '.').split(deli)

            # get rid of trailing newline characters
            if split_line[-1] == '\n':
                split_line[-1] = ''

            # check if first character is not (digit or minus symbol (hyphen))
            # to identify non-data lines. skip non-data lines.
            if not (line[0].isdigit() or line[0] == '-') or len(split_line) < 1:
                if split_line != ['']:
                    prev_line = split_line
                continue

            # if line contains data, split line into data fields, fill empty
            # fields with 'nan'
            split_line[:] = (item or 'nan' for item in split_line)
            # if this is the first data-containing line...
            if not len(col_headers):
                # get the column headers
                col_headers = __process_header(prev_line)
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
    output_dict = __write_to_out_dict(num_dat, col_headers)

    return output_dict


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
    args = parser.parse_args()
    return args


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
    in_files: list of strings
        The file handles to all delimited files that are to be imported.

    """
    if not recursive:
        in_files = sorted(glob.glob('{}{}*.{}'.format(in_dir, os.sep, ext)))
    else:
        in_files = []
        dir_list = [x[0] + os.sep for x in os.walk(in_dir)]
        for directory in dir_list:
            in_files.extend(sorted(glob.glob('{}*.{}'.format(directory, ext))))
        in_files = [f.replace(in_dir, '').lstrip(os.sep) for f in in_files]
    return in_files


def __save_out_file(f_no_ext, output_dict, out_ext, out_dir):
    """

    Save the imported data to an output database, either in Numpy or Matlab
    format.

    Parameters
    ----------
    f_no_ext: str
        The import file name without file extension.
    output_dict: dict
        The output data stored in a dictionary where keys correspond to column
        headers, values correspond to data.
    out_ext: str
        The file extension (format) of the output file. Options are :code:`npz`
        for Numpy format and :code:`mat` for Matlab database format.
    out_dir: str
        The path to the directory in which to save the output file.

    Returns
    -------
    out_file: str
        A handle to the output file that was generated after import.

    """
    out_file = ''
    if out_ext == 'mat':
        out_file = '{}/{}.mat'.format(out_dir, f_no_ext)
        scipy.io.savemat(out_file, output_dict)
    elif out_ext == 'npz':
        out_file = '{}/{}.npz'.format(out_dir, f_no_ext)
        np.savez(out_file, **output_dict)
    return out_file


def import_txt(in_file, force=False, deli='\t', dec_mark='.', out_ext='npz',
               out_dir=''):
    """

    Import a delimited data file into Numpy or Matlab database format.

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
    retdat: bool, optional
        If :code:`True`, a dictionary containing all data is returned.

    Returns
    -------
    out_file: str
        A handle to the output file that was generated during import.
    import_status: str
        The import status of :code:`in_file`. If :code:`True`, the file was
        successfully imported. If :code:`False`, file import was attempted and
        failed. If :code:`None`, file import was not attempted (most likely
        because an output file with the same name already exists).
    output_dict: dict
        The data that was imported from :code:`in_file`.

    """
    if out_dir == '':
        out_dir = os.getcwd()
        file_no_ext = os.path.splitext(in_file)[0]
    else:
        file_no_ext = os.path.splitext(in_file)[0].split(os.sep)[-1]

    out_file = os.sep.join([out_dir, ".".join([file_no_ext, out_ext])])
    out_file_exists = (os.path.isfile(out_file))
    import_status = None
    out_file = None
    output_dict = None

    if (not out_file_exists) or (force is True):
        try:
            output_dict = __process_file(in_file, dec_mark, deli)
            out_file = __save_out_file(file_no_ext, output_dict, out_ext,
                                       out_dir)
            import_status = True
        except ValueError:
            import_status = False

    return out_file, import_status, output_dict


def __print_import_stats(in_file, status):
    """

    Print the import status to the console.

    Parameters
    ----------
    in_file: str
        The file name of the file for which to print the status.
    status: True, False, None
        The import status of :code:`in_file`.

    """
    if status is False:
        out_col = _Colors.FAIL
    elif status is True:
        out_col = _Colors.OKGREEN
    else:
        out_col = _Colors.WARNING

    out_str = ': '.join([str(status), str(in_file)])
    __print_status(out_str, out_col)


def import_dir(in_dir, in_ext='txt', recursive=False, force=False, deli='\t',
               dec_mark='.', out_ext='npz', out_dir='', print_stat=False):
    """

    Import all delimited data files in a directory into Numpy or Matlab
    database format. Optionally, all data files in a directory and all its
    child directories can be imported.

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


    Returns
    -------
    in_files: list of strings
        The file handles of all files for which import was attempted.
    out_files: list of strings
        The file handles of all output files that were generated during the
        import process.
    import_status: list of bools
        The import status of each file in :code:`in_files`. If :code:`True`,
        the file was successfully imported. If :code:`False`, file import was
        attempted and failed. If :code:`None`, file import was not attempted
        (most likely because an output file with the same name already exists).

    """
    in_files = __get_file_handles(in_dir, in_ext, recursive)
    out_files = []
    import_status = []

    if print_stat:
        print('importing {} files:'.format(len(in_files)))

    for in_file in in_files:

        out_file, status, _ = import_txt(in_file, force=force, deli=deli,
                                         dec_mark=dec_mark, out_ext=out_ext,
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
    import_dir(os.getcwd(), ARGS.extension, ARGS.recursive,
               ARGS.force, ARGS.delimiter, ARGS.mark,
               ARGS.outformat, os.getcwd(), True)
