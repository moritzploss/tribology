"""

Import delimited data files.

"""

import argparse
import glob
import os
import re

import numpy as np
import scipy.io


class __bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def __print_status(string, status_color=__bcolors.ENDC):
    print(status_color + string + __bcolors.ENDC)


def __is_floatable(num):
    """
    check if 'num' can be converted to float
    """
    try:
        float(num)
        return True
    except ValueError:
        return False


def __to_float(num):
    """
    try to convert 'num' to float, return 'num' if it's not possible
    """
    try:
        float(num)
        return float(num)
    except ValueError:
        return num


def __assemble_data_table(num_data_tables, max_num_data_length):
    # re-assemble complete data table from list of data tables
    num_data = np.zeros((
        (len(num_data_tables) - 1) * max_num_data_length +
        num_data_tables[-1].shape[0],
        num_data_tables[-1].shape[1]),
        dtype=object)
    for idx, data_table in enumerate(num_data_tables):
        # do this for all but the last data table
        if idx + 1 < len(num_data_tables):
            num_data[idx * max_num_data_length: (idx + 1) * max_num_data_length,
            :] = data_table
        # do this for the last data table
        else:
            num_data[idx * max_num_data_length:, :] = data_table
    return num_data


def __write_to_out_dict(num_data, column_headers):
    """

    extract the data columns from the num_data array and write them to a
    dictionary

    Parameters
    ----------
    num_data
    column_headers

    Returns
    -------

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


def __process_header(col_headers, split_line, prev_line):
    # replace non-alphanumeric characters and trailing underscores
    # in column headers
    col_headers[:] = (re.sub("\W+", '_', item.lower()).strip('_')
                      for item in prev_line)
    # convert data type for easy matlab export
    col_headers = np.asarray(col_headers, dtype='object')
    # write first line to data table. note that "to_float" is not
    # necessary, but it's much faster
    num_dat = np.asarray(
        [__to_float(item.rstrip('\n'))
         for item in split_line]).reshape((1, len(split_line)))

    return col_headers, num_dat


def __process_data(split_line, num_dat, max_len, num_data_tables):
    # if data table becomes large, make new data table and add old
    # table to table list (for speed)
    if num_dat.shape[0] == max_len:
        num_data_tables.append(num_dat)
        num_dat = np.asarray(
            [__to_float(item.rstrip('\n')) for item in
             split_line]).reshape(
            (1, len(split_line)))
    # else simply append to data table
    else:
        num_dat = np.append(num_dat, np.asarray(
            [__to_float(item.rstrip('\n')) for item in split_line])
                            .reshape((1, len(split_line))), axis=0)
    return num_dat


def __process_file(in_file, decimal_mark, deli, padding=0):
    """
    extract data from tab-separated text file and return dictionary containing
    all data. the 'padding' parameter allows to ignore leading columns
    (from the left), i.e., if padding = 8, the first 8 columns are ignored.
    """
    max_len = 1000
    num_dat = []
    col_headers = []
    num_data_tables = []
    prev_line = ''
    with open(in_file) as dat_file:
        for line in dat_file:
            split_line = line.replace(decimal_mark, '.').split(deli)
            # get rid of trailing newline characters
            if split_line[-1] == '\n':
                split_line[-1] = ''
            # check if first character is not (digit or minus symbol (hyphen))
            # to identify non-data lines
            if not (line[0].isdigit() or line[0] == '-') or len(split_line) < 1:
                if split_line != ['']:
                    prev_line = split_line
                continue
            # if line contains data, split line into data fields, fill empty
            # fields with 'nan'
            split_line[:] = (item or 'nan' for item in split_line)
            # if this is the first data-containing line...
            if not len(col_headers):
                col_headers, num_dat = __process_header(col_headers, split_line,
                                                        prev_line)
            else:
                num_dat = __process_data(split_line, num_dat, max_len,
                                         num_data_tables)

    num_data_tables.append(num_dat)
    num_dat = __assemble_data_table(num_data_tables, max_len)
    output_dict = __write_to_out_dict(num_dat, col_headers)

    return output_dict


def __parse_args():
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
    parser.add_argument('-r', '--recursive', action="store_true", default=True,
                        help='recursively walk through all sub-directories of'
                             ' current working directory')
    args = parser.parse_args()
    return args


def __get_file_handles(directory, ext, recursive=False):
    if not recursive:
        in_files = sorted(glob.glob('*.{}'.format(ext)))
    else:
        in_files = []
        dir_list = [x[0] + os.sep for x in os.walk(directory)]
        for dir in dir_list:
            in_files.extend(sorted(glob.glob('{}*.{}'.format(dir, ext))))
        in_files = [f.replace(directory, '').lstrip(os.sep) for f in in_files]
    return in_files


def __save_out_file(f_no_ext, output_dict, out_ext, out_dir):
    out_file = ''
    if out_ext == 'mat':
        out_file = '{}/{}.mat'.format(out_dir, f_no_ext)
        scipy.io.savemat(out_file, output_dict)
    elif out_ext == 'npz':
        out_file = '{}/{}.npz'.format(out_dir, f_no_ext)
        np.savez(out_file, **output_dict)
    return out_file


def import_txt(in_file, force=False, deli='\t', dec_mark='.', out_ext='npz',
               out_dir=os.getcwd()):
    """
    orchestrate data import
    """

    file_no_ext = os.path.splitext(in_file)[0]
    out_file = os.sep.join([out_dir, ".".join([file_no_ext, out_ext])])
    out_file_exists = (os.path.isfile(out_file))

    import_status = None
    out_file = None

    if (not out_file_exists) or (force is True):
        try:
            output_dict = __process_file(in_file, dec_mark, deli)
            out_file = __save_out_file(file_no_ext, output_dict, out_ext,
                                       out_dir)
            import_status = True
        except ValueError:
            import_status = False
    return out_file, import_status


def __print_import_stats(print_stat, in_file, status):
    if print_stat:
        out_str = ': '.join([str(status), str(in_file)])
        if status is False:
            out_col = __bcolors.FAIL
        elif status is True:
            out_col = __bcolors.OKGREEN
        else:
            out_col = __bcolors.WARNING
        __print_status(out_str, out_col)


def import_dir(in_dir, in_ext, recursive=False, force=False, deli='\t',
               dec_mark='.',  out_ext='npz', out_dir=os.getcwd(),
               print_stat=False):
    """


    Parameters
    ----------
    in_dir: str
        Path to directory (top-level directory in case of `recursive=True`)
    in_ext: str
        File extension of files to import (without dot).
    recursive: bool, optional
        Default `False`. If `True`, all files in `in_dir` and all its child
        directories are imported. Output files are saved to the same directory
        as source files irrespective of `out_dir` value.
    force
    deli
    dec_mark
    out_ext
    out_dir
    print_stat

    Returns
    -------

    """
    in_files = __get_file_handles(in_dir, in_ext, recursive)
    out_files = []
    import_status = []

    if print_stat:
        print('importing {} files:'.format(len(in_files)))

    for in_file in in_files:

        out_file, status = import_txt(in_file, force=force, deli=deli,
                                      dec_mark=dec_mark, out_ext=out_ext,
                                      out_dir=out_dir)
        out_files.append(out_file)
        import_status.append(status)
        __print_import_stats(print_stat, in_file, status)

    return out_files, import_status


if __name__ == "__main__":
    args = __parse_args()
    import_dir(os.getcwd(), args.extension, args.recursive, args.force,
               args.delimiter, args.mark, args.outformat, os.getcwd(), True)
