"""

Import delimited .txt files.

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


def __process_file(file, decimal_mark, padding=0):
    """
    extract data from tab-separated text file and return dictionary containing
    all data. the 'padding' parameter allows to ignore leading columns
    (from the left), i.e., if padding = 8, the first 8 columns are ignored.
    """
    max_num_data_length = 1000
    num_data = []
    column_headers = []
    num_data_tables = []
    with open(file) as mtm_file:
        for line in mtm_file:
            # split line and ignore first 8 columns if the line is long (indicates 2D-type PCS file)
            if decimal_mark != '.':
                split_line = line.replace(decimal_mark, '.').split('\t')
            else:
                split_line = line.split('\t')
            # get rid of trailing newline characters
            if split_line and split_line[-1] == '\n':
                split_line[-1] = ''
            # check if first character is not (digit or minus symbol (hyphen)) to identify non-data lines
            if not (line[0].isdigit() or line[0] == '-') or len(split_line) < 1:
                if split_line != ['']:
                    previous_line = split_line
                continue
            # if line contains data, split line into data fields, fill empty fields with 'nan'
            split_line[:] = (item or 'nan' for item in split_line)
            # if this is the first data-containing line...
            if not len(column_headers):
                # replace non-alphanumeric characters and trailing underscores in column headers
                column_headers[:] = (re.sub("\W+", '_', item.lower()).strip('_')
                                     for item in previous_line)
                # convert data type for easy matlab export
                column_headers = np.asarray(column_headers, dtype='object')
                # write first line to data table. note that "to_float" is not necessary, but it's much faster
                num_data = np.asarray([__to_float(item.rstrip('\n')) for item in
                                       split_line]).reshape(
                    (1, len(split_line)))
            # ...else append current line to data table
            else:
                # if data table becomes large, make new data table and add old table to table list (for speed)
                if num_data.shape[0] == max_num_data_length:
                    num_data_tables.append(num_data)
                    num_data = np.asarray(
                        [__to_float(item.rstrip('\n')) for item in
                         split_line]).reshape(
                        (1, len(split_line)))
                # else simply append to data table
                else:
                    num_data = np.append(num_data, np.asarray(
                        [__to_float(item.rstrip('\n')) for item in split_line])
                                         .reshape((1, len(split_line))), axis=0)
    num_data_tables.append(num_data)
    # re-assemble complete data table from list of data tables
    num_data = np.zeros((
        (len(num_data_tables) - 1) * max_num_data_length +
        num_data_tables[-1].shape[0],
        num_data_tables[-1].shape[1]),
        dtype=object
    )
    for idx, data_table in enumerate(num_data_tables):
        # do this for all but the last data table
        if idx + 1 < len(num_data_tables):
            num_data[idx * max_num_data_length: (idx + 1) * max_num_data_length,
            :] = data_table
        # do this for the last data table (because it might be smaller than the others)
        else:
            num_data[idx * max_num_data_length:, :] = data_table

    # extract the data columns from the num_data array and write them to a dictionary
    output_dict = {'column_headers': column_headers}
    for idx, column in enumerate(column_headers):
        # take care of all other columns
        if column:
            # if empty data columns are not padded with tabs
            if idx >= num_data.shape[1]:
                output_dict[column] = np.zeros(num_data.shape[1]) * float('nan')
            else:
                # if data is of numeric type
                if __is_floatable(num_data[0, idx]):
                    output_dict[column] = num_data[:, idx].astype(float)[:,
                                          None]
                # if data is of other type (string)
                else:
                    output_dict[column] = num_data[:, idx].astype(object)[:,
                                          None]
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
    parser.add_argument('-r', '--recursive', action="store_true", default=False,
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
            in_files.extend(
                sorted(glob.glob('{}*.{}'.format(dir, ext))))
        in_files = [f.replace(directory, '').lstrip(os.sep) for f in in_files]
    return in_files


def __save_out_file(f_no_ext, output_dict, out_ext):
    if out_ext == 'mat':
        scipy.io.savemat('{}/{}.mat'.format(os.getcwd(), f_no_ext),
                         output_dict)
    elif out_ext == 'npz':
        np.savez('{}/{}'.format(os.getcwd(), f_no_ext), **output_dict)


def import_txt(in_files, force=False, deli='\t', dec_mark='.',
               out_ext='npz', recursive=False, status=False):
    """
    orchestrate data import
    """
    out_files = []
    imported = []

    for file in in_files:
        f_no_ext = os.path.splitext(file)[0]

        out_file_exists = (f_no_ext + '.{}'.format(out_ext)) \
                          in os.path.dirname(file)

        if out_file_exists and force is False:
            __print_status("file '{}.{}' already exists. did not overwrite {} "
                           "database".format(f_no_ext,
                                             args.outformat,
                                             args.outformat),
                           __bcolors.ENDC)
        else:
            try:
                output_dict = __process_file(file, dec_mark)
                __save_out_file(f_no_ext, output_dict, out_ext)

                print_string = "file '{}' imported. {} database".format(file,
                                                                        out_ext)
                if out_file_exists:
                    __print_status("{} overwritten".format(print_string),
                                   __bcolors.OKGREEN)
                else:
                    __print_status("{} created".format(print_string),
                                   __bcolors.OKGREEN)

            except (ValueError):
                __print_status(
                    "### did not import file '{}'. file does not seem to be a "
                    "simple tab-separated file ###"
                    .format(file), __bcolors.FAIL, )
    print('{}\nfinished all imports'.format(120 * '-'))
    return


if __name__ == "__main__":
    args = __parse_args()
    in_files = __get_file_handles(os.getcwd(), args.ext, args.recursive)
    import_txt(in_files, force=args.force, deli=args.delimiter,
               dec_mark=args.mark, out_ext=args.outformat)
