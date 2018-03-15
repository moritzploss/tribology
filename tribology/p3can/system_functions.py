import math
import os
import sys

import scipy.io

from Constants import PrintOpts, SubDir


def in_d(key, dic, return_if_false=None):
    """Checks if key is in dic. returns corresponding value if key is found;
    returns None (default) or return_if_false (if specified) otherwise"""
    if key in dic:
        return dic[key]
    else:
        if return_if_false is None:
            return None
        else:
            return return_if_false


def exit_program(reason):
    """Exits program and raises exception with text 'reason'"""
    print("\np3can error:\n")
    raise SystemExit(reason)


def print_it(message, level=PrintOpts.lvl0.value):
    """Pre-pends string 'level' to message and prints combined string to
    terminal"""
    if os.environ["AUTO_PRINT"] == "True":
        print("{}{}".format(level, message))


def make_directory(base_dir, directory_name):
    """Makes sub-directory of base_dir with name directory_name. Checks if
    directory exists already, then creates it if
    necessary. Returns path of new directory."""
    dir_to_create = os.sep.join([base_dir, directory_name])
    if os.path.isdir(dir_to_create) is False:
        os.mkdir(dir_to_create)
    return dir_to_create


def print_progress(current, to_reach):
    """Print progress bar to terminal based on ratio of (current:to_reach)"""
    current += 1
    prog_bar_len = 25
    approximate_progress = int(math.floor(current / to_reach * prog_bar_len))
    progress_bar = approximate_progress * '#' + '' + \
                   (prog_bar_len - approximate_progress) * ' '
    if os.environ["AUTO_PRINT"] == "True":
        sys.stdout.write(
            "\r{}progress: |{}| ({} %){}".format(PrintOpts.lvl1.value,
                                                 progress_bar,
                                                 round(
                                                     current / to_reach * 100),
                                                 '' if current != to_reach
                                                 else '\n'))
        sys.stdout.flush()
    return current


def print_rel_prox_to_sol(resulting_force, normal_load, threshold_factor):
    """Print progress bar to indicate proximity of half space solution to
    acceptable solution"""
    prog_bar_len = 25
    rel_dist = abs((resulting_force - normal_load) / normal_load)
    to_reach = 1
    proximity = 1 - rel_dist
    try:
        approximate_progress = int(round(proximity / to_reach * prog_bar_len))
    except ValueError:
        approximate_progress = 1
    progress_bar = approximate_progress * '#' + '' + min(
        (prog_bar_len - approximate_progress), prog_bar_len - 1) * ' '
    if abs((resulting_force - normal_load) / normal_load) < abs(
            threshold_factor):
        append_str = '\n'
    else:
        append_str = ''
    if os.environ["AUTO_PRINT"] == "True":
        sys.stdout.write(
            "\r{}progress: |{}| ({} N away from solution){}".format(
                PrintOpts.lvl1.value, progress_bar,
                round(abs((resulting_force - normal_load)), 1), append_str))
        sys.stdout.flush()
    return proximity


def get_p3can_version(plot_level=PrintOpts.lvl0.value, raise_dir=0):
    """Reads program version number from build.txt file and returns it. Prints
    warning message to terminal if build.txt
    not found"""
    try:
        with open(('..{}'.format(os.sep)) * raise_dir + '..{}build.txt'.format(
                os.sep)) as build_file:
            return build_file.readline()
    except FileNotFoundError:
        print_it(
            "warning: could not find file 'build.txt'. continue running unknown"
            " p3can version.",
            plot_level)
        return "<unknown>"


def save_to_matlab(dat_dict, directory, file_name):
    """Saves all variables in dat_dict to one common matlab database with name
    file_name. the variable names in the matlab database correspond to the key
    in dat_dict; the variable values corresponds to the value in dat_dict"""
    scipy.io.savemat(os.sep.join(
        [make_directory(directory, SubDir.matlab.value), file_name]) + '.mat',
                     dat_dict)


def to_preci(x, p):
    """
    returns a string representation of x formatted with a precision of p

    Originally based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/
                            JavaScriptCore/kjs/number_object.cpp

    I found the below code here:
    http://randlet.com/blog/python-significant-figures-format/

    """

    x = float(x)

    if x == 0.:
        return "0." + "0" * (p - 1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x / tens)

    if n < math.pow(10, p - 1):
        e = e - 1
        tens = math.pow(10, e - p + 1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens - x):
        n = n + 1

    if n >= math.pow(10, p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p - 1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e + 1])
        if e + 1 < len(m):
            out.append(".")
            out.extend(m[e + 1:])
    else:
        out.append("0.")
        out.extend(["0"] * -(e + 1))
        out.append(m)

    return "".join(out)
