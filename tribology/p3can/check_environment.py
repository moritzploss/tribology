import sys

from Constants import PrintOpts
from system_functions import print_it, exit_program


def check_lapack():
    """Check if BLAS/LAPACK distribution is installed"""
    import numpy as np
    if not np.__config__.lapack_opt_info or not np.__config__.blas_opt_info:
        print_it('could not find BLAS/LAPACK install\n'
                 'install BLAS/LAPACK to improve p3can performance',
                 PrintOpts.lvl1.value)
        print_it('did not check for ATLAS', PrintOpts.lvl1.value)
    else:
        print_it('found BLAS/LAPACK install', PrintOpts.lvl1.value)
        print_it('did not check for ATLAS', PrintOpts.lvl1.value)


def check_python_environment():
    """Check if all required python environments are installed before actual
    calculation is started"""
    print_it('python: {}'.format(sys.version.replace('\n', ' ')),
             PrintOpts.lvl1.value)
    try:
        import matplotlib as mpl
        print_it('matplotlib: {}'.format(mpl.__version__),
                 PrintOpts.lvl1.value)
    except ImportError:
        print_import_error('matplotlib', PrintOpts.lvl1.value)
        exit_program('module not found')
    try:
        import numexpr
        print_it('numexpr: {}'.format(numexpr.__version__),
                 PrintOpts.lvl1.value)
    except ImportError:
        print_import_error('numexpr', PrintOpts.lvl1.value)
        exit_program('module not found')
    try:
        import numpy
        print_it('numpy: {}'.format(numpy.__version__),
                 PrintOpts.lvl1.value)
    except ImportError:
        print_import_error('numpy', PrintOpts.lvl1.value)
        exit_program('module not found')
    try:
        import scipy
        print_it('scipy: {}'.format(scipy.__version__),
                 PrintOpts.lvl1.value)
    except ImportError:
        print_import_error('scipy', PrintOpts.lvl1.value)
        exit_program('module not found')
    try:
        import jinja2
        print_it('jinja2: {}'.format(jinja2.__version__),
                 PrintOpts.lvl1.value)
    except ImportError:
        print_import_error('jinja2', PrintOpts.lvl1.value)
        print_it('this will cause an error if auto_report is True')


def print_import_error(module_name, level):
    """Print import error message"""
    print_it("could not find (import) module '{}'".format(module_name), level)
    print_it("the module is missing or unknown to the python interpreter",
             level)
