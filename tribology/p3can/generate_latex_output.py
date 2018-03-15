import getpass
import math
import os
import platform
import re
import shutil
import subprocess
import sys

import jinja2
import matplotlib as mpl
import numexpr
import numpy as np
import scipy as sc

from Constants import SubDir, MasterDir, PrintOpts, TexTempl
from system_functions import make_directory, print_it, get_p3can_version


def generate_latex_output(simulation, tribo_system, ui=None, res_dir=None):
    """Orchestrate the LaTeX report generation. Load templates, set-up jinja2
    environment, define jinja2 variables, write LaTeX code to output file and
    trigger call to pdflatex"""
    print_it('generating simulation report by calling pdflatex')
    tex_file_path = make_directory(res_dir, SubDir.latex_files.value)
    make_directory(res_dir, SubDir.tex_figs.value)

    tex_filename = 'p3can_simulation_report.tex'
    tex_file_handle = tex_file_path + '/' + tex_filename
    calc_spec_tex_filename = 'p3can_simulation_report_calc_specific.tex'
    calc_spec_tex_file_handle = tex_file_path + '/' + calc_spec_tex_filename

    simulation.latex_jinja_env = jinja2.Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=jinja2.FileSystemLoader(os.path.abspath('.'))
    )

    template_report = simulation.latex_jinja_env.get_template(
        '/'.join([MasterDir.latex_templates.value,
                  TexTempl.Generic.value]).replace(os.getcwd(), ''))
    print_it('generating calculation-specific latex output',
             PrintOpts.lvl1.value)
    tribo_system.generate_latex_output(calc_spec_tex_file_handle, simulation)
    print_it('generating report figures', PrintOpts.lvl1.value)
    tribo_system.generate_latex_figures(res_dir=res_dir)

    with open(calc_spec_tex_file_handle, 'r') as calc_specific_latex_file:
        calc_specific_output = calc_specific_latex_file.read()

    latex_variables = {'section1': 'Long Form',
                       'section2': 'Short Form',
                       'sim_name': simulation.simulation_name,
                       'p3can_version': get_p3can_version(PrintOpts.lvl1.value),
                       'user_name': getpass.getuser(),
                       'python_version': re.sub('_', '\_',
                                                re.sub('#', '\#', sys.version)),
                       'platform_version': re.sub('_', '\_',
                                                  re.sub('#', '\#',
                                                         platform.version())),
                       'mpl_version': mpl.__version__,
                       'numexpr_version': numexpr.__version__,
                       'np_version': np.__version__,
                       'sc_version': sc.__version__,
                       'machine': re.sub('_', '\_', platform.machine()),
                       'sys_platform': platform.system(),
                       'run_time': (
                       str(math.ceil(simulation.get_time_elapsed())) + ' s'),
                       'user_input': generate_user_input_table(ui),
                       'calc_specific_output': str(calc_specific_output)}
    with open(tex_file_handle, 'w') as f:
        f.write(template_report.render(latex_variables))
    latex_it(tex_filename, res_dir)


def latex_it(tex_file, res_dir):
    """Compile tex_file by calling pdflatex. Raise compiler error in case
    something goes wrong. Delete LaTeX files after successful compilation
    (except for pdf file)"""
    print_it('calling pdflatex two times', PrintOpts.lvl1.value)
    pdf_filename = 'p3can_simulation_report.pdf'
    os.chdir(str(res_dir + os.sep + 'latex-files'))
    try:
        for _ in range(2):
            cmd = ['pdflatex', '-interaction', 'nonstopmode', tex_file]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            out, err = proc.communicate()
            raise_pdflatex_compiler_error(proc, pdf_filename, out, err, cmd)

        move_to_directory = re.sub("latex-files", "", os.getcwd())
        os.rename((os.sep.join([os.getcwd(), 'p3can_simulation_report.pdf'])),
                  (os.sep.join(
                      [move_to_directory, 'p3can_simulation_report.pdf'])))

        print_it('deleting latex files', PrintOpts.lvl1.value)
        os.chdir('..{}'.format(os.sep))
        shutil.rmtree('latex-files')
        os.chdir('..{}'.format(os.sep))
        print_it('{}report written to: {}'
                 .format(PrintOpts.lvl1.value, os.sep.join(
            (os.getcwd().strip("Results{}".format(os.sep)), res_dir,
             'p3can_simulation_report.pdf'))))

    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print_it('no pdflatex installation found', PrintOpts.lvl1.value)
            print_it('proceeding without generating report',
                     PrintOpts.lvl1.value)
            print_it(
                'check out https://www.latex-project.org/ '
                'if you want to install latex',
                PrintOpts.lvl1.value)
        else:
            raise


def raise_pdflatex_compiler_error(proc, pdf_filename, out, err, cmd):
    """Print LaTeX compiler error message to console"""
    retcode = proc.returncode
    if not retcode == 0:
        os.unlink(pdf_filename)
        print('{}pdflatex encountered an error. here is what it says:\n'.format(
            PrintOpts.lvl1.value),
              out.decode('ascii'), err.decode('ascii'))
        raise ValueError(
            'Error {} executing command: {}'.format(retcode, ' '.join(cmd)))


def generate_user_input_table(ui):
    """Generate user input table in a format that allows for copy-paste from
    LaTeX pdf to python editor"""
    calculation_input = []
    for key, value in ui.items():
        if key != 'parameter_id':
            calculation_input.append(key)
    calculation_input = sorted(calculation_input)
    table_user_input = []
    for i in calculation_input:
        if isinstance(ui[i], str):
            table_user_input.append((i, ("'" + str(ui[i]) + "'")))
        else:
            table_user_input.append((i, str(ui[i])))
    return table_user_input


def get_calc_specific_latex_template(latex_template, simulation):
    """Return the calculation-specific LaTeX template"""
    return simulation.latex_jinja_env.get_template(
        '/'.join([MasterDir.latex_templates.value, latex_template])
            .replace(os.getcwd(), ''))
