import os
from shutil import copyfile
import numpy as np
import itertools
import collections


def make_directory(dir_location, directory_name):
    """Make directory if it doesn't exist already"""
    if os.path.isdir(dir_location + directory_name) is False:
        os.mkdir(dir_location + directory_name)
    return os.sep.join([dir_location, directory_name])


def backup_current_user_input():
    """Make a copy of the USER_INPUT.py file that can be restored after all integration tests were run"""
    copyfile('{}/USER_INPUT.py'.format(os.getcwd()), '{}/USER_INPUT_BEFORE_TEST.py'.format(os.getcwd()))


def get_input_template(template):
    """"Load the current input template"""
    with open(template, 'r') as temp:
        return temp


def create_tests(res_dir, template, template_dir, value_pairs, keys):
    """Create integration test files based on the user input template specified"""
    make_directory(res_dir, '')
    in_file = '{}/{}'.format(template_dir, template)
    out_files = []
    for values in value_pairs:
        out_file = res_dir + '/{}--'.format(template.strip('.py'))
        out_file += '.py'
        for idx, value in enumerate(values):
            out_file = out_file.strip('.py') + '--{}-{}'.format(keys[idx], value).replace('.', ',')
            out_file += '.py'
        out_files.append(out_file)

    for f_idx, values in enumerate(value_pairs):
        lines_replacements = []
        input_lines = []
        lines_to_trash = []
        out_file = out_files[f_idx]
        for idx, value in enumerate(values):

            with open(in_file, 'r') as input_file:
                for line in input_file:

                    if keys[idx] in line.strip():
                        lines_to_trash.append(line)
                        lines_replacements.append('{} = {}\n'.format(keys[idx], value))
                    if line not in input_lines:
                        input_lines.append(line)
        for x in lines_to_trash:
            input_lines.remove(x)

        with open(out_file, 'w') as output_file:
            for line in input_lines:
                output_file.write(line)

            output_file.write('# test parameters\n')
            for line in lines_replacements:
                output_file.write(line)
            output_file.close()
    print('created {} test templates in {}'.format(len(value_pairs), res_dir))


def get_test_permutations(parameters):
    """Generate permutations of all parameter options specified"""
    parameters = collections.OrderedDict(sorted(parameters.items()))
    value_pairs = list(itertools.product(*parameters.values()))
    keys = [*parameters.keys()]
    return value_pairs, keys


def main():
    """Orchestrate integration test generation based on 'parameters' specified below. For each permutation of parameters
    in the 'parameters' dictionary, one unique p3can input file (test template) will be generated in the 'run' directory
    """
    # test template output directory
    res_dir = '{}/run'.format(os.getcwd())
    # directory that contains p3can user input templates
    template_dir = '../../UserInputTemplates'
    # user input template that should be used as a template for the integration test generation
    template = 'Template05_PinOnDisk.py'
    # parameter values that should be used for the integration test generation
    parameters = {'global_force': np.linspace(5, 70, 14),
                  # 'radial_clearance': np.linspace(-0.5, 0.5, num=4),
                  # 'number_rollers':  np.arange(2, 28, 7),
                  # 'type_profile_cb1': ["'File'", "'ISO'"],
                  # 'type_profile_cb2': ["'File'", "'ISO'"],
                  # 'type_profile_cb3': ["'File'", "'ISO'"],
                  # 'rot_velocity1': [-300, 300]
                  }
    os.chdir(template_dir)
    template_dir = os.getcwd()
    value_pairs, keys = get_test_permutations(parameters)
    create_tests(res_dir, template, template_dir, value_pairs, keys)

if __name__ == "__main__":
    main()
