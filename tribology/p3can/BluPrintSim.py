import copy
import os
import re
import time

from Constants import PrintOpts
from check_environment import check_lapack, check_python_environment
from import_user_input import get_user_input
from system_functions import print_it, get_p3can_version, in_d


class Sim:
    # TODO: add class attribute describtion
    """Simulation information"""

    def __init__(self, in_file, out_dir):

        self.ui, self.invalid_input = get_user_input(in_file)

        os.environ["AUTO_PRINT"] = str(self.ui['AUTO_PRINT'])

        print_it("p3can simulation started")
        print_it("this is p3can version {}".format(get_p3can_version()))

        print_it("checking python and fortran libraries")

        check_python_environment()
        check_lapack()

        print_it("setting up simulation")
        print_it("validating user input", PrintOpts.lvl1.value)

        for key, value in self.invalid_input.items():
            print_it(
                "parameter '{} = {}' seems to be invalid. "
                "trying to go on anyway.".format(key, value),
                PrintOpts.lvl2.value)

        self.simulation_type = self.ui['simulation_type']
        self.simulation_name = self.ui['simulation_name']
        self.auto_print = self.ui['auto_print']
        self.auto_plot = self.ui['auto_plot']
        self.auto_report = self.ui['auto_report']
        self.start_time = time.time()

        self.calculation_input = ["simulation_type"]
        self.required_user_input = None
        self.optional_user_input = None
        self.positive_vars = None
        self.res_dir = None
        self.odd_vars = None
        self.alphanumeric_vars = None
        self.parameter_id = None
        self.infl_db_file_hand = None

        self.mk_uniq_res_dir(self.ui['simulation_name'], out_dir)
        print_it("results folder created: {}".format(self.res_dir),
                 PrintOpts.lvl1.value)

        self.input_file_name = self.save_user_input()
        print_it(
            "saved calculation input to file: {}".format(self.input_file_name),
            PrintOpts.lvl1.value)

    def print_time_elapsed(self):
        """Prints time elapsed since start of program"""
        print_it(("time elapsed: " "%d s" % (time.time() - self.start_time)))

    def get_time_elapsed(self):
        """Gets time elapsed since start of program"""
        return time.time() - self.start_time

    def finished(self):
        """Prints total run time at end of program execution"""
        print_it("\np3can finished")
        print_it("you can find all calculation outputs in {}".format(
            self.res_dir.replace('..{}'.format(os.sep), '')))
        self.print_time_elapsed()

    def mk_uniq_res_dir(self, simulation_name, out_dir):
        """Makes unique results folder based on 'simulation_name' and current
        date and time"""
        out_path = os.sep.join([out_dir, 'results/'])
        if os.path.isdir(out_path) is False:
            os.mkdir(out_path)
        valid_file_name = "".join(
            [c for c in simulation_name if re.match(r'\w|-', c)])
        results_folder = out_path + valid_file_name + '--' + time.strftime(
            "%Y-%m-%d--%H-%M-%S")
        unique_results_folder = copy.copy(results_folder)
        counter = 1
        while os.path.isdir('/'.join([os.getcwd(), unique_results_folder])):
            unique_results_folder = results_folder + '-{}'.format(counter)
            counter += 1
        self.res_dir = unique_results_folder
        os.mkdir(self.res_dir)

    def save_user_input(self):
        """Creates text file in results folder that contains the user input
        required for the current simulation. variables that are not required
        for the simulation but were defined by the user are not saved to the
        file"""
        input_file_name = os.sep.join([self.res_dir, "simulation_input.txt"])
        with open(input_file_name, "a") as results_file:
            print_str = "# list of relevant user input.\n" \
                        "# variables that were not used for the calculation " \
                        "are not listed here.\n\n"
            results_file.write(print_str)
            for key, value in self.ui.items():
                if isinstance(value, str):
                    print_str = key + " = '" + value + "'\n"
                else:
                    print_str = key + " = {}\n".format(value)
                results_file.write(print_str)
        return input_file_name

    def mk_uniq_parameter_id(self):
        """Make a unique parameter id (basically one very long string
        containing all inputs) based on the user input. this is used to
        identify unique calculation inputs (if the very long strings are
        identical) and allows to reuse
        influence matrices"""
        self.parameter_id = '-'.join(
            [str(in_d('simulation_type', self.ui, 'x')),
             str(in_d('radial_clearance', self.ui, 'x')),
             str(in_d('radial_clearance', self.ui, 'x')),
             str(in_d('number_balls', self.ui, 'x')),
             str(in_d('number_rollers', self.ui, 'x')),
             str(in_d('number_pins', self.ui, 'x')),
             str(in_d('number_planets', self.ui, 'x')),
             str(in_d('e_cb1', self.ui, 'x')),
             str(in_d('ny_cb1', self.ui, 'x')),
             str(in_d('diameter_cb1', self.ui, 'x')),
             str(in_d('length_cb1', self.ui, 'x')),
             str(in_d('type_profile_cb1', self.ui, 'x')),
             str(in_d('path_profile_cb1', self.ui, 'x')),
             str(in_d('profile_radius_cb1', self.ui, 'x')),
             str(in_d('e_cb2', self.ui, 'x')),
             str(in_d('ny_cb2', self.ui, 'x')),
             str(in_d('diameter_cb2', self.ui, 'x')),
             str(in_d('type_profile_cb2', self.ui, 'x')),
             str(in_d('profile_radius_cb2', self.ui, 'x')),
             str(in_d('length_cb2', self.ui, 'x')),
             str(in_d('path_profile_cb2', self.ui, 'x')),
             str(in_d('e_cb3', self.ui, 'x')),
             str(in_d('ny_cb3', self.ui, 'x')),
             str(in_d('diameter_cb3', self.ui, 'x')),
             str(in_d('type_profile_cb3', self.ui, 'x')),
             str(in_d('length_cb3', self.ui, 'x')),
             str(in_d('profile_radius_cb3', self.ui, 'x')),
             str(in_d('path_profile_cb3', self.ui, 'x')),
             str(in_d('global_force', self.ui, 'x')),
             str(in_d('res_x', self.ui, 'x')),
             str(in_d('res_y', self.ui, 'x')),
             str(in_d('res_pol', self.ui, 'x'))
             ]).replace('/', 'Slash')
        self.ui.update({'parameter_id': self.parameter_id})
