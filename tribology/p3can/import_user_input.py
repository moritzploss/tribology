import os
import re
import copy

import numpy as np

from Constants import VarDicts, VarProps, OddEvenCheck, VarType, SimType, \
    Profiles
from system_functions import exit_program


def str_to_bool(s):
    """convert string representation of boolean to boolean"""
    if s == 'True':
        return True
    elif s == 'False':
        return False
    else:
        return s


def import_user_input(in_file):
    """read file input file line by line and try to extract user input.
    return dictionary with user input"""
    user_input = {}
    invalid_input = {}
    try:
        with open(in_file) as user_input_file:
            for line in user_input_file:
                # if line not empty
                if line.strip():
                    # if line is comment, ignore
                    if line.strip()[0] == '#':
                        continue
                    else:
                        # remove trailing comments and whitespace
                        value_pair = line.split('#')[0].strip(' ')
                        # split at '=' and remove whitespace
                        key = value_pair.split('=')[0].strip(' ')
                        value = value_pair.split('=')[1].strip('\n').strip(' ')
                        # value = value.strip("'").strip('"').strip("'")
                        user_input.update({key: value})
                        if value.startswith("'") or value.startswith('"'):
                            value = value[1:-1]
                            if re.search('[a-zA-Z]', value):
                                user_input.update({key: value})
                        elif value == 'True':
                            user_input.update({key: True})
                            if key == 'auto_print':
                                user_input.update({'AUTO_PRINT': 'True'})
                        elif value == 'False':
                            user_input.update({key: False})
                            if key == 'auto_print':
                                user_input.update({'AUTO_PRINT': 'False'})
                        else:
                            try:
                                user_input.update({key: float(value)})
                                if int(user_input[key]) == user_input[key]:
                                    user_input[key] = int(value)
                            except ValueError:
                                invalid_input.update({key: value})
        return user_input, invalid_input
    except FileNotFoundError:
        exit_program("can't find file {}".format(in_file))


def check_base_required(user_input):
    """checks if all basic (mandatory) variables are defined"""
    required_inputs = copy.deepcopy(VarDicts.required_user_input.value)
    for key, value in required_inputs.items():
        if key not in user_input:
            exit_program(
                "variable '{}' not defined in input file".format(key))
    check_vars(user_input)
    return required_inputs


def check_all_required(user_input, required_inputs):
    """checks if all calculation-specific (mandatory) variables are defined"""
    for key, value in required_inputs.items():
        if key not in user_input:
            exit_program(
                "variable '{}' not defined in input file".format(key))
    required_inputs.update(get_sim_specific_input(user_input)[0])
    check_vars(user_input)


def check_base_optional(user_input):
    """checks if all basic optional variables are defined properly"""
    to_check = {}
    optional_inputs = VarDicts.optional_user_input.value
    for key, value in user_input.items():
        if key in optional_inputs:
            to_check[key] = value
    check_vars(to_check)
    return optional_inputs


def check_all_optional(user_input, optional_inputs):
    """checks if all calculation-specific optional variables are defined
    properly"""
    to_check = {}
    optional_inputs.update(get_sim_specific_input(user_input)[1])
    for key, value in user_input.items():
        if key in optional_inputs:
            to_check[key] = value
    check_vars(to_check)


def check_sign(key, value, sign_var):
    """checks if the sign of a numeric variable is as expected"""
    if (isinstance(value, int) and not isinstance(value, bool)) is False \
            and isinstance(value, float) is False:
        exit_program(
            "variable '{}' with value '{}' should not be of type {}"
            .format(key, value, type(value)))
    if np.sign(value) != sign_var:
        exit_program("variable '{}' with value '{}' should be {} 0"
                     .format(key, value, {1.0: '<', -1.0: '>'}[np.sign(value)]))


def check_alpha(key, value):
    """Checks if the string has only alphanumeric characters, except for '_'
    and '-' """
    try:
        if value[0].isalnum() is False or value[-1].isalnum() is False:
            if value[0] == '-' or value[-1] == '-':
                exit_program(
                    "variable '" + key +
                    "' should not start or end with a hyphen '-' ")
            else:
                exit_program(
                    "variable '{}' with value '{}' should only contain "
                    "alphanumeric characters and hyphens '-'"
                    .format(key, value))
    except TypeError:
        exit_program(
            "variable '{}' with value '{}' should be of type string. strings "
            "need to be enclosed in quotes.".
            format(key, value))
    for character in value:
        if character.isalnum() is False and character != '-':
            exit_program(
                "variable '{}' with value '{}' should only contain "
                "alphanumeric characters and hyphens '-'".
                format(key, value))


def check_odd_even(key, value, should_be_odd_or_even=OddEvenCheck.Undefined):
    """Checks if the variable is odd or even"""
    if should_be_odd_or_even.value == OddEvenCheck.Odd.value or \
                    should_be_odd_or_even.value == OddEvenCheck.Even.value:
        if value % 2 == OddEvenCheck.Odd.value and \
                        should_be_odd_or_even.value != OddEvenCheck.Odd.value:
            exit_program("variable '" + value + "' should be an even number")
        elif value % 2 == OddEvenCheck.Even.value and \
                        should_be_odd_or_even.value != OddEvenCheck.Even.value:
            exit_program("variable '{}' should be an odd number".format(key))


def check_bool(key, value, isbool=True):
    """Checks if the variable is boolean"""
    if isinstance(value, bool) is not True and isbool is True:
        exit_program(
            "variable '{}' with value '{}' should be of type boolean".format(
                key, value))
    elif isinstance(value, bool) is True and isbool is False:
        exit_program(
            "variable '{}' with value '{}' (boolean) should not "
            "be of type boolean".format(key, value))


def check_int(key, value):
    """Checks if the variable is int"""
    if isinstance(value, int) is not True:
        exit_program(
            "variable '{}' with value '{}' should be of type int".format(key,
                                                                         value))


def check_vars(var_dict):
    """"Orchestrates variable checks"""
    for key, value in var_dict.items():
        if key in VarProps.bool_vars.value:
            check_bool(key, value, True)
        else:
            check_bool(key, value, False)
            if key in VarProps.positive_vars.value:
                check_sign(key, value, 1)
            if key in VarProps.negative_vars.value:
                check_sign(key, value, -1)
            if key in VarProps.odd_vars.value:
                check_odd_even(key, value, OddEvenCheck.Odd)
            if key in VarProps.alphanumeric_vars.value:
                check_alpha(key, value)
            if key in VarProps.int_vars.value:
                check_int(key, value)


def check_exists(key, user_input, else_quit=False):
    """Checks if the variable was defined by user, i.e., if it's in the
    global variables"""
    if key not in user_input and else_quit is True:
        exit_program("variable '{}' not defined".format(key))
    else:
        return user_input[key]


def get_sim_specific_input(user_input):
    """returns additional simulation-specific required and optional user
    input variables"""
    simulation_type = user_input['simulation_type']
    required_user_input = {}
    optional_user_input = {}
    if simulation_type == SimType.cyl_rol_bearing.value:
        required_user_input = {'length_cb1': VarType.real_number,
                               'e_cb3': VarType.real_number,
                               'ny_cb3': VarType.real_number,
                               'diameter_cb3': VarType.real_number,
                               'type_profile_cb1': VarType.string,
                               'type_profile_cb2': VarType.string,
                               'type_profile_cb3': VarType.string,
                               'number_rollers': VarType.integer,
                               'diameter_cb2': VarType.real_number}
        optional_user_input = {'radial_clearance': VarType.real_number,
                               'res_pol': VarType.integer,
                               'rot_velocity2': VarType.real_number,
                               'rms_roughness_cb3': VarType.real_number}
    elif simulation_type == SimType.deep_gro_ball_bearing.value:
        pass
    elif simulation_type == SimType.cyl_rol_thrust_bearing.value:
        required_user_input = {'length_cb1': VarType.real_number,
                               'e_cb3': VarType.real_number,
                               'ny_cb3': VarType.real_number,
                               'mean_diameter': VarType.real_number,
                               'type_profile_cb1': VarType.string,
                               'number_rollers': VarType.integer}
        optional_user_input = {'rot_velocity2': VarType.real_number}
    elif simulation_type == SimType.ball_on_disk.value:
        required_user_input = {'number_balls': VarType.integer,
                               'sliding_diameter': VarType.real_number,
                               'diameter_cb2': VarType.real_number,
                               'rot_velocity2': VarType.real_number}
    elif simulation_type == SimType.pin_on_disk.value:
        check_vars({'diameter_cb1': user_input['diameter_cb1'],
                    'type_profile_cb1': user_input['type_profile_cb1']})
        if user_input[
            'type_profile_cb1'] is Profiles.Circle.value and check_exists(
                'profile_radius_cb1') is True and \
                        user_input['diameter_cb1'] != 2 * user_input[
                    'profile_radius_cb1']:
            required_user_input.update({'length_cb1': VarType.real_number})
        required_user_input.update({'number_pins': VarType.integer,
                                    'sliding_diameter': VarType.real_number})
    elif simulation_type == SimType.four_ball.value:
        required_user_input = {'diameter_cb2': VarType.real_number}
    elif simulation_type == SimType.ball_on_three_plates.value:
        optional_user_input = {'plate_angle': VarType.real_number}
    elif simulation_type == SimType.ring_on_ring.value:
        required_user_input = {'length_cb1': VarType.real_number,
                               'type_profile_cb1': VarType.string,
                               'type_profile_cb2': VarType.string,
                               'number_planets': VarType.integer,
                               'rot_velocity2': VarType.real_number}
    else:
        exit_program(
            "'simulation_type' option '{}' not defined".format(simulation_type))
    return required_user_input, optional_user_input


def validate_file_path(path_variable):
    """Validates file path"""
    if os.path.isfile(path_variable) is False:
        exit_program("cannot find file '" + eval(path_variable) + "'"
                     "\nplease specify a valid file path (including file "
                     "extension) in variable '" + path_variable + "'" +
                     "\nthe file path needs to be defined relative to "
                     "'the input file'")


def validate_contact_body(user_input):
    """Checks if the user's contact body definition is consistent"""
    bodies_to_validate = {SimType.cyl_rol_bearing.value: ['cb1', 'cb2', 'cb3'],
                          SimType.cyl_rol_thrust_bearing.value: ['cb1'],
                          SimType.pin_on_disk.value: ['cb1'],
                          SimType.ring_on_ring.value: ['cb1', 'cb2']}
    for contact_body in bodies_to_validate.get(user_input['simulation_type'],
                                               []):
        type_profile = "type_profile_{}".format(contact_body)
        if user_input[type_profile] == Profiles.NoProfile.value:
            pass
        elif user_input[type_profile] == Profiles.ISO.value:
            check_vars({"length_" + contact_body: user_input[
                "length_" + contact_body]})
        elif user_input[type_profile] == Profiles.Circle.value:
            check_vars(
                {"length_" + contact_body: user_input["length_" + contact_body],
                 "profile_radius_" + contact_body: user_input[
                     "profile_radius_" + contact_body]})
        elif user_input[type_profile] == Profiles.File.value:
            validate_file_path(user_input["path_profile_" + contact_body])
        else:
            exit_program("'type_profile' option '{}' is undefined".format(
                user_input[type_profile]))


def get_user_input(in_file):
    """Orchestrates the process of user input handling"""
    user_input, invalid_input = import_user_input(in_file)
    required_inputs = check_base_required(user_input)
    optional_inputs = check_base_optional(user_input)
    spec_required_inputs, spec_optional_inputs = get_sim_specific_input(
        user_input)
    required_inputs.update(spec_required_inputs)
    optional_inputs.update(spec_optional_inputs)
    check_all_required(user_input, required_inputs)
    check_all_optional(user_input, optional_inputs)
    validate_contact_body(user_input)
    return user_input, invalid_input
