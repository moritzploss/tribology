import copy
import math
import warnings

import numpy as np

from Constants import Profiles, PrintOpts
from system_functions import exit_program, print_it


def create_2d_profile2(diameter, length, type_profile, x_axis, res_x,
                       path_profile=None, circle_radius=None):
    """creates a 2d profile along x-axis (body length) based on user input"""
    x_profile = np.zeros(res_x)
    file_x_axis = None
    file_x_profile = None
    if type_profile == Profiles.ISO.value:
        # problem is that ISO profile is not defined at the end points of the
        # roller, hence the profile at the end points is calculated based on a
        # slightly shorter x_axis here.
        x_axis_slightly_shorter = copy.deepcopy(x_axis)
        x_axis_slightly_shorter[[0, -1]] = np.multiply([x_axis[0], x_axis[-1]],
                                                       0.99999999)
        x_profile = -0.00035 * diameter * np.log(
            1 / (1 - np.power(((2 * x_axis_slightly_shorter) / length), 2)))
        # x_profile[0] = x_profile[1] * 1.5
        # x_profile[-1] = x_profile[-2] * 1.5

    elif type_profile == Profiles.Circle.value:
        x_profile = np.array(abs(circle_radius) - np.sqrt(
            np.power(circle_radius, 2) - np.power(x_axis, 2)))
        sign_circle_radius = np.sign(circle_radius)
        x_profile = -np.multiply(x_profile, sign_circle_radius)

    elif type_profile == Profiles.File.value:
        content_text_file = []
        try:
            content_text_file = np.loadtxt(path_profile, delimiter="\t")
        except ValueError:
            exit_program(
                "\nexpecting the file '{}' to be a two column tab-delimited "
                ".txt-file"
                "\nwith first column corresponding to length of body and "
                "second column corresponding to "
                "profile data".format(path_profile))

        # get original x-axis and profile data from file and normalize it
        x_axis = np.array(content_text_file[:, 0]) - max(
            np.array(content_text_file[:, 0])) / 2
        file_x_axis = x_axis
        x_profile = np.array(content_text_file[:, 1]) - min(
            content_text_file[:, 1])
        file_x_profile = (x_profile - max(x_profile)) / 1000
        x_profile = x_profile - max(x_profile)
        step = max(1, round(len(content_text_file) / 15000))

        # deactivate filter warnings for polynomial fits
        warnings.simplefilter('ignore', np.RankWarning)

        try:
            print_it("fitting measured profile with 3rd grade polynomial",
                     PrintOpts.lvl1.value)
            start_point_left, start_point_right = fit_profile(x_axis, x_profile,
                                                              step, 3)
            x_axis_new = np.array(
                content_text_file[start_point_left:-start_point_right + 1, 0])
            x_profile_new = np.array(
                x_profile[start_point_left:-start_point_right + 1])
            x_axis_new = x_axis_new - min(x_axis_new)
            x_axis_new -= max(x_axis_new) / 2
            print_it("3rd grade polynomial fit successful",
                     PrintOpts.lvl1.value)
        except ValueError:
            print_it("3rd grade polynomial fit failed", PrintOpts.lvl1.value)
            print_it("fitting measured profile with 2nd grade polynomial",
                     PrintOpts.lvl1.value)
            start_point_left, start_point_right = fit_profile(x_axis, x_profile,
                                                              step, 2)
            x_axis_new = np.array(
                content_text_file[start_point_left:-start_point_right + 1, 0])
            x_profile_new = np.array(
                x_profile[start_point_left:-start_point_right + 1])
            x_axis_new = x_axis_new - min(x_axis_new)
            x_axis_new -= max(x_axis_new) / 2
            print_it("2nd grade polynomial fit successful",
                     PrintOpts.lvl1.value)
        warnings.simplefilter('always', np.RankWarning)

        # generate profile representation and corresponding x-axis based on
        # polynomial fit
        short_x_axis = np.linspace(x_axis_new[0], x_axis_new[-1], res_x)
        z = np.polyfit(x_axis_new, x_profile_new,
                       find_poly(x_axis_new, x_profile_new))
        p = np.poly1d(z)
        interpolated_x_profile = p(short_x_axis)
        x_axis = short_x_axis

        # normalize profile and convert units to millimeter
        x_profile = (interpolated_x_profile - min(interpolated_x_profile))
        x_profile = (x_profile - max(x_profile)) / 1000

    delta_x = abs(x_axis[0] - x_axis[1])
    return x_axis, x_profile, file_x_axis, file_x_profile, delta_x


def revolve_profile2(diameter, x_profile, y_axis, res_x, res_y):
    """creates a 3d profile by revolving self.x_profile around the central
    axis of the body"""
    profile = np.zeros((res_x, res_y))
    x_profile = copy.copy(-x_profile)
    sign_diameter = np.sign(diameter)
    for i in range(res_x):
        for j in range(res_y):
            bracket = math.pow((diameter / 2 - sign_diameter * x_profile[i]),
                               2) - math.pow(y_axis[j], 2)
            if bracket <= 0:
                profile[i, j] = abs(diameter) / 2
            else:
                profile[i, j] = abs(diameter) / 2 - sign_diameter * math.sqrt(
                    bracket)
    if diameter > 0:
        min_body_profile = profile.min()
    else:
        min_body_profile = profile[round(res_x / 2) - 1, round(res_y / 2) - 1]
    profile = -(profile - min_body_profile)
    y_profile = profile[:, math.floor(res_y / 2)]
    return profile, y_profile


def find_poly(x_axis_new, x_profile_new):
    """find the best polynomial fit, either degree 2, 3 or 4"""
    pol_grade = [2, 3, 4]
    z = (np.polyfit(x_axis_new, x_profile_new, item) for item in pol_grade)
    p = list((np.poly1d(item) for item in z))
    profiles = (item(x_axis_new) for item in p)
    diffs = list(sum(abs(item - x_profile_new)) for item in profiles)
    return pol_grade[diffs.index(min(diffs))]


def fit_profile(x_axis, x_profile, step, polyfit_grade):
    """ This function is mainly concerned with correcting the profile to get
    rid of the parts of the profile that do not belong to the raceway, but were
    measured anyway. To see what the function does, remove the comments below;
    you will see a matplotlib plot that is continuously updated while the
    function is running"""
    start_point = -step
    length = len(x_axis)
    half_length = round(length / 2)
    trigger = 0
    diff_triangle = 0
    diff_polynomial = 1
    start_point_left = None
    start_point_right = None
    while diff_polynomial > np.sum(diff_triangle):
        start_point += step
        seventh_length = round((half_length - start_point) / 3.5) + start_point

        # polynomial coefficients for a parabolic fit through the points
        # start_point + seventh_length + half_length
        z_parabolic = np.polyfit(
            x_axis[[start_point, seventh_length - 1, half_length - 1]],
            x_profile[[start_point, seventh_length - 1, half_length - 1]],
            polyfit_grade)
        # polynomial coefficients for a linear fit through the points
        # start_point + seventh_length
        z_triangle1 = np.polyfit(x_axis[[start_point, seventh_length - 1]],
                                 x_profile[[start_point, seventh_length - 1]],
                                 1)
        # polynomial coefficients for a linear fit through the points
        # start_point + half_length
        z_triangle2 = np.polyfit(
            x_axis[[seventh_length, 2 * seventh_length - 1]],
            x_profile[[seventh_length, 2 * seventh_length - 1]], 1)
        z_triangle3 = np.polyfit(x_axis[[2 * seventh_length, half_length - 1]],
                                 x_profile[
                                     [2 * seventh_length, half_length - 1]], 1)

        # evaluate the polynomials
        p_parabolic = np.poly1d(z_parabolic)
        [p_triangle1, p_triangle2, p_triangle3] = (np.poly1d(z) for z in
                                                   [z_triangle1, z_triangle2,
                                                    z_triangle3])

        # create shorter x_axis segments
        short_x_axis = x_axis[start_point:half_length]
        split_short_x_axis_1 = x_axis[start_point:seventh_length - 1]
        split_short_x_axis_2 = x_axis[seventh_length:2 * seventh_length - 1]
        split_short_x_axis_3 = x_axis[2 * seventh_length + 1:half_length - 1]

        # create interpolation data
        fitted_profile = p_parabolic(short_x_axis)
        triangle1 = p_triangle1(split_short_x_axis_1)
        triangle2 = p_triangle2(split_short_x_axis_2)
        triangle3 = p_triangle3(split_short_x_axis_3)

        # calculate difference between original profile and parabolic
        # interpolation
        diff_polynomial = np.sum(
            abs(fitted_profile - x_profile[start_point:half_length]))

        # calculate difference between original profile and triangular
        # interpolation
        combined_triang = np.concatenate((triangle1, triangle2, triangle3),
                                         axis=0)
        if len(combined_triang) > len(
                x_profile[start_point:start_point + len(combined_triang)]):
            diff_triangle = np.absolute(
                combined_triang -
                x_profile[start_point:start_point + len(combined_triang) + 1]
            )
        elif len(combined_triang) < len(
                x_profile[start_point:start_point + len(combined_triang)]
        ):
            diff_triangle = np.absolute(
                combined_triang -
                x_profile[start_point:start_point + len(combined_triang) - 1])
        else:
            diff_triangle = np.absolute(
                combined_triang -
                x_profile[start_point:start_point + len(combined_triang)])

        # uncomment these lines to see what happens during the fitting process:
        # print("diff_poly: " + str(diff_polynomial))
        # print("diff_triangles: " + str(diff_triangle))
        """
        import matplotlib.pyplot as plt
        plt.ion()
        if start_point == 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            line1, = ax.plot(x_axis, x_profile, 'k-', 
                             label='original profile data')
            line2, = plt.plot(short_x_axis, fitted_profile, '-.', 
                              label='poly 3rd grade approx.')
            plt.setp(line2, color='r', linewidth=2.0)
            line3, = plt.plot(split_short_x_axis_1, triangle1, '-', 
                              label='linear triangles approx.')
            plt.setp(line3, color='k', linewidth=3.0)
            line4, = plt.plot(split_short_x_axis_2, triangle2, '-')
            plt.setp(line4, color='k', linewidth=3.0)
            line5, = plt.plot(split_short_x_axis_3, triangle3, '-')
            plt.setp(line5, color='k', linewidth=3.0)
            dots, = plt.plot(x_axis[[start_point, seventh_length - 1, 
                                     half_length - 1]], 
                             x_profile[[start_point, seventh_length - 1, 
                                        half_length - 1]],
                             'o', label='interpolation nodes')
            plt.setp(dots, color='k', markersize=10)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc=4)
        else:
            line1.set_xdata(x_axis)
            line1.set_ydata(x_profile)
            line2.set_xdata(short_x_axis)
            line2.set_ydata(fitted_profile)
            line3.set_xdata(split_short_x_axis_1)
            line3.set_ydata(triangle1)
            line4.set_xdata(split_short_x_axis_2)
            line4.set_ydata(triangle2)
            line5.set_xdata(split_short_x_axis_3)
            line5.set_ydata(triangle3)
            dots.set_xdata(x_axis[[start_point, seventh_length - 1, 
                                   half_length - 1]])
            dots.set_ydata(x_profile[[start_point, seventh_length - 1, 
                                      half_length - 1]])
            import time
            time.sleep(0.02)
            fig.canvas.draw()
        """
        if diff_polynomial < np.sum(diff_triangle):
            if trigger == 1:
                start_point_right = start_point
                x_profile = np.flipud(x_profile)
                trigger = 2
            if trigger == 0:
                start_point_left = start_point
                x_profile = np.flipud(x_profile)
                diff_polynomial = 1
                diff_triangle = 0
                trigger = 1
                start_point = -step

    return start_point_left, start_point_right
