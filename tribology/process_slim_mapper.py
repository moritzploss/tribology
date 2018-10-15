# -*- coding: utf-8 -*-

"""

The below functions can be used to post-process SLIM (spacer layer imaging method)
mapper bitmap files created by MTM or EHD test rigs by PCS Instruments.

"""

import copy
import math
import os
import os.path
import sys

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import cv2


def __print_progress(current, to_reach):
    """

    Print progress bar to console based on ratio of current / to_reach.

    Parameters
    ----------
    current: int
        Current iteration number.
    to_reach: int
        Final iteration number to reach.

    Returns
    -------
    current: int
        The next iteration number.

    """
    current += 1
    prog_bar_len = 25
    approximate_progress = int(math.floor(current / to_reach * prog_bar_len))
    progress_bar = approximate_progress * '#' + '' + \
        (prog_bar_len - approximate_progress) * ' '
    sys.stdout.write("\rprogress: |{}| ({} %){}"
                     .format(progress_bar,
                             round(current / to_reach * 100),
                             '' if current != to_reach else '\n'))
    sys.stdout.flush()
    return current


def load_npz(file):
    """

    Load an npz database file into a mutable dictionary.

    Parameters
    ----------
    file: str
        Path to npz database file.

    Returns
    -------
    dat_cop: dict
        Loaded data.

    """
    dat = np.load(file)
    dat_cop = {}
    for var in dat.files:
        dat_cop[var] = copy.deepcopy(dat[var])
    return dat_cop


def __show_circles(img, x, y, rads, img_path):
    """

    Plot an image where the automatically detected circles are overlayed over
    the original bitmap file.

    Parameters
    ----------
    img: ndarray
        Array with image data.
    x: list
        Pixel coordinates of radii centers in x-direction
    y: list
        Pixel coordinates of radii centers in y-direction
    rads: list
        Radii of detected circles in units of pixels
    img_path:
        Path of original bitmap file

    """
    cv2.circle(img, (int(np.mean(x)), int(np.mean(y))), int(np.mean(rads)),
               (0, 255, 0), 2)
    # draw center of circle
    cv2.circle(img, (int(np.mean(x)), int(np.mean(y))), 2, (0, 255, 0), 3)
    cv2.imshow('detected circles ' + img_path.split(os.sep)[-1], img)
    cv2.waitKey(20000)
    cv2.destroyAllWindows()


def __find_mean_circles(img_gray):
    """

    Call the HoughCircles function using a range of fit parameters, then
    average over all detected circles to find mean center (x, y) and radius.

    Parameters
    ----------
    img_gray: ndarray
        Array containing grayscale image data.

    Returns
    -------
    x: list
        List of x-coordinates of radii centers
    y: list
        List of y-coordinates of radii centers
    rads: list
        list of radii

    """
    rads = []
    x = []
    y = []

    for par in np.linspace(1.1, 50, 40):
        circles = cv2.HoughCircles(img_gray, cv2.HOUGH_GRADIENT, par, 1500,
                                   param1=50, param2=30, minRadius=200,
                                   maxRadius=0)
        # pylint: disable=unsubscriptable-object
        new_circles = np.int16(np.around(circles))
        for circ in new_circles[0, :]:
            rads.append(circ[2])
            x.append(circ[0])
            y.append(circ[1])

    return x, y, rads


def __dist_2d(x_1, y_1, x_2, y_2):
    """

    Calculate the absolute distance between two points in the x-y plane.

    Parameters
    ----------
    x_1: float
        x-coordinate of first point.
    y_1: float
        y-coordinate of first point.
    x_2: float
        x-coordinate of second point.
    y_2: float
        y-coordinate of second point.

    Returns
    -------
    distance: float

    """
    return np.sqrt((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2)


def __plot_thick(file, thick, x_vals, y_vals, rgb, skip, crop, out_dir):
    """

    Plot the evaluated film thickness data in a 2D plot, using the extracted
    rgb values.

    Parameters
    ----------
    file: str
        Path to the bitmap file.
    thick: ndarray
        The thickness values for each pixel.
    x_vals: tuple
        Minimum and maximum evaluated pixel coordinates in x-direction.
    y_vals: tuple
        Minimum and maximum evaluated pixel coordinates in y-direction.
    rgb: ndarray
        Array of rgb tuples, containing color information for each pixel.
    skip: int
        The skip value provided by the user.
    crop: float
        The crop value provided by the user.
    out_dir: str
        Path to output directory.

    Returns
    -------
    out_file: str
        Path to plot file.

    """
    x_grid, y_grid = np.meshgrid(
        np.linspace(x_vals[0], x_vals[1], np.shape(thick)[0]),
        np.linspace(y_vals[0], y_vals[1], np.shape(thick)[1]),
        sparse=False, indexing='ij')

    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(111)
    ar_shape = x_grid.shape
    nu_shape = (ar_shape[0] * ar_shape[1], 1)
    x = x_grid.reshape(nu_shape)
    y = y_grid.reshape(nu_shape)
    z = np.fliplr(thick).reshape(nu_shape)
    idx = np.where(np.isnan(z) == False)
    ax.scatter(x[idx], y[idx],
               c=list(list(l) for l in (np.fliplr(rgb).reshape(nu_shape))[idx]),
               marker='s', s=0.5 * skip)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(file.split('/')[-1], fontsize=6)
    ax.set_aspect('equal', adjustable='box')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file = out_dir + '/' + file.rstrip('.bmp').split(os.sep)[-1] \
        + '-crop-{}-skip-{}.png'.format(crop, skip)
    plt.savefig(out_file)
    # plt.show()
    plt.close(fig)
    return out_file


def __adjust_rad_for_skip(diam, skip):
    """

    If necessary, increase the crop circle diameter (radius) for the chosen
    skip value. This is required because the diameter must be divisible by
    skip without rest.

    Parameters
    ----------
    diam: int
        Initial diameter of crop circle (in pixels).
    skip: int
        The numbers of pixels to skip when evaluating film thickness
        information.

    Returns
    -------
    diam: int
        Diameter of circle that is divisible by skip without rest.

    """
    while diam % skip != 0:
        diam += 1
    return diam


def __distance(c_1, c_2):
    """

    Get the absolute distance between two rgb color tuples in the 3D rgb space.

    Parameters
    ----------
    c_1: tuple
        First set of rgb values.
    c_2: tuple
        Second set of rgb values.

    Returns
    -------
    dist: float
        Absolute distance between color tuples.

    """
    (r_1, g_1, b_1) = c_1
    (r_2, g_2, b_2) = c_2
    dist = math.sqrt((r_1 - r_2) ** 2 + (g_1 - g_2) ** 2 + (b_1 - b_2) ** 2)
    return dist


def __rgb_to_thick(point, rgb_map):
    """

    Find rgb value that matches pixels best, then extract corresponding film
    thickness value.

    Parameters
    ----------
    point: tuple
        rgb value of point as extracted from bitmap file
    rgb_map: dict
        rgb map created from SLIM Spacer Calibration file, with rgb-tuples as
        keys and film thickness values as values.

    Returns
    -------
    film_thick: float
        Film thickness value for current pixel
    rgb_norm: tuple
        rgb value for pixel as extracted from color map, normalized to range 0-1
    dist: float
        Distance in RGB space between actual color and closest color

    """
    colors = list(rgb_map.keys())
    closest_colors = sorted(colors, key=lambda color: __distance(color, point))
    dist = __distance(colors[0], point)
    rgb = closest_colors[0]
    film_thick = rgb_map[rgb][0]
    rgb_norm = (c/255 for c in rgb)
    return film_thick, rgb_norm, dist


def __set_aperture(ap_in):
    """

    Make sure that aperture dictionary has correct format.

    Parameters
    ----------
    ap_in: dict
        Aperture as defined by user.

    Returns
    -------
    ap_out: dict
        Correctly formated and complete aperture dictionary

    """
    sides = ['top', 'right', 'bottom', 'left']
    ap_out = {}
    for side in sides:
        if ap_in and side in ap_in:
            ap_out[side] = ap_in[side]
        else:
            ap_out[side] = 0
    return ap_out


def __get_thick(file, x, y, rads, x_vals, y_vals, r_mean, rgb_map, skip=1,
                crop=0.25, aperture=None):
    """

    Get the film thickness values for each pixel of the bitmap image file.

    Parameters
    ----------
    file: str
        Path to bitmap image.
    x: list
        Centers of automatically detected circle centers in x-direction.
    y: list
        Centers of automatically detected circle centers in y-direction.
    rads: ndarray
        Radii of automatically detected circle centers.
    x_vals: tuple
        Minimum and maximum evaluated pixel values in x-direction.
    y_vals: tuple
        Minimum and maximum evaluated pixel values in y-direction.
    r_mean: int
        Radius of automatically detected circles.
    rgb_map: dict
        Dictionary mapping rgb tuples (keys) to thickness values (val).
    skip: int
        Skip value as defined by user.
    crop: float
        Crop value as defined by user.
    aperture: dict
        Aperture as defined by user.

    Returns
    -------
    thick: ndarray
        Thickness values for each pixel of bitmap image.
    rgb: ndarray
        Array of rgb tuples for each pixel of bitmap image.
    dist: ndarray
        Array of distance values between actual color and closest known color
        in the RGB space.

    """
    img = Image.open(file)
    pixels = img.load()

    diam = __adjust_rad_for_skip(2 * r_mean, skip)
    mat_size = int(diam / skip)
    thick = np.ones((mat_size, mat_size)) * float('nan')
    dist = np.ones((mat_size, mat_size)) * float('nan')
    rgb = np.zeros((mat_size, mat_size), dtype=object)
    rgb[:][:] = "white"
    x_mean = np.mean(x)
    y_mean = np.mean(y)

    aperture = __set_aperture(aperture)

    r_mean = np.mean(rads)
    for idx_1, x_idx in enumerate(range(x_vals[0], x_vals[1], skip)):
        for idx_2, y_idx in enumerate(range(y_vals[0], y_vals[1], skip)):
            position = __dist_2d(x_idx, y_idx, x_mean, y_mean)
            if position <= r_mean * (1 - crop) and \
                    y_idx <= y_mean + r_mean * (1 - aperture['bottom']) and \
                    y_idx >= y_mean - r_mean * (1 - aperture['top']) and \
                    x_idx <= x_mean + r_mean * (1 - aperture['right']) and \
                    x_idx >= x_mean - r_mean * (1 - aperture['left']):
                pix = pixels[x_idx, y_idx]
                thick[idx_1, idx_2], rgb[idx_1, idx_2], dist[idx_1, idx_2] = \
                    __rgb_to_thick(pix, rgb_map)

    return thick, rgb, dist


def __get_extrema(x, y, r_mean):
    """

    Get the minima and maxima of x and y circle center coordinates.

    Parameters
    ----------
    x: list
        x-coordinates
    y: list
        y-coordinates
    r_mean: int
        mean radius

    Returns
    -------
    x_vals: tuple
        Containing minimum and maximum x-value
    y_vals: tuple
        Containing minimum and maximum y-value

    """
    x_mean = int(np.mean(x))
    y_mean = int(np.mean(y))

    x_min = x_mean - r_mean
    x_max = x_mean + r_mean
    y_min = y_mean - r_mean
    y_max = y_mean + r_mean

    x_vals = (x_min, x_max)
    y_vals = (y_min, y_max)
    return x_vals, y_vals


def __load_grayscale_img(image):
    """

    Load an image file into an array.


    Parameters
    ----------
    image: str
        The path to a PIL-compatible image file in BMP format.

    Returns
    -------
    img_gray: ndarray
        Array containing grayscale image data.
    img_orig: ndarray
        Array containing original image data.

    """
    img = cv2.imread(image)
    img_orig = img.copy()
    img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    return img_gray, img_orig


def __load_rgb_map(file):
    """

    Load a npz database containing an imported PCS Spacer Calibration data set.

    To create an npz database from a PCS Spacer Calibration file, the
    `import_del` function of the `data_import` module should be used.

    Parameters
    ----------
    file: str
        Path to npz database.

    Returns
    -------
    rgb_dat: dict
        Dictionary mapping rgb tuples (keys) to thickness values (val).

    """
    spacer_calib = load_npz(file)
    rgb_dat = {}
    for idx, film in enumerate(spacer_calib['film']):
        rgb_dat.update(
            {(spacer_calib['r'][idx][0],
              spacer_calib['g'][idx][0],
              spacer_calib['b'][idx][0]): film})
    return rgb_dat


def __get_data_at_step(file, mtm_dat, var):
    """

    Find the data step that corresponds to an BPM image file and get the
    corresponding data (key `var`) from an MTM data file.

    Parameters
    ----------
    file: str
        The path to a BMP image file. The file name given by the PCS instrument
        may not be changed, otherwise the step extraction may fail.
    mtm_dat: dict
        Instrument output data loaded from npz database file.
    var: str
        A key in the npz database.

    Returns
    -------
    dat: ndarray
        An array containing a single data point (float, int)

    """
    if '_ZERO.bmp' in file:
        step = 1
    else:
        step = int(file.split('_')[-1].replace('.bmp', ''))
    step_idx = mtm_dat['step_start'][step-1]
    dat = mtm_dat[var][step_idx]
    return dat


def __calc_mean_thick(thick):
    """

    Calculate the mean value of the thickness data.

    Parameters
    ----------
    thick: ndarray
        Pixel-by-pixel thickness data.

    Returns
    -------
    Mean thickness: float

    """
    return np.nanmean(thick)


def slim2thick(file, rgb_map, rads=None, skip=1, crop=0.0, aperture=None):
    """

    Extract film thickness data from a PCS Instruments SLIM image file in BMP
    format. To obtain the thickness data, the color information of the BMP is
    converted to thickness information according to the rgb color map (Spacer
    Calibration file) provided by the instrument.

    The center point of the relevant (circular) contact area between the ball
    and the glass window is automatically detected. The radius of the contact
    area is automatically detected if `rads` is not provided. Thickness
    information outside of the detected area is discarded. The `crop` and
    `aperture` parameters can be used to further constrain the area of interest.

    The `skip` parameter determines the number of data points that is evaluated.
    For example, for `skip = 1`, every pixel within the area of interest is
    evaluated; for `skip = 10` every 10th pixel is evaluated.

    Parameters
    ----------
    file: str
        Path to a PCS bitmap 3D SLIM mapper output file
    rgb_map: str
        Path to a PCS Spacer Calibration file in npz format. The original
        instrument output file needs to be imported into npz database format
        using the `import_del` function of the `data_import` module first.
    rads: int, float, list, ndarray, optional
        The contact radius (in pixels) of the circular contact area between
        ball and glass window. If several radii are provided, the mean value
        will be used.
    skip: positive int, optional
        The number of pixels to skip (in both x and y direction). Higher skip
        numbers will lead to faster processing; lower skip numbers to higher
        accuracies.
    crop: float, optional
        Determines by how much the radius of the contact area is cropped
        during the data processing step. Value needs to be in the range [0, 1].
        A crop value of 0.5 will reduce the radius by 50 %.
    aperture: dict, optional
        A dictionary defining a rectangular area relative to the crop radius.
        Data points outside the rectangular area are not evaluated. The borders
        of the rectangle are defined relative to the crop radius. Values in the
        range [0, 1] may be defined for the following keys:

        - top
        - right
        - bottom
        - left

        If the aperture value is 0, the borders of the rectangle
        intersect with the crop radius on all four sides (i.e., the aperture
        has no effect). The following will reduce the crop radius on the top and
        bottom by 30 %:

        aperture = {'top': 0.3, 'bottom': 0.3}

    Returns
    -------
    thick: ndarray
        Film thickness in nanometer for each evaluated pixel.
    rgb: ndarray
        Best matching RGB tuple for each evaluated pixel obtained by comparison
        of bitmap image and SLIM Spacer Calibration file.
    rads: ndarray
        The automatically detected radii of the contact between ball and window.
        If `rads` is provided as an argument, `rads` is simply returned.
    x_vals: tuple
        Minimum and maximum evaluated pixel coordinate (in x-direction)
    y_vals: tuple
        Minimum and maximum evaluated pixel coordinate (in y-direction)
    dist: ndarray
        Array of distance values between actual color and closest known color
        in the RGB space, normalized to 255.

    """
    img_gray, _ = __load_grayscale_img(file)
    rgb_dat = __load_rgb_map(rgb_map)

    if not rads:
        x, y, rads = __find_mean_circles(img_gray)
    else:
        x, y, _ = __find_mean_circles(img_gray)

    # __show_circles(img_orig, x, y, rads, file)

    r_mean = int(np.mean(rads))
    x_vals, y_vals = __get_extrema(x, y, r_mean)
    thick, rgb, dist = __get_thick(file, x, y, rads, x_vals, y_vals, r_mean,
                                   rgb_dat, skip, crop, aperture)

    return thick, rgb, rads, x_vals, y_vals, dist


def slim2thick_batch(bitmaps, zero_bmp, rgb_map, mtm_file,
                     relative=False, pcs_vars=('time_accumulated_s',),
                     plot=False, plt_dir='', print_prog=True, skip=1,
                     crop=0.0, aperture=None):
    """

    Batch process a number of PCS Spacer Layer image files that share the same
    zero step, instrument output file and spacer layer calibration file.

    This function is a wrapper function for :code:`slim2thick`. See docstring
    for more information.

    Returns a dictionary containing mean film thickness values and, depending
    on user input, corresponding MTM variable values at the time of the film
    thickness measurement step.

    Parameters
    ----------
    bitmaps: tuple, list
        Paths to spacer layer bitmap files as list/tuple of strings.
    zero_bmp: str
        Path to bitmap image that corresponds to initial zero step of
        experiment.
    rgb_map: str
        Path to the shared spacer layer calibration file, imported into npz
        database format using the `import_del` function of the `import_data`
        module.
    mtm_file: str
        Path to the shared instrument output file, imported into npz database
        format using the `import_del` function of the `import_data` module.
    relative: bool, optional
        If True, the mean film thickness value is calculated relative to the
        mean film thickness of the zero step; if False, the absolute mean
        film thickness is calculated.
    pcs_vars: list, optional
        A list of MTM variables. Variable names must match those used in npz
        database created using the `import_del` function of the `import_data`
        module. By default, values for variable `time_accumulated_s` are
        included.
    plot: bool, optional
        If True, film thickness plots will be created for each image.
    plt_dir: str, optional
        Path to the output directroy in which to store plots.
    print_prog: bool
        If true, status updates are printed to the command line.
    skip: positive int, optional
        See function :code:`slim2thick`
    crop: float, optional
        See function :code:`slim2thick`
    aperture: dict, optional
        See function :code:`slim2thick`

    Returns
    -------
    outdict: dict
        The mean film thickness data for each bitmap file, along with extracted
        data from instrument output files (depending on user inputs).

    """
    mtm_dat = load_npz(mtm_file)

    # initialize dictionary that holds outputs
    out_dict = {
        'mean_thickness_nm': [],
        'mean_color_error': [],
        'skip': skip,
        'crop': crop,
        'aperture': aperture,
        'plots': []
    }
    for var in pcs_vars:
        out_dict.update({var: []})

    files = [zero_bmp] + list(bitmaps)
    zero_thick = 0
    rads = None

    if print_prog:
        print('analyzing {} images'.format(len(files)))

    # loop over images and extract film thickness for all steps
    for idx, file in enumerate(files):

        # calculate mean film thickness from image data
        if "_ZERO" in file:
            thick, rgb, rads, xtrem, ytrem, dist = slim2thick(
                zero_bmp, rgb_map, skip=skip, crop=crop, aperture=aperture)
            zero_thick = __calc_mean_thick(thick)
        else:
            thick, rgb, _, xtrem, ytrem, dist = slim2thick(
                file, rgb_map, rads=rads, skip=skip, crop=crop,
                aperture=aperture)
        mean_thick = __calc_mean_thick(thick)

        if relative:
            mean_thick -= zero_thick

        # store data in output dictionary
        out_dict['mean_thickness_nm'].append(mean_thick)
        out_dict['mean_color_error'].append(np.nanmean(dist))
        for var in pcs_vars:
            out_dict[var].append(__get_data_at_step(file, mtm_dat, var))
        if plot:
            out_dict['plots'].append(
                __plot_thick(file, thick, xtrem, ytrem, rgb, skip, crop,
                             plt_dir))
        if print_prog:
            __print_progress(idx, len(files))

    return out_dict
