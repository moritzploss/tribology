"""

Short demonstration of how to use the `slim2thick_batch` function to
batch-process MTM/EHD SLIM mapper bitmap files in order to extract numeric film
thickness data.

"""

import glob
import tribology as tr
import matplotlib.pyplot as plt


def plot_results(dat, mtm_file):
    mtm_dat = tr.load_npz(mtm_file)
    var_x = 'time_accumulated_s'
    var_y = 'traction_coeff'
    var_yy = 'mean_thickness_nm'

    fig, ax_y = plt.subplots(dpi=300)
    ax_yy = ax_y.twinx()
    ax_y.plot(mtm_dat[var_x], mtm_dat[var_y], '-')
    ax_yy.plot(dat[var_x], dat[var_yy], '-.')

    ax_y.set_xlabel(var_x)
    ax_y.set_ylabel(var_y)
    ax_yy.set_ylabel(var_yy)

    plt.tight_layout()
    plt.show()


def batch_process_mapper():
    # import MTM output file and spacer layer calibration file into npz
    # database format, save file handles
    mtm_file = tr.import_pcs('demo-mtm-output.txt')[0]
    slim_calib_file = tr.import_pcs('demo-3D_SpacerCalibration.txt')[0]

    # save file handles to zero bitmap file and all other bitmaps
    bmp_dir = 'demo mapper images'
    zero_bmp = bmp_dir + '/demo-3D_ZERO.bmp'
    bitmaps = sorted(glob.glob(bmp_dir + '/*3D_Step_*.bmp', recursive=True))

    # define where to save output plots
    out_dir = 'demo plots'

    # reduce the automatically detected contact area radius by 20 % to
    # exclude data close to the edges of the contact area
    crop = 0.2

    # only evaluate every 5th pixel (for speed)
    skip = 5

    # in the vertical direction (top, bottom), limit the evaluated area to the
    # central third of the reduced contact area. don't reduce the evaluated
    # area in the horizontal direction (right, left).
    aperture = {
        'top': 0.667,
        'right': 0,         # not required, only specified for clarity
        'bottom': 0.667,
        'left': 0           # not required, only specified for clarity
    }

    # define a list of MTM data variables that should be extracted together with
    # the film thickness. the variable names need to match those in the
    # `mtm_file` database
    pcs_vars = [
        'time_accumulated_s',
    ]

    # evaluate the mean film thickness for all bitmaps. the thickness will be
    # given relative to the mean thickness of `zero_bmp` since relative=True.
    dat = tr.slim2thick_batch(
        bitmaps, zero_bmp, slim_calib_file, mtm_file,
        plot=True, plt_dir=out_dir, skip=skip, crop=crop,
        relative=True, aperture=aperture, pcs_vars=pcs_vars
    )

    # print film thickness data to console (unit nanometer)
    print(dat['film_thickness_nm'])

    # plot traction coefficient and film thickness as a function of time
    plot_results(dat, mtm_file)


if __name__ == "__main__":
    batch_process_mapper()
