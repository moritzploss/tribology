import gc
import os

import matplotlib.pyplot as plt
import numpy as np

from Constants import CMaps, PltOpts, SubDir, PrintOpts
from ContactBodies import ContactBody
from system_functions import make_directory, print_it, print_progress


def plt_profile_approx(res_dir, sub_dir):
    """If contact body profile was generated from txt file input, plot the
    approximated profile together with the profile loaded from the file"""
    for body in gc.get_objects():
        if isinstance(body, ContactBody):
            if body.type_profile == "File":
                fname = "profile-approximation-{}".format(body.name)
                plt_2d_scatt_line(body.file_x_axis, body.file_x_profile,
                                  body.x_axis, body.x_profile,
                                  body.x_label, body.z_label, fname, res_dir,
                                  sub_dir, fname)


def plt_surf_roughness(res_dir):
    """If contact body profile is rough, plot surface roughness"""
    for body in gc.get_objects():
        if isinstance(body, ContactBody):
            if body.roughness_mat is not None:
                fname = "surface-roughness-{}".format(body.name)
                print(fname)
                plt_3d(body.x_axis, body.y_axis, body.roughness_mat,
                       body.x_label, body.y_label,
                       'roughness height in mm',
                       'surface roughness {}'.format(body.name), res_dir,
                       sub_dir=SubDir.profiles.value, fname=fname)


def save_fig(res_dir, sub_dir, fname):
    """Save figure with name 'fname' to png format"""
    if sub_dir is None:
        fhandle = os.sep.join([res_dir, fname + '.png'])
    else:
        fhandle = os.sep.join(
            [make_directory(res_dir, sub_dir), fname + '.png'])
    plt.savefig(fhandle,
                bbox_inches='tight',
                dpi=PltOpts.res.value)
    plt.close()
    return fhandle


def plt_3d(x_ax, y_ax, z_data, x_label, y_label, z_label, title, res_dir,
           sub_dir=None, fname='file_name',
           z_lim=None, azim=None):
    """Generate 3D surface plot"""
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = Axes3D(fig)
    x_grid, y_grid = np.meshgrid(x_ax, y_ax, sparse=False, indexing='ij')
    ax.plot_surface(x_grid, y_grid, z_data,
                    rstride=1,
                    cstride=1,
                    cmap=CMaps.default.value,
                    linewidth=0.5,
                    antialiased=True,
                    edgecolors='k')
    ax.set_xlabel('\n' + x_label)
    ax.set_ylabel('\n' + y_label)
    ax.set_zlabel('\n\n' + z_label)
    if z_lim is not None:
        ax.set_zlim(z_lim)
    if azim is not None:
        ax.view_init(azim=azim)
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    plt.title(title, y=1.08)
    return save_fig(res_dir, sub_dir, fname)


def plt_2d(x_ax, y_data, x_label, y_label, title, res_dir, sub_dir=None,
           fname='file_name'):
    """Generate 2D scatter plot"""
    plt.figure()
    plt.scatter(x_ax, y_data,
                s=PltOpts.msize.value,
                c=y_data,
                cmap=CMaps.markers.value,
                edgecolors='k')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title, y=1.08)
    return save_fig(res_dir, sub_dir, fname)


def plt_contact(body1, body2, profile_type, res_dir, sub_dir=None, fname=None):
    """Generate contact plot, showing two un-deformed contact bodies"""
    if fname is None:
        fname = '-'.join([body2.name, body1.name, 'contact', profile_type])
    if profile_type == PltOpts.DD.value:
        plt.figure()
        plt.scatter(body2.x_axis, body2.x_profile,
                    s=PltOpts.msize.value,
                    c=body2.x_profile,
                    cmap=CMaps.markers.value,
                    edgecolors='k')
        plt.scatter(body1.x_axis, -body1.x_profile,
                    s=PltOpts.msize.value,
                    c=body1.x_profile,
                    cmap=CMaps.markers.value,
                    edgecolors='k')
        plt.xlabel(body2.x_label)
        plt.ylabel(body2.z_label)
    elif profile_type == PltOpts.DDD.value:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = Axes3D(fig)
        x_grid, y_grid = np.meshgrid(body2.x_axis, body2.y_axis, sparse=False,
                                     indexing='ij')
        ax.plot_surface(x_grid, y_grid, body2.profile,
                        rstride=1,
                        cstride=1,
                        cmap=CMaps.contacts.value,
                        linewidth=0.5,
                        antialiased=True)
        ax.plot_surface(x_grid, y_grid, -body1.profile,
                        rstride=1,
                        cstride=1,
                        cmap=CMaps.contacts.value,
                        linewidth=0.5,
                        antialiased=True)
        ax.set_xlabel('\n' + body2.x_label)
        ax.set_ylabel('\n' + body2.y_label)
        ax.set_zlabel('\n\n' + 'profile in mm')
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    plt.title('contact {} and {}'.format(body1.name, body2.name), y=1.08)
    return save_fig(res_dir, sub_dir, fname)


def plt_profile(body, profile_type, res_dir, sub_dir=None, fname=None):
    """Plot contact body profile in 2D or 3D"""
    if fname is None:
        fname = body.name + '-profile-' + profile_type
    # 2D
    if profile_type == PltOpts.DD.value:
        plt_2d(body.x_axis, body.x_profile, body.x_label, body.y_label,
               (body.name + '-profile'), res_dir,
               sub_dir, fname)
    # 3D
    elif profile_type == PltOpts.DDD.value:
        plt_3d(body.x_axis, body.y_axis, body.profile, body.x_label,
               body.y_label, 'profile in mm',
               (body.name + '-profile'), res_dir, sub_dir, fname)


def plt_2d_scatt_line(x_ax1, y_data1, x_ax2, y_data2, x_label, y_label, title,
                      res_dir, sub_dir=None,
                      fname='file_name'):
    """Plot 2D scatter and line plot"""
    plt.figure()
    plt.plot(x_ax1, y_data1, '--', color='0')
    plt.scatter(x_ax2, y_data2, s=PltOpts.msize.value, c=y_data2,
                cmap=CMaps.markers.value, edgecolors='k')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title, y=1.08)
    return save_fig(res_dir, sub_dir, fname)


def plt_2d_2y_ax(x_ax1, y_data1, x_ax2, y_data2, x_label, title, leg1, leg2,
                 res_dir, sub_dir=None, fname='file_name'):
    """Plot 2D line plot with two y-axis"""
    fig, ax1 = plt.subplots()
    plot1 = ax1.plot(x_ax1, y_data1, '-', color='0', label=leg1)
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(leg1)
    plt.title(title, y=1.08)
    ax2 = ax1.twinx()
    plot2 = ax2.plot(x_ax2, y_data2, '--', color='0', label=leg2)
    ax2.set_ylabel(leg2)
    plots = plot1 + plot2
    labels = [l.get_label() for l in plots]
    legend = ax1.legend(plots, labels, loc=4)
    legend.get_frame().set_facecolor('#FFFFFF')
    return save_fig(res_dir, sub_dir, fname)


def plt_energy_ring_on_ring(tribo_system, res_dir, sub_dir, fname):
    """Plot energy figures for ring-on-ring tribosystems. This is a legacy
    function and should be replaced to follow the more general plotting approach
    that other plot functions use"""
    fig, ax1 = plt.subplots()
    plot1 = ax1.plot(tribo_system.sun.x_axis, tribo_system.sun.e_akin, '-',
                     color='0.75', label='e_a,kin sun')
    plot2 = ax1.plot(tribo_system.planet.x_axis, tribo_system.planet.e_akin,
                     '--', color='0', label='e_a,kin planet')
    ax1.set_xlabel(tribo_system.sun.x_label)
    ax1.set_ylabel('e_akin in W $\mathregular{m^{-2}}$')
    plt.title('e_akin vs pv_rel', y=1.08)
    ax2 = ax1.twinx()
    plot3 = ax2.plot(tribo_system.sun.x_axis, tribo_system.pv, '-.',
                     color='0.75', label='pv_rel sun/planet')
    ax2.set_ylabel('pv_rel in W $\mathregular{m^{-2}}$')
    plots = plot1 + plot2 + plot3
    labels = [l.get_label() for l in plots]
    legend = ax1.legend(plots, labels, loc=4)
    legend.get_frame().set_facecolor('#FFFFFF')
    return save_fig(res_dir, sub_dir, fname)


def plt_cum_dist(data, x_label, res_dir, sub_dir, fname):
    """Plot 2D cumulative distribution of data"""
    values, base = np.histogram(data, bins=200)
    cumulative = np.cumsum(values)
    plt.plot(base[:-1], cumulative / max(cumulative), c='k')
    plt.xlim(0, np.amax(data))
    plt.xlabel(x_label)
    plt.ylabel('cumulative distribution')
    return save_fig(res_dir, sub_dir, fname)


def plt_profiles(res_dir):
    """plot profiles of all contact bodies using python's garbage collector to
    find all objects of type ContactBody"""
    print_it("plotting profiles", PrintOpts.lvl1.value)
    num_cbs = list(
        isinstance(obj, ContactBody) for obj in gc.get_objects()).count(True)
    counter = 0
    for body in gc.get_objects():
        if isinstance(body, ContactBody):
            if body.type_profile == "File":
                fname = "profile-approximation-{}".format(body.name)
                plt_2d_scatt_line(body.file_x_axis, body.file_x_profile,
                                  body.x_axis, body.x_profile,
                                  body.x_label, body.z_label, fname, res_dir,
                                  SubDir.profiles.value, fname)

            if body.roughness_mat is not None:
                fname = "surface-roughness-{}".format(body.name)
                plt_3d(body.x_axis, body.y_axis, body.roughness_mat,
                       body.x_label, body.y_label,
                       'roughness height in mm',
                       'surface roughness {}'.format(body.name), res_dir,
                       sub_dir=SubDir.profiles.value, fname=fname)

            plt_profile(body, PltOpts.DD.value, res_dir, SubDir.profiles.value)
            plt_profile(body, PltOpts.DDD.value, res_dir, SubDir.profiles.value)
            counter = print_progress(counter, num_cbs)
