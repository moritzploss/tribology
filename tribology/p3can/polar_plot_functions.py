import math

import matplotlib.pyplot as plt

from Constants import CMaps, PltOpts
from cartesian_plot_functions import save_fig


def polplt_line(theta_ax, r_data, r_label, title, res_dir, sub_dir=None,
                fname='file_name', r_max=None):
    """Polar 2D line plot"""
    ax, fig = setup_polplt()
    ax.plot(theta_ax, r_data, 'k-', color='0', markersize='10')
    format_polplt(ax, r_max, r_label, title)
    return save_fig(res_dir, sub_dir, fname)


def polplt_line_line(theta_ax1, r_data1, theta_ax2, r_data2, r_label, title,
                     label1, label2, res_dir, sub_dir=None,
                     fname='file_name', r_max=None):
    """Polar 2D line plot with 2 lines"""
    ax, fig = setup_polplt()
    ax.plot(theta_ax1, r_data1, 'k-', label=label1)
    ax.plot(theta_ax2, r_data2, 'k--', label=label2)
    format_polplt(ax, r_max, r_label, title)
    return save_fig(res_dir, sub_dir, fname)


def polplt_scatt_line(theta_ax1, r_data1, theta_ax2, r_data2, r_label, title,
                      res_dir, sub_dir=None, fname='file_name',
                      r_max=None):
    """Polar 2D scatter and line plot"""
    ax, fig = setup_polplt()
    ax.scatter(theta_ax1, r_data1, s=PltOpts.msize.value, c=r_data1,
               cmap=CMaps.markers.value)
    ax.plot(theta_ax2, r_data2, 'k--', color='0')
    format_polplt(ax, r_max, r_label, title)
    return save_fig(res_dir, sub_dir, fname)


def setup_polplt():
    """Generic figure setup for 2D polar plot"""
    fig = plt.figure(facecolor='w')
    ax = plt.subplot(111, polar=True)
    ax.grid(True)
    return ax, fig


def format_polplt(ax, r_max, r_label, title):
    """Format 2D polar plot labels and figure layout"""
    ax.set_theta_zero_location("S")
    ax.set_rlabel_position(112.5)
    label_position = ax.get_rlabel_position()
    if r_max is None:
        ax.text(math.radians(label_position + 15), ax.get_rmax() / 2, r_label,
                rotation=22.5, ha='center', va='center')
    else:
        ax.text(math.radians(label_position + 15), r_max / 2, r_label,
                rotation=22.5, ha='center', va='center')
    plt.title(title, y=1.08)
    if r_max is not None:
        ax.set_rmax(r_max)
