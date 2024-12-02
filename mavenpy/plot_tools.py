from collections.abc import Iterable

import matplotlib.dates as mdates
from dateutil.parser import parse


def get_ax(ax, index):

    '''Get an axis from a list or singleton of pyplot
    axes (e.g. from plt.subplots() or plt.subplots(nrows=2)).'''

    if isinstance(ax, Iterable):
        ax = ax[index]
    else:
        ax = ax
    return ax


def iter_ax(ax):

    '''Returns an iterable list of axes, given
    an axis that may be a singleton or tuple
    (e.g. from plt.subplots() or plt.subplots(nrows=2)).'''

    if isinstance(ax, Iterable):
        ax = ax
    else:
        ax = (ax,)
    return ax


def add_colorbar_outside(im, fig, ax, colorbar_width=0.01,
                         margin=0.01, **kwargs):

    '''Makes a colorbar that is located just outside of an existing axis
    (useful for multipanel comparison)'''

    # Based on StackOverflow post here:
    # https://stackoverflow.com/questions/71500468/positioning-multiple-colorbars-outside-of-subplots-matplotlib

    # Get bounding box that marks the boundaries of the axis:
    # [x0 (left), y0 (bottom), x1 (right), y1 (top)] of the axis.
    bbox = ax.get_position()

    # [left most position, bottom position, width, height] of color bar.
    cax = fig.add_axes(
        [bbox.x1 + margin, bbox.y0, colorbar_width, bbox.height])

    # Add colorbar using the new axis:
    cbar = fig.colorbar(im, cax=cax, **kwargs)

    return cbar


def format_xaxis(ax, start_date, end_date):

    start_dt = parse(start_date)
    end_dt = parse(end_date)
    n_days = (end_dt - start_dt).days

    if n_days > 2:
        minor_locator = mdates.DayLocator()
        major_locator = mdates.DayLocator(interval=1)
        major_formatter = mdates.DateFormatter("%b %-d")

    else:
        minor_locator = mdates.MinuteLocator(byminute=[0, 15, 30, 45])
        major_locator = mdates.HourLocator(interval=3)
        major_formatter = mdates.DateFormatter("%b %-d\n%H:%M")

    ax.xaxis.set_major_formatter(major_formatter)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.xaxis.set_major_locator(major_locator)

    # ax.set_xlabel("Time")
    # fig.autofmt_xdate()
    ax.set_xlim(start_dt, end_dt + dt.timedelta(days=1))
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
