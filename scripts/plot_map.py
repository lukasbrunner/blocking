#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Copyright (C) 2018 Lukas Brunner (Wegener Center/University of Graz)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Abstract: A simple plot script for longitude-latitude resolved plots on a map.
TODO:
- add the option to draw a contour again

"""
import logging
import os
import argparse
import xarray
import numpy as np
from collections import OrderedDict

import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import blocking.utils as ut
mpl.rc('font', **{'size': 11})

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

PLOTPATH = './../plots'


def get_season(months, str_='{}'):
    """
    Creates a season string.

    Parameters:
    - months (list of int)
    - str_ (str, optional): Formatter string, should contain exactly one {}
      at the position where the season substring is included.

    Returns:
    str
    """
    if months is None:
        return ''
    elif len(set(months).difference([1, 2, 12])) == 0:
        return str_.format('DJF')
    elif len(set(months).difference([3, 4, 5])) == 0:
        return str_.format('MAM')
    elif len(set(months).difference([6, 7, 8])) == 0:
        return str_.format('JJA')
    elif len(set(months).difference([9, 10, 11])) == 0:
        return str_.format('SON')
    elif len(set(months).difference([11, 12, 1, 2, 3])) == 0:
        return str_.format('NDJFM')
    elif len(set(months).difference([5, 6, 7, 8, 9])) == 0:
        return str_.format('MJJAS')
    else:
        return str_.format('-'.join(map(str, months)))


def get_title(args, ds):
    """
    Creates a title or returns user set title.

    Parameters:
    - args (parser object)
    - ds (xarray dataset)

    Returns:
    str
    """
    if args.title is not None:
        return args.title

    title = args.varn
    if args.period is not None:
        title += ' {}'.format(
            ' to '.join(args.period))
    title += get_season(args.months, ' ({})')
    return title


def get_filename(args):
    """
    Creates a filename or returns user set filename.

    Parameters:
    - args (parser object)

    Returns:
    str
    """
    if args.plotname is not None:
        return args.plotname

    filename = 'Mapplot_{}'.format(args.varn)
    if args.period is not None:
        filename += '_{}'.format('-'.join(args.period))
    filename += get_season(args.months, '_{}') + '.png'
    return filename


def get_kwargs_from_config(config):
    """
    Sets default settings or imports settings from config file.

    Parameters:
    - config (str): path/to/config (no extension).

    Returns:
    config object
    """

    if config is None:
        class kwargs():
            subplots = dict(
                figsize = (8, 4),
                subplot_kw = dict(
                    projection =
                    # ccrs.Robinson(central_longitude=-90)
                    # ccrs.Mollweide()
                    ccrs.PlateCarree(central_longitude=-90)
                ))
            subplots_adjust = dict(
                left=.1,
                right=.87,
                bottom=0.1,
                top=.95)
            pcolormesh = dict()
    else:
        kwargs = __import__(config)

    return kwargs


def draw_polygon(ax, kwargs):
    """
    Draws user defined polygons on the map.

    Parameters:
    - ax (axis object)
    - kwargs.polygons (config object): List of polygons to draw.

    Returns:
    None
    """
    if getattr(kwargs, 'polygons', False):
        if not isinstance(kwargs.polygons, list):
            kwargs.polygons = list(kwargs.polygons)
        for poly in kwargs.polygons:
            ax.add_patch(poly)


def draw_dots(ax, ds, args, kwargs_config):
    """
    Draws symbols on the map.

    Parameters:
    - ax (axis object)
    - ds (xarray dataset)
    - args (parser object)
    - kwargs.dots (config object): ax.scatter kwargs.

    Returns:
    None
    """
    # default settings
    kwargs = dict(
        s = .1,
        color = 'black',
        marker = 'o',
        transform=ccrs.PlateCarree(),
        # edgecolor='none',  # suppresses strange behaviour for s<1
        # alpha = None
    )
    # update settings from config if applicable
    if getattr(kwargs_config, 'dots', False):
        kwargs.update(kwargs_config.dots)

    if args.indicator is not None:
        ds = ds.transpose(ut.get_latitude_name(ds), ut.get_longitude_name(ds))
        var = ds[args.indicator]
        [lats, lons] = [ds[dimn].data for dimn in var.dims]
        idxs = np.where(var.data != 0)

        if args.indicator_size is not None:
            if args.indicator_size > 0 :
                ms = args.indicator_size
            else:
                ms = var.data[idxs] / float(var.max()) * -args.indicator_size

            # update settings from parser if applicable
            kwargs.update(dict(s=ms))

        # plot symbol to all coordinates where var != 0
        sp = ax.scatter(lons[idxs[1]], lats[idxs[0]], **kwargs)

        # make labels
        if len(np.unique(kwargs['s'])) == 1:  # only one dot size
            sp = ax.scatter(999, 999, **kwargs)
            legend = ax.legend([sp], [args.indicator], scatterpoints=1,
                               loc='lower right', fontsize='xx-small')
        else:  # multiple dot sizes
            sp = OrderedDict()
            for pp in [25, 50, 75, 100]:  # percentiles
                kwargs.update(dict(s=np.percentile(ms, pp)))
                ll = '{}p: {:.2%}'.format(pp, np.percentile(var.data[idxs], pp))
                sp[ll] = ax.scatter(999, 999, **kwargs)
            # /make labels
            legend = ax.legend(sp.values(), sp.keys(),
                               title=args.indicator,
                               scatterpoints=1,
                               loc='lower right',
                               fontsize='xx-small')
        plt.setp(legend.get_title(), fontsize='x-small')  # adjust title size

def plot(ds, args):

    kwargs = get_kwargs_from_config(args.config)

    lon_name = ut.get_longitude_name(ds)
    lat_name = ut.get_latitude_name(ds)
    lons = range(-180, 180+1, 60)
    lats = range(-90, 90+1, 30)

    # create figure and draw map
    fig, ax = plt.subplots(**kwargs.subplots)
    fig.subplots_adjust(**kwargs.subplots_adjust)

    ax.set_global()
    ax.coastlines()
    ax.gridlines(xlocs=lons, ylocs=lats, **dict(linestyle = ':'))

    # create axis for colorbar and adjust it
    # posn = ax.get_position()
    # cbar_ax = fig.add_axes([posn.x0 + posn.width + .01, .14, .04, .77])

    pc = ds[args.varn].plot.pcolormesh(
        x=lon_name,
        y=lat_name,
        ax=ax,
        transform=ccrs.PlateCarree(),
        # cbar_ax=cbar_ax,
        **kwargs.pcolormesh)

    # set labels and ticks
    ax.set_xlabel('Longitude')
    try:  # ticks only work for rectangular coordinate systems
        ax.set_xticks(lons, crs=ccrs.PlateCarree())
    except RuntimeError:
        pass
    longitude_formatter = LongitudeFormatter()
    ax.xaxis.set_major_formatter(longitude_formatter)

    ax.set_ylabel('Latitude')
    try:
        ax.set_yticks(lats, crs=ccrs.PlateCarree())
    except RuntimeError:
        pass
    latitude_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(latitude_formatter)

    # optional: draw polygons & dots on the map
    draw_polygon(ax, kwargs)
    draw_dots(ax, ds, args, kwargs)

    title = get_title(args, ds)
    ax.set_title(title, y=1.01)

    filename = get_filename(args)
    plt.savefig(os.path.join(PLOTPATH, filename), dpi=300)
    logging.info('Saved. {}'.format(os.path.join(PLOTPATH, filename)))

    plt.close('all')


def read_input():
    """Parses input"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        dest='filename', metavar='FILENAME', type=str,
        help='Valid /path/to/file.nc')

    parser.add_argument(
        '--variable-name', '-varn', dest='varn', metavar='VARIABLE_NAME',
        type=str, default=None, help='Variable name included in file.nc')

    parser.add_argument(
        '--months', '-m', dest='months', default=None,
        type=lambda x: list(map(int, x.split(','))),
        help='Months as comma-separated integers')
    parser.add_argument(
        '--period', '-p', dest='period', default=None,
        # ['2006-09-01', '2016-08-31'],
        type=lambda x: list(map(lambda y: str(y.strip()), x.split(','))),
        help='Time period in the form "yyyy-mm-dd, yyyy-mm-dd".')
    parser.add_argument(
        '--title', '-t', dest='title', type=str, default=None,
        help='Title string of the plot.')
    parser.add_argument(
        '--plot-filename', '-pf', dest='plotname', default=None, type=str,
        help='Name of the output file with extension (sets filetype)')

    parser.add_argument(
        '--config-filename', '-c', dest='config', default=None,
        type=str, help='Filename of a valid config file')

    parser.add_argument(
        '--indicator', '-i', dest='indicator', default=None, type=str,
        help=' '.join([
            'A valid variable name. Values should be [0, 1] otherwise a',
            'warning is logged and behaviour might be unexpected. Can be',
            'used to indicate blocking, significance.']))
    parser.add_argument(
        '--indicator-size', '-is', dest='indicator_size', default=None,
        type=float, help=' '.join([
            'Size of the indicator. Special convention: if < 0 size will',
            'be set depending in the value of the indicator variable',
            'with the maximum value shown as abs(indicator_size)']))

    args = parser.parse_args()

    logmsg = 'Read parser input: \n\n'
    for ii, jj in sorted(vars(args).items()):
        logmsg += '  {}: {}\n'.format(ii, jj)
    logging.info(logmsg)
    return args


def get_variable_name(ds, args):
    """Tries to automatically set the variable to plot."""
    if args.varn is None:
        varns = set(ds.variables.keys()).difference(ds.coords)
        if len(varns) == 1:
            args.varn = varns.pop()
        else:
            errmsg = 'More than one variable in data set! Specify varn'
            raise IOError(errmsg)


def get_config_filename(args):
    """Tries to automatically set the config to use."""
    if args.config is None:
        try:
            if args.varn == 'Blocking':
                args.config = 'plot_map_config_blocking'
            if args.varn == 'GeopotentialHeight':
                args.config = 'plot_map_config_gph'
            # add config names here
        except IOError:
            pass


def main():
    """Runs the program"""
    logging.info('Running program ' + __file__)
    args = read_input()

    ds = xarray.open_dataset(args.filename)

    get_variable_name(ds, args)
    get_config_filename(args)

    time_name = ut.get_time_name(ds)
    if time_name is not None:
        ds = ut.get_time_subset(
            ds, time_name, period=args.period, months=args.months)
        period = ds[time_name][[0,-1]].data
        if period[0] == period[-1]: period = [period[0]]
        args.period = map(lambda x: str(x)[:10], period)
        ds = ds.mean(time_name, keep_attrs=True).squeeze()

    plot(ds, args)


if __name__ == '__main__':
    main()
