#!/usr/bin/env python
# -*- coding: utf-8 -*-
# a test
"""
Author:
- Lukas Brunner || lukas.brunner@uni-graz.at

Abstract:

"""
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
import os
import argparse
import xarray
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import seaborn.apionly as sns
import seaborn as sns
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import blocking.utils as ut
mpl.rc('font', **{'size': 11})

from matplotlib.ticker import MaxNLocator

import plot_map_config as config

PLOTPATH = './../plots'


def get_season(months, str_='{}'):
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


def get_title(args):
    if args.title is not None:
        return args.title

    return '{} {} to {}{}'.format(
        args.varn, args.period[0], args.period[1],
        get_season(args.months, ' ({})'))


def get_filename(args):
    if args.plotname is not None:
        return args.plotname

    filename = 'Mapplot_{}'.format(args.varn)
    if args.period is not None:
        filename += '_{}'.format('-'.join(args.period))
    filename += get_season(args.months, '_{}') + '.png'
    return filename


def plot(ds, args):

    kwargs = __import__(args.config)

    fig, ax = plt.subplots(
        figsize=(8, 4),
        subplot_kw={
            'projection': ccrs.PlateCarree(central_longitude=0)})
    fig.subplots_adjust(**kwargs.subplots_adjust)

    ax.set_global()
    ax.coastlines()
    cmap = mpl.pyplot.cm.get_cmap('YlOrRd')

    lon_name = ut.get_longitude_name(ds)
    lat_name = ut.get_latitude_name(ds)
    ds = ds.transpose(lon_name, lat_name)
    import ipdb; ipdb.set_trace()
    pc = ax.pcolormesh(ds[lon_name], ds[lat_name], ds[args.varn],
                       **kwargs.pcolormesh)

    # --- make colorbar same size as map ---
    def resize_colorbar(event):
        plt.draw()
        posn = ax.get_position()
        cbar_ax.set_position(
            [posn.x0 + posn.width + .01, posn.y0, .04, posn.height])

    cbar_ax = fig.add_axes([0, 0, .1, .1])
    fig.canvas.mpl_connect('resize_event', resize_colorbar)
    plt.colorbar(
        pc, cax=cbar_ax, **kwargs.colorbar)
    resize_colorbar(None)
    # --- done ---


    ax.set_xlabel('Longitude')
    ax.set_xticks(range(-180, 180+1, 60), crs=ccrs.PlateCarree())
    longitude_formatter = LongitudeFormatter()
    ax.xaxis.set_major_formatter(longitude_formatter)

    ax.set_ylabel('Latitude')
    ax.set_yticks(range(-90, 90+1, 30), crs=ccrs.PlateCarree())
    latitude_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(latitude_formatter)

    title = get_title(args)
    ax.set_title(title, y=1.01)

    filename = get_filename(args)
    plt.savefig(os.path.join(PLOTPATH, filename), dpi=300)
    logging.info('Saved. {}'.format(os.path.join(PLOTPATH, filename)))

    plt.close('all')


def read_input():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        dest='filename', metavar='FILENAME', type=str,
        help='Valid /path/to/file.nc')
    parser.add_argument(
        dest='varn', metavar='VARIABLE_NAME', type=str,
        help='Variable name included in file.nc')

    parser.add_argument(
        '--months', '-m', dest='months', default=None,
        type=lambda x: list(map(int, x.split(','))),
        help='Months as comma-separated integers')
    parser.add_argument(
        '--period', '-p', dest='period', default=None,
        #['2006-09-01', '2016-08-31'],
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

    args = parser.parse_args()

    if args.config is None:
        if args.varn == 'Blocking':
            args.config = 'plot_map_config_blocking'
        else:
            raise NotImplementedError('Write config first!')

    logmsg = 'Read parser input: \n\n'
    for ii, jj in sorted(vars(args).iteritems()):
        logmsg += '  {}: {}\n'.format(ii, jj)
    logging.info(logmsg)
    return args


def main():
    logging.info('Running program ' + __file__)
    args = read_input()

    ds = xarray.open_dataset(args.filename)

    time_name = ut.get_time_name(ds)
    if time_name is not None:
        ds = ut.get_time_subset(
            ds, time_name, period=args.period, months=args.months)
        ds = ds.mean(time_name).squeeze()

    plot(ds, args)


if __name__ == '__main__':
    main()
