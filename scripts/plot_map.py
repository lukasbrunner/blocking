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
import seaborn.apionly as sns
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import blocking.utils as ut
mpl.rc('font', **{'size': 11})


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


def plot_title(args):
    if args.title is not None:
        return args.title
    str_ = '{} {} to {}{}'.format(
        '{}', args.period[0], args.period[1],
        get_season(args.months, ' ({})'))
    if args.varn == 'GeopotentialHeight':
        return str_.format('Mean $500\,$hPa geopotential height [m]')
    elif args.varn == 'GeopotentialHeightGradient_north':
        return str_.format('Mean $500\,$hPa geopotential height gradient north [m/deg.]')
    elif args.varn == 'GeopotentialHeightGradient_south':
        return str_.format('Mean $500\,$hPa geopotential height gradient south [m/deg.]')
    elif args.varn == 'GeopotentialHeightGradient_equator':
        return str_.format('Mean $500\,$hPa geopotential height gradient equator [m/deg.]')
    elif args.varn == 'Blocking':
        return str_.format('Blocking frequency')
    else:
        return str_.format(args.varn)


def plot(ds, args):
    if args.varn == 'Blocking':
        levels = np.arange(0., .03+.0001, .001)
        ticks = np.arange(0, .03+.001, .005)
    elif args.varn == 'IB':
        levels = np.arange(0., .08+.001, .002)
        ticks = np.arange(0, .08+.001, .01)
    elif args.varn == 'GeopotentialHeight':
        levels = np.arange(4800, 5900+1, 50)
        ticks = levels[::2]
    def plot_properties(varn):
        kwargs = {'extend': 'both'}
        if 'Gradient' in varn:
            kwargs['levels'] = np.arange(-20, 20+1, 2)
            # kwargs['colors'] = sns.color_palette('RdBu_r', 21)
        elif varn == 'GeopotentialHeight':
            kwargs['levels'] = levels
        elif varn in ['Blocking', 'IB']:
            # kwargs['levels'] = levels
            cmap = plt.cm.get_cmap('YlOrRd')
            colors = cmap(np.linspace(0, len(levels)/(len(levels)+1.), cmap.N))
            # https://stackoverflow.com/questions/40982050/matplotlib-how-to-cut-the-unwanted-part-of-a-colorbar#40983666
            cmap_new = mpl.colors.LinearSegmentedColormap.from_list('temp', colors)
            cmap_new.set_under('w')
            cmap_new.set_over(cmap(np.linspace(0, 1, len(levels)+10)[-1]))
            kwargs['extend'] = 'max'
            kwargs['cmap'] = cmap_new
            kwargs['vmin'] = levels[1]
            # http://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html
            kwargs['norm'] = mpl.colors.BoundaryNorm(
                levels, ncolors=cmap_new.N, clip=False)
        return kwargs

    fig, ax = plt.subplots(
        figsize=(8, 4),
        subplot_kw={
            'projection': ccrs.PlateCarree(central_longitude=0)})
    fig.subplots_adjust(hspace=0, wspace=0,
                        left=0.1, bottom=0.05,
                        right=1.06, top=.99)

    ax.set_global()
    ax.coastlines()

    p = ds[args.varn].plot.pcolormesh(
        ax=ax,
        transform=ccrs.PlateCarree(),
        **plot_properties(args.varn))

    # Plot contour at 5800 m (nice Omega shape)
    color = sns.color_palette('colorblind', 6)[2:]
    ds_sel = ds.sel(**{'Longitude': np.arange(0, 87.5+.1, 2.5),
                       'Latitude': np.arange(30, 80+.1, 2.5)})
    ax.contour(ds_sel[args.varn].Longitude,
               ds_sel[args.varn].Latitude,
               ds_sel[args.varn], [5800],
               linewidths=2.,
               colors=color
    )

    ax.set_xlabel('Longitude')
    ax.set_xticks(range(-180, 180+1, 60), crs=ccrs.PlateCarree())
    longitude_formatter = LongitudeFormatter()
    ax.xaxis.set_major_formatter(longitude_formatter)

    ax.set_ylabel('Latitude')
    ax.set_yticks(range(-90, 90+1, 30), crs=ccrs.PlateCarree())
    latitude_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(latitude_formatter)


    poly = mpl.patches.Polygon([[-10, 30], [-10, 75], [100, 75], [100, 30]],
                               closed=True, facecolor='none', edgecolor='black', lw=2, ls='-')
    ax.add_patch(poly)

    # For plot (a)
    # poly = mpl.patches.Polygon([[-180, 55], [180, 55], [180, 65], [-180, 65]],
    #                            closed=True, facecolor='none', edgecolor='black', lw=1)
    # ax.add_patch(poly)
    # poly = mpl.patches.Polygon([[-180, -55], [180, -55], [180, -65], [-180, -65]],
    #                            closed=True, facecolor='none', edgecolor='black', lw=1)
    # ax.add_patch(poly)

    # For plot (b)
    # poly = mpl.patches.Polygon([[-180, 50], [180, 50]], closed=False,
    #                            facecolor='none', edgecolor='black', lw=1)
    # ax.add_patch(poly)

    # poly = mpl.patches.Polygon([[-180, -50], [180, -50]], closed=False,
    #                            facecolor='none', edgecolor='black', lw=1)
    # ax.add_patch(poly)

    def resize_colorbar(event):
        plt.draw()
        posn = ax.get_position()
        cbar_ax.set_position(
            [posn.x0 + posn.width + .01, posn.y0, .04, posn.height])
    p.colorbar.remove()  # plot my own colorbar
    cbar_ax = fig.add_axes([0, 0, .1, .1])
    fig.canvas.mpl_connect('resize_event', resize_colorbar)
    plt.colorbar(
        p, cax=cbar_ax,
        ticks=ticks)
    resize_colorbar(None)

    title = plot_title(args)
    ax.set_title(title, y=1.01)

    filename = 'Mapplot_{}'.format(args.varn)
    if args.period is not None:
        filename += '_{}'.format('-'.join(args.period))
    filename += get_season(args.months, '_{}') + '.png'
    plt.savefig(os.path.join('./', filename), format='png', dpi=300)
    logging.info('Saved. {}'.format(os.path.join('./', filename)))
    plt.close('all')


def read_input():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--filename', '-f', dest='filename', type=str, required=True)
    parser.add_argument(
        '--variable-name', '-v', dest='varn', type=str, required=True)
    parser.add_argument(
        '--months', '-m', dest='months', default=None,
        type=lambda x: list(map(int, x.split(','))),
        help='Months as comma-separated integers (default=None)')
    parser.add_argument(
        '--period', '-p', dest='period', default=['2006-09-01', '2016-08-31'],
        type=lambda x: list(map(lambda y: str(y.strip()), x.split(','))))
    parser.add_argument(
        '--title', '-t', dest='title', default=None)

    args = parser.parse_args()
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
    lon_name = ut.get_longitude_name(ds)
    lat_name = ut.get_latitude_name(ds)

    if time_name is not None:
        ds = ut.get_time_subset(ds, time_name, period=args.period, months=args.months)
        ds = ds.mean(time_name).squeeze()
    ds = ds.transpose(lat_name, lon_name)

    plot(ds, args)


if __name__ == '__main__':
    main()
