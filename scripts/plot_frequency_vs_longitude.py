#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Copyright (C) 2017 Lukas Brunner (Wegener Center/University of Graz)

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

Abstract: TODO

"""

import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
import argparse
import os
import xarray
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn.apionly as sns
from cartopy.mpl.ticker import LongitudeFormatter

import blocking.utils as ut

matplotlib.rc('font', **{'size': 12.5})

OUTPATH = './../plots/'


def roll_longitude(ds, lon_name, lon_start=-180):
    """
    """

    # map longitudes to the range [-180, 180[
    ds = ds.assign(
        **{lon_name: map(lambda x: x - 360 if x > 180 else x, ds[lon_name])})

    # roll longitude dim so that lon_start is at pos 0
    [idx], = np.where(ds[lon_name] == lon_start)
    ds = ds.roll(**{lon_name: ds.dims[lon_name]-idx})

    # make values increase for plot
    if lon_start >= 0:
        # map to [lon_start-360, lon_start[
        ds = ds.assign(
            **{lon_name: map(lambda x: x - 360 if x >= lon_start else x, ds[lon_name].data)})
    else:
        # map to [lon_start, lon_start+360[
        ds = ds.assign(
            **{lon_name: map(lambda x: x + 360 if x < lon_start else x, ds[lon_name].data)})
    return ds


def get_plot_kwargs(varn, memorize=[]):

    def get_color(varn):
        colors = sns.color_palette('colorblind', 6)
        if 'Lejenaes' in varn: return colors[0]
        elif 'Tibaldi' in varn: return colors[1]
        elif 'Scherrer' in varn: return colors[2]
        elif 'Davini' in varn: return colors[3]
        elif 'Brunner' in varn: return colors[4]
        else: return colors[5]

    def get_linestyle(varn):
        varn = varn.split('_')[0]
        memorize.append(varn)
        linestyles = ['-', '--', ':']
        return linestyles[(np.array(memorize)==varn).sum()-1]

    def get_label(varn):
        if varn == 'Lejenaes': return u'Lejenäs'
        elif varn == 'Tibaldi': return 'Tibaldi'
        elif varn == 'Scherrer': return 'Scherrer'
        elif varn == 'Scherrer_IB': return 'Scherrer IB'
        elif varn == 'Davini': return 'Davini'
        elif varn == 'Davini_eIB': return 'Davini eIB'
        elif varn == 'Brunner': return 'Brunner'
        elif varn == 'Brunnner_eIB': return 'Brunner eIB'
        elif varn == 'Brunner_IB': return 'Brunner IB'
        else: return varn

    return {'color': get_color(varn),
            'linestyle': get_linestyle(varn),
            'label': get_label(varn)}


def plot(ds_list, args):
    logging.info('Plot data')

    colors = sns.color_palette('colorblind', len(ds_list))
    fig = plt.figure(figsize=(9, 4.5))
    # fig.subplots_adjust(left=.07, bottom=0.12,
    #                     right=.99, top=.93)

    for idx, (varn, ds) in enumerate(ds_list):
        lon_name = ut.get_longitude_name(ds)
        ds = roll_longitude(ds, lon_name, args.lon_start)

        var, lons = ds[varn].data, ds[lon_name].data
        plt.plot(lons, var, color=colors[idx], label=args.labels[idx])
        # ds[varn].plot(**get_plot_kwargs(varn))

    xticks = plt.gca().get_xticks()
    if all(map(lambda x: float(x).is_integer(), xticks)):
        xticks = xticks.astype(int)

    xticklabels = []
    for tick in xticks:
        if tick >= 180:
            tick -= 360
        if tick == 0:
            xticklabels.append(u'0\xb0')
        elif tick == 180 or tick == -180:
            xticklabels.append(u'180\xb0')
        elif tick > 0:
            xticklabels.append(u'{}\xb0E'.format(tick))
        else:
            xticklabels.append(u'{}\xb0W'.format(-tick))
    plt.gca().set_xticklabels(xticklabels)

    plt.gca().set_xlim(lons[0], lons[-1])
    plt.gca().xaxis.grid()
    plt.gca().set_xlabel(u'Longitude')

    yticks = plt.gca().get_yticks()
    yticks = yticks[yticks>=0]
    yticklabels = np.round(100*yticks, 10)
    if all(map(lambda x: float(x).is_integer(), yticklabels)):
        yticklabels = yticklabels.astype(int)

    plt.gca().set_ylim(-.01, None)
    plt.gca().set_yticks(yticks)
    plt.gca().set_yticklabels(yticklabels)
    plt.gca().yaxis.grid()
    plt.gca().set_ylabel('Frequency (%)')

    if args.title is None:
        title = '1D blocking frequency'
        if args.season is not None:
            title += ' {}'.format(args.season)
    else:
        title = args.title

    plt.gca().set_title(title)
    plt.gca().legend(fontsize='small', ncol=2)

    filename = 'BlockingLongitude'
    if args.season is not None:
        filename += '_{}'.format(args.season)
    plt.savefig(os.path.join(OUTPATH, filename + '.png'), dpi=300)
    plt.close()
    logging.info('Saved. {}'.format(os.path.join(path, filename + '.png')))


def read_input():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--filenames', '-f', dest='filenames',
        nargs='*', type=str, required=True,
        help='At least one valid filename with the blocking information')
    parser.add_argument(
        '--variable-names', '-v', dest='varns',
        nargs='*', type=str, required=True,
        help='List of variable names with same length as filenames')
    parser.add_argument(
        '--labels', '-l', dest='labels',
        nargs='*', type=str, default=None,
        help='List of labels with same length as filenames')
    parser.add_argument(
        '--season', '-s', dest='season', type=str, default=None,
        choices=['JJA', 'SON', 'DJF', 'MAM'])
    parser.add_argument(
        '--longitude-start', '-lon-start', dest='lon_start',
        type=float, default=-90.,
        help='Start value for the x-axis [-180, 180[ (default=-90).')
    parser.add_argument(
        '--title', '-t', dest='title', type=str, default=None,
        help='Plot title')

    args = parser.parse_args()

    if args.labels is None:
        args.labels = [None] * len(args.filenames)

    logmsg = 'Read parser input: \n\n'
    for ii, jj in sorted(vars(args).iteritems()):
        logmsg += '  {}: {}\n'.format(ii, jj)
    logging.info(logmsg)

    return args


def main():
    args = read_input()
    ds_list = []
    for varn, filename in zip(args.varns, args.filenames):
        ds = xarray.open_dataset(filename)

        # ---
        time_name = ut.get_time_name(ds)
        if time_name is not None:
            if args.season is not None:
                ds = ut.get_season(ds, time_name, args.season)
                ds.load()
            ds = ds.mean(time_name)
        # ---

        ds_list.append((varn, ds))

    plot(ds_list, args)


if __name__ == '__main__':
    main()
