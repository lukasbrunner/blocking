#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
utils provides some utility functions

Copyright (c) 2017 Lukas Brunner (Wegener Center/University of Graz)

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
"""
from numpy.core import datetime64
import xarray

def calculate_daily_mean(ds, time_name='time'):
    """
    Calculates the daily mean.
    """
    return ds.resample(
        '1d', dim=time_name, how='mean', keep_attrs=True)


def calculate_gph_from_gp(ds,
                          gp_name='z',
                          gp_unit='m**2 s**-2',
                          gph_name='GeopotentialHeight'):
    """
    Creates a new variable with name gph_name from the variable gp_name
    by dividing it through the mean gravitational accelerating g=9.80665
    m s**-2.

    Parameters:
    - gp_name='z' (str, optional): Name of the variable
      containing the geopotential
    - gp_unit='m**2 s**-2' (str, optional): Unit of gp_name
    - gph_name='GeopotentialHeight' (str, optional): Name of the newly
      created variable containing the geopotential height
    """
    g = 9.80665  # m s**-2
    # https://dx.doi.org/10.6028/NIST.SP.330e2008
    if ds[gp_name].attrs['units'] != gp_unit:
        errmsg = 'Geopotential unit should be {} not {}'.format(
            gp_unit, ds[gp_name].attrs['units'])
        raise ValueError(errmsg)

    ds[gph_name] = xarray.Variable(
        ds.variables[gp_name].dims,
        ds.variables[gp_name].data / g,
        attrs={
            'units': 'm',
            'history': 'Calculated from {} with g={}'.format(gp_name, g)})
    return ds


def reduce_to_1D(ds, latitude_range, latitude_name='Latitude',
                 time_mean=True, time_name='Time'):
    """
    TODO
    """
    if isinstance(latitude_range, (int, float)):
        latitude_range = [latitude_range, latitude_range]
    elif len(latitude_range) == 1:
        latitude_range = [latitude_range[0], latitude_range[0]]
    elif len(latitude_range) != 2:
        errmsg = ' '.join(['latitude_range has to be float or list of',
                           'one or two floats and not {}']).format(
                               latitude_range)
        raise ValueError(errmsg)
    lats = ds[latitude_name].data
    lats = lats[(lats >= latitude_range[0]) & (lats <= latitude_range[1])]
    ds = ds.sel(**{latitude_name: lats})
    ds = (ds.sum(latitude_name) > 0).astype(int)

    if time_mean:
        ds = ds.mean(time_name)
    return ds


def get_longitude_name(ds):
    for varn in ds.variables.keys():
        if ('units' in ds[varn].attrs and
            ds[varn].attrs['units'] in ['degree_east', 'degrees_east']):
            return varn
    return None


def get_latitude_name(ds):
    for varn in ds.variables.keys():
        if ('units' in ds[varn].attrs and
            ds[varn].attrs['units'] in ['degree_north', 'degrees_north']):
            return varn
    return None


def get_time_name(ds):
    for varn in ds.variables.keys():
        if (('units' in ds[varn].attrs and
             'since' in ds[varn].attrs['units']) or
            ('units' in ds[varn].encoding and
             'since' in ds[varn].encoding['units'])):
            return varn
    for varn in ds.variables.keys():
        try:
            var = ds[varn].data[0]
        except IndexError:
            var = ds[varn].data
        if isinstance(var, datetime64):
            return varn
    return None


def get_season(ds, time_name, season):
    ds = ds.sel(**{time_name: ds['{}.season'.format(time_name)]==season})
    return ds


def get_months(ds, time_name, months):
    if isinstance(months, int):
        months = [months]
    ds = ds.sel(**{time_name: map(
        lambda x: x in months, ds['{}.month'.format(time_name)])})
    return ds


def get_period(ds, time_name, period):
    if len(period) == 1:
        ds = ds.sel(**{time_name: period[0]})
    elif len(period) == 2:
        ds = ds.sel(**{time_name: slice(*period)})
    return ds


def get_time_subset(ds, time_name, period=None, months=None, season=None):
    if period is not None:
        ds = get_period(ds, time_name, period)
    if months is not None:
        ds =  get_months(ds, time_name, months)
    if season is not None:
        ds = get_season(ds, time_name, season)
    if None in [period, months, season]:
        ds.load()
    return ds
