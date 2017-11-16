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
