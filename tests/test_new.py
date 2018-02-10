#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author:
- Lukas Brunner || lukas.brunner@live.at

Abstract:

"""


import os
import xarray as xr
import copy as cp
import numpy as np
import unittest
import numpy.ma.testutils as np_test
import warnings

import matplotlib.pyplot as plt

from blocking.BlockingDetection import Blocking



# create a GPH for the southern hemisphere (mirrored to NH later)
var = np.tile(np.linspace(5000, 6000, 36), 144*5).reshape(5, 144, -1).swapaxes(1, 2)
# var(time, latitude, longitude)

# set artificial 'blocking' cases:

# special case: -180/180 degree border
var[:, 11:12, 0] = var[:, 11+6:12+6, 0]
var[:, 11:12, -1] =  var[:, 11+6:12+6, -1]

var[:, 13:15, 0] = var[:, 13+6:15+6, 0]
var[:, 13:15, -1] = var[:, 13+6:15+6, -1]

# 'ladder' of 1x1 boxes
for lat, lon in zip(range(30), range(2, 35+2, 1)):
    var[:, lat:lat+1, lon:lon+1] = var[:, lat+6:lat+7, lon:lon+1]

# 'fence' of 1x20 stripes
for lat, lon in zip(range(10), range(38, 38+24, 2)):
    var[:, lat:lat+21, lon:lon+1] = np.tile(var[:, lat+26, lon:lon+1], 21).reshape(5, 21, -1)

# 'horizontal fence of ?x1 stripes with ? increasing
lon = 63
for idx, lat in enumerate(range(6, 30, 3)):
    lon += idx
    var[:, lat:lat+1, lon:lon+1+idx] = var[:, lat+6:lat+7, lon:lon+1+idx]

# 'blocks' too short for extended IB (interrupted by one box)
var[:, 11:15, 100:103] = var[:, 11+6:15+6, 100:103]
var[:, 11:15, 104:107] = var[:, 11+6:15+6, 104:107]

var[:, 17:21, 108:111] = var[:, 17+6:21+6, 108:111]
var[:, 17:21, 112:115] = var[:, 17+6:21+6, 112:115]

# create xarray
data = np.concatenate((var, var[:, ::-1]), axis=1)
ds = xr.Dataset({
    'time': ('time', range(5), {'units': 'days since 2000-01-01'}),
    'longitude': ('longitude', np.arange(0, 360, 2.5), {'units': 'degrees_east'}),
    'latitude': ('latitude', np.arange(-88.75, 90, 2.5), {'units': 'degree_north'}),
    'GeopotentialHeight': (('time', 'latitude', 'longitude'), data, {'units': 'm'})})

ds = xr.decode_cf(ds)

blk = Blocking()
blk.import_xarray(ds)

blk.set_up()

blk.calculate_gradients(
    delta_index=6)
# equator gradient for 'blocking' is exactly zero

blk.calculate_ib(
    gradient_equator_below=0,
    gradient_pole_below=None,
    gradient_equator2_above=None,
    ib_name='no_ib')

blk.calculate_ib(
    gradient_equator_below=.1,
    gradient_pole_below=-10,
    gradient_equator2_above=5,
    ib_name='ib_all')

blk.calculate_ib(
    gradient_equator_below=.1,  # gradient is exactly 0
    gradient_pole_below=None,
    gradient_equator2_above=None,
    ib_name='ib_eq_only')

blk.calculate_ib(
    gradient_equator_below=.1,
    gradient_pole_below=-10,
    gradient_equator2_above=None,
    ib_name='ib_eq_pole')

blk.save('testfile1.nc')
