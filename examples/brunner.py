#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculates a blocking index following Brunner et al. (2017)

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

"""

from blocking.BlockingDetection import Blocking

blk = Blocking()
blk.read('./../data/erai_geopotential_500hpa_20160101_20161231.nc')

blk.get_time_subset(period=['2016-03-01', '2016-06-01'], season='MAM')
blk.calculate_daily_mean()
blk.calculate_gph_from_gp()

blk.set_up()

blk.calculate_gradients(delta_degree=15)

blk.calculate_ib(
    gradient_equator_below=0,
    gradient_pole_below=-10,
    gradient_equator2_above=5)

blk.calculate_eib(min_extent_degree=15)

blk.calculate_blocking(
    stationary_pm_days=2,
    longitude_pm_degree=7.5,
    latitude_pm_degree=2.5)
blk.reduce_to_1D([50, 75])

blk.save('./../data/erai_blocking_20160101-20161231_Brunner1D.nc', 'Blocking')
