#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculates a blocking index following Tibaldi and Molteni (1990)

Copyright (C) 2017 Lukas Brunner

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from BlockingDetection import Blocking

INPATH = './../data/erai_geopotential_500hpa_20160101_20161231.nc'
OUTPATH = INPATH.replace(
    'geopotential', 'blocking').replace('.nc', '_Tibaldi1D.nc')

blk = Blocking()

blk.read(INPATH)

blk.set_up()

blk.calculate_gradients(delta_degree=20)

blk.calculate_ib(
    gradient_equator_below=0,
    gradient_pole_below=-10,
    gradient_equator2_above=None,
    ib_name='Blocking')

blk.reduce_to_1D([55., 65.])

blk.save(OUTPATH, 'Blocking')
