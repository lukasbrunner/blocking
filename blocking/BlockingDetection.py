#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BlockingDetection is a blocking detection algorithm based on xarray.

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

import logging
import numpy as np
import xarray
from scipy.ndimage import measurements

import blocking.utils as ut

class Blocking(object):
    """
    An xarray-based blocking detection algorithm. For a theoretical
    discussion see Brunner, L. (PhD thesis, in prep.).
    """

    def __repr__(self):
        return 'Class {}: \n{}'.format(self.__class__.__name__, self.ds)

    def __str__(self):
        return self.__repr__()

    def __init__(self, verbose=True, debug=False):
        """
        Creates an empty Blocking instance.

        Parameters:
        - verbose, debug (bool, optional): Set logging level (default is verbose)
        """
        if debug:
            logging.basicConfig(
                level=logging.DEBUG,
                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            logging.debug('Initialized  empty Blocking() instance in debug mode')
        elif verbose:
            logging.basicConfig(
                level=logging.INFO,
                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        else:
            logging.basicConfig(
                level=logging.WARNING,
                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        self.ds = None
        self._grid = {}

    def import_xarray(self, ds):
        """
        Imports an existing xarray data set.

        Parameters:
        - ds (data set): Valid xarray data set.
        """
        if self.ds is None:
            if not isinstance(ds, xarray.core.dataset.Dataset):
                errmsg = 'ds has to be a xarray data set!'
                raise ValueError(errmsg)
            self.ds = ds
            logging.debug('import_xarray: {}'.format(self.__str__))
        else:
            errmsg = 'Blocking() is already set!'
            raise ValueError(errmsg)

    def read(self, filename):
        """
        Reads a file into a xarray dataset.

        Parameters:
        - filename (str): Valid path + filename
        """
        if self.ds is None:
            self.ds = xarray.open_dataset(filename)
            logging.debug('read: {}'.format(self.__str__))
        else:
            errmsg = 'Blocking() is already set!'
            raise ValueError(errmsg)

    def save(self, path, save_vars=None, **kwargs):
        """
        Saves data set as netCDF.

        Parameters:
        - path (str): Valid path and filename to save.
        - save_vars=None ([str | list], optional): If not None only variables
          in save_vars will be saved.
        - **kwargs (dict, optional): xarray.dataset.to_netcdf kwargs
        """
        if save_vars is not None:
            if isinstance(save_vars, str):
                save_vars = [save_vars]
            ds = xarray.Dataset(attrs=self.ds.attrs)
            for varn in save_vars:
                ds[varn] = self.ds[varn]
                ds.to_netcdf(path, **kwargs)
        else:
            self.ds.to_netcdf(path, **kwargs)
        logging.info('Saved {}'.format(path))

    def set_up(self,
               time_name=None,
               longitude_name=None,
               latitude_name=None,
               pressure_name=None,
               pressure_level=50000.,
               pressure_unit='pa',
               variable_name='GeopotentialHeight',
               variable_unit='m',
               keep_vars=None):
        """
        Prepares the dataset for blocking detection. Does consistency checks
        and tests if all required information is available. Sets internal
        variables.

        Parameters:
        - *_name (str, optional): Name of required variables in the dataset.
          If None will try to identify them by unit. Exception: pressure has
          to be set if it is contained in dataset. TODO: update this
        - pressure_level=50000. (float, optional): Pressure level to select in
          units of pressure_unit. Has to exist in dataset. Can be set to None if
          pressure is not contained in the dataset.
        - pressure_unit='pa' (str, optional): Unit of pressure_level.
        - variable_unit='m', (str, optional): Unit of variable.
        """

        self._set_time_grid(time_name)
        self._set_longitude_grid(longitude_name)
        self._set_latitude_grid(latitude_name)
        self._set_pressure(pressure_name, pressure_unit, pressure_level)

        self._set_variable(variable_name, variable_unit)
        self._select_variables(keep_vars)

    def calculate_gradients(self,
                            delta_degree=15,
                            delta_index=None,
                            gradient_name='GeopotentialHeightGradient'):
        """
        Calculates discrete gradients over a given latitude delta.
        For more information refer to Brunner and Steiner (2017)

        Parameters:
        - delta_degree (float, optional): Delta in degree latitude (default is 15)
        - delta_index (int, optional): Delta as number of gird points
        - gradient_name (str, optional): Basename of the newly created variables

        Note:
        Either delta_degree or delta_index have to be give. If both are given
        they have to lead to the same delta_index or a ValueError will be raised
        """
        logging.info('Calculate gradients as {}...'.format(gradient_name))
        didx = self._convert_dimension_to_index(
            self._latitude_name, delta_degree, delta_index)
        dgrad = abs(self.ds[self._latitude_name][0] - self.ds[self._latitude_name][didx]).data

        var, lat_axis = self._swap_to_front(
            self.ds[self._variable_name], self._latitude_name)
        var = var.data
        nan_array = np.empty((didx,) + var.shape[1:]) * np.nan
        attrs = {'units': 'm degree**-1',
                 'latitude_delta': didx,
                 'latitude_delta_unit': 'grid_points'}

        gradient_north = np.concatenate(
            (var[didx:] - var[:-didx], nan_array)).swapaxes(0, lat_axis) / dgrad
        self.gphgn_name = '{}_north'.format(gradient_name)
        self.ds[self.gphgn_name] = xarray.Variable(
            self.ds[self._variable_name].dims,
            gradient_north, attrs=attrs)
        del gradient_north

        gradient_south = np.concatenate(
            (nan_array, var[:-didx] - var[didx:])).swapaxes(0, lat_axis) / dgrad
        self.gphgs_name = '{}_south'.format(gradient_name)
        self.ds[self.gphgs_name] = xarray.Variable(
            self.ds[self._variable_name].dims,
            gradient_south, attrs=attrs)
        del gradient_south

        idx_eq = np.abs(self.ds.variables[self._latitude_name].data).argmin()
        gradient_equator = np.concatenate((
            var[2*didx:idx_eq+didx] - var[didx:idx_eq],
            nan_array, nan_array,
            var[idx_eq-didx:-2*didx] - var[idx_eq:-didx])).swapaxes(0, lat_axis) / dgrad

        self.gphge_name = '{}_equator'.format(gradient_name)
        self.ds[self.gphge_name] = xarray.Variable(
            self.ds[self._variable_name].dims,
            gradient_equator, attrs=attrs)
        del gradient_equator, var
        # setattr(self, '{}_equator'.format(gradient_name),
        #         self.ds['{}_equator'.format(gradient_name)])
        logging.info('Calculating gradients... DONE')

    def calculate_ib(self,
                     gradient_equator_below=0,
                     gradient_pole_below=-10,
                     gradient_equator2_above=5,
                     ib_name='IB',
                     pressure_level=50000):
        """
        Calculates instantaneous blocking (IB) from gradients.

        Parameters:
        - gradient_* (float, optional): Gradient criteria.
          gradient_equator2_above and/or gradient_pole_below can be set to None
          to be omitted.
        - ib_name (str, optional): Name of the newly created variable
        - pressure_level (float, optional): If given and if a pressure level
          was set in set_up() has to match this level or a ValueError will be
          raised. Will be omitted if not pressure exists in the dataset.

        """
        logging.info('Calculating IB as {}...'.format(ib_name))
        if (self.gphgn_name is None or
            self.gphgs_name is None or
            self.gphge_name is None or
            self.gphgn_name not in self.ds or
            self.gphgs_name not in self.ds or
            self.gphge_name not in self.ds):
            errmsg = ' '.join([
                'Gradients not defined or not found.',
                'Run calculate_gradients first?'])
            raise ValueError(errmsg)
        if self._pressure_name is not None and \
           self._pressure_name in self.ds.dims:

            if self._pressure_level is not None and \
               self._pressure_level != pressure_level:
                errmsg = ' '.join([
                    'Pressure level was level was set to {} {} in set_up()',
                    'setting a different pressure level ({}) here is not',
                    'allowed!']).format(
                        self._pressure_level,
                        self._pressure_unit,
                        pressure_level)
                raise ValueError(errmsg)

            ds = self.ds.sel(**{self._pressure_name:pressure_level})
        else:
            ds = self.ds

        lats = ds.variables[self._latitude_name].data
        ds_nh = ds.sel(**{self._latitude_name:lats[lats>=0]})
        ds_sh = ds.sel(**{self._latitude_name:lats[lats<0]})

        def _test_thresholds(g1, g2, g3, t1, t2, t3):
            if t3 is None and t1 is None:
                return np.where(g2 < t2, 1, 0)
            elif t3 is None:
                return np.where((g1 < t1) & (g2 < t2), 1, 0)
            elif t1 is None:
                return np.where((g1 < t1) & (g3 > t3), 1, 0)
            return np.where((g1 < t1) & (g2 < t2) & (g3 > t3), 1, 0)

        old = np.seterr(all='ignore')
        ib_nh = _test_thresholds(
            ds_nh.variables[self.gphgn_name],
            ds_nh.variables[self.gphgs_name],
            ds_nh.variables[self.gphge_name],
            gradient_pole_below,
            gradient_equator_below,
            gradient_equator2_above)

        ib_sh = _test_thresholds(
            ds_sh.variables[self.gphgs_name],
            ds_sh.variables[self.gphgn_name],
            ds_sh.variables[self.gphge_name],
            gradient_pole_below,
            gradient_equator_below,
            gradient_equator2_above)
        np.seterr(**old)

        self.ib_name = ib_name
        if gradient_equator2_above is None:
            gradient_equator2_above = 'not set'
        if gradient_pole_below is None:
            gradient_pole_below = 'not set'
        dims = self.ds[self.gphgn_name].dims
        self.ds[self.ib_name] = xarray.Variable(
            dims, np.concatenate(
                (ib_sh, ib_nh),
                axis=dims.index(self._latitude_name)),
            attrs={
                'gradient_pole_below': gradient_pole_below,
                'gradient_equator_below': gradient_equator_below,
                'gradient_equator2_above': gradient_equator2_above})
        del ib_nh, ib_sh, ds_nh, ds_sh, ds
        logging.info('Calculating IB... DONE')

    def calculate_eib(self,
                      min_extent_degree=15,
                      min_extent_index=None,
                      eib_name='ExtendedIB'):
        """
        Calculates extended IB from IB.

        Parameters:
        - min_extent_degree (float, optional): Minimum extent in degree
          (default is 15)
        - min_extent_index (int, optional): Minimum extent as number of
          grid points
        - eib_name (str, optioinal): Name of the newly created variable

        Note: Either min_extent_degree or min_extend_index have to be given.
        If both are given they have to be consistent.
        """
        logging.info('Calculating extended IB as {}...'.format(eib_name))
        didx = self._convert_dimension_to_index(
            self._longitude_name, min_extent_degree, min_extent_index)

        var, lon_axis = self._swap_to_front(
            self.ds[self.ib_name], self._longitude_name)
        var = var.data # TODO: do empty like before and eib will become a xarray automatically
        eib = np.zeros_like(var)
        var_add = var[:didx]
        var = np.concatenate((var, var_add))  # allow structures to cross date border

        structure = np.zeros((3,) * len(var.shape))
        # TODO: does this have to be so complicated?
        idx = (slice(None, None, None),) + (1,)*(len(var.shape) - 1)
        structure[idx] = 1  # find structures along lon_axis (=first axis)

        labels, _ = measurements.label(var, structure=structure)
        slices = measurements.find_objects(labels)
        for slice_ in slices:
            if (slice_[0].stop - slice_[0].start) >= didx:
                eib[slice_] = 1  # set Extended IB

                # NOTE (Info for line above & date border crossing):
                # If slice_[0].stop is larger then len of axis it will
                # just take the axis till the end. Shift it to cover the
                # beginning here:
                if slice_[0].stop > self.ds.dims[self._longitude_name]:
                    slice_shifted = (
                        (slice(0, slice_[0].stop -
                               self.ds.dims[self._longitude_name], 1),) + slice_[1:])
                    eib[slice_shifted] = 1

        eib = eib.swapaxes(0, lon_axis)
        self.eib_name = eib_name
        self.ds[self.eib_name] = xarray.Variable(
            self.ds[self.ib_name].dims,
            eib, attrs={
                'min_longitude_extent':didx,
                'min_longitude_extent_unit': 'grid_points'})
        del var, eib, labels, slices
        logging.info('Calculating extended IB... DONE')

    def calculate_blocking(self,
                           stationary_pm_days=2,
                           stationary_pm_index=None,
                           longitude_pm_degree=7.5,
                           longitude_pm_index=None,
                           latitude_pm_degree=2.5,
                           latitude_pm_index=None,
                           blocking_name='Blocking'):
        """
        Calculates Blocking from ExtendedIB.

        Parameters:
        - stationary_pm_days (float, optional): Number of days for
          persistency criterion: central days plus/minus
          (default=2 means 5 days in total)
        - stationary_pm_index (int, optional): Number of time steps
          for persistency criterion
        - longitude_pm_degree, latitude_pm_degree (float, optional):
          Plus/minus degree longitude and latitude around the center that
          the block is allowed to move during the given time range
          (default=7.5/2.5 lon/lat)
        - longitude_pm_index, latitude_pm_index (int optional):
          Like above but for grid points instead of degree
        - blocking_name (str, optional): Name of the newly created variable

        Note: One of each index/value pair has to be given. If both are given
        they have to be consistent.
        """
        logging.info('Calculating blocking as {}...'.format(blocking_name))
        didx_time = self._convert_dimension_to_index(
            self._time_name, stationary_pm_days, stationary_pm_index)
        didx_lon = self._convert_dimension_to_index(
            self._longitude_name, longitude_pm_degree, longitude_pm_index)
        didx_lat = self._convert_dimension_to_index(
            self._latitude_name, latitude_pm_degree, latitude_pm_index)

        # Here finally I need to know exactly where which dimension is
        dims = self.ds[self.eib_name].dims
        sort = [dims.index(dim) for dim in [self._time_name,
                                            self._longitude_name,
                                            self._latitude_name]]
        var = self.ds.transpose(dims[sort[0]], dims[sort[1]], dims[sort[2]])
        var = var[self.eib_name]
        blk = np.zeros_like(var)
        for i_time in range(self.ds.dims[self._time_name]):
            if 10.*i_time % self.ds.dims[self._time_name] < 10:
                logging.info('Calculating time step... {:.0%}'.format(
                    float(i_time) / self.ds.dims[self._time_name]))

            eib = var.isel(**{self._time_name:i_time})

            # Weak border criterion: if neighbor is out of range remove it
            idx = [i_time - dd for dd in range(1, didx_time+1) if i_time >= dd]
            idx += [i_time + dd for dd in range(1, didx_time+1)
                    if i_time + dd < self.ds.dims[self._time_name]]
            eib_neighbours = var.isel(**{self._time_name:idx})

            # TODO: implement strong border criterion: if neighbor is out of range
            # shift it to the other side

            for i_lon, i_lat in zip(*np.where(eib.data==1)):
                lons = map(lambda x: x - self.ds.dims[self._longitude_name]
                           if x >= self.ds.dims[self._longitude_name] else x,
                           range(i_lon - didx_lon, i_lon + didx_lon + 1))

                # TODO: latitude could theoretically be out of range with a
                # very large value of latitude_pm -> check here
                lats = range(i_lat - didx_lat, i_lat + didx_lat + 1)

                # ---
                # TODO: could allow pm 4 days (instead of pm 2) and demand
                # eib at at least 4 of them (that could then be as before pm 2
                # or  -1, 1, 2, 3 or 1, 2, 3, 4 etc)
                # TODO: should probably still keep it consecutive
                # ---

                box = eib_neighbours.isel(**{
                    self._longitude_name:lons,
                    self._latitude_name:lats})
                if np.all(box.sum(dim=[self._longitude_name,
                                       self._latitude_name])) > 0:
                    # Weak blocking criterion
                    blk[i_time, i_lon, i_lat] = 1

                    # TODO: stronger blocking criterion: also set all the eibs on other
                    # days to blocked
        self.blocking_name = blocking_name
        self.ds[self.blocking_name] = xarray.Variable(
            dims, blk.transpose(sort),
            attrs={
                'stationary_pm_index': didx_time,
                'longitude_pm_index': didx_lon,
                'latitude_pm_index': didx_lat})
        del blk
        logging.info('Calculating blocking... DONE')

    def reduce_to_1D(self, latitude_range, time_mean=True, inplace=True):
        ds = ut.reduce_to_1D(self.ds,
                             latitude_range=latitude_range,
                             latitude_name=self._latitude_name,
                             time_mean=time_mean,
                             time_name=self._time_name)
        if not inplace:
            return ds
        self.ds = ds

    def calculate_gph_from_gp(self,
                              gp_name='z',
                              gp_unit='m**2 s**-2',
                              gph_name='GeopotentialHeight'):
        self.ds = ut.calculate_gph_from_gp(self.ds, gp_name, gp_unit, gph_name)
        logging.info('Calculated GPH from GP')

    def get_time_subset(
            self, time_name=None, period=None, months=None, season=None):
        if time_name is None:
            time_name = ut.get_time_name(self.ds)
        self.ds = ut.get_time_subset(
            self.ds, time_name=time_name, period=period,
            months=months, season=season)
        logmsg = 'Selected time subset:'
        if period is not None:
            logmsg += ' {}'.format(' to '.join(period))
        if months is not None:
            logmsg += ' months: {}'.format(', '.join(months))
        if season is not None:
            logmsg += ' {}'.format(season)
        logging.info(logmsg)

    def calculate_daily_mean(self, time_name=None):
        if time_name is None:
            time_name = ut.get_time_name(self.ds)
        self.ds = ut.calculate_daily_mean(self.ds, time_name)
        logging.info('Calculated daily mean')

    def _convert_dimension_to_index(self, dimn, dim=None, idx=None):
        if dim is None and idx is None:
            errmsg = ' '.join([
                'Either a value or a index (or both) have to be given',
                'for dimension {}']).format(dimn)
            raise ValueError(errmsg)
        elif dim is not None:
            delta_index = dim / float(self._grid[dimn])
            if not delta_index.is_integer():
                errmsg = '{}={} is not dividable by grid delta={}'.format(
                    dimn, dim, self._grid[dimn])
                raise ValueError(errmsg)
            if idx is not None and idx != delta_index:
                errmsg = ' '.join([
                    'Value and index are given for dimension {} but they are',
                    'not unambiguous: {} / {} != {}']).format(
                        dimn, dim, idx, self._grid[dimn])
            return int(delta_index)
        return int(idx)

    def _set_time_grid(self, time_name):
        if time_name is None:
            self._time_name = ut.get_time_name(self.ds)
        else:
            self._time_name = time_name
        if self._time_name is None:
            errmsg = 'Name of time dimension was not given and could be found!'
            raise ValueError(errmsg)
        self._set_grid(self._time_name)


    # TODO: delete
    # def _get_variable(self, possible_units, varn=None):
    #     """
    #     Return a variable name based on given dimension choices. Can be used
    #     to detect standard dimension like time, longitude, latitude.
    #     NOTE: Since this class works on a xarray.dataset with decode_cf=True
    #     the time dimension does not have an unit attribute -> time unit is
    #     stored in .encoding

    #     Parameters:
    #     - possible units ([str | list]): Unit string identifying a dimension.
    #       Possible options: 'since' (time), degree(s)_east (longitude),
    #       degree(s)_north (latitude)
    #     - varn=None (str, optional): If given varn has to be a valid dimension
    #       and possible_units will be ignored.
    #     """
    #     if isinstance(possible_units, str):
    #         possible_units = [possible_units]
    #     if varn is not None and varn not in self.ds:  # if given has to be in ds
    #         errmsg = '{} not found in dataset'.format(varn)
    #         raise ValueError(errmsg)
    #     else:  # try to find by unit
    #         for vv in self.ds.variables.keys():
    #             if (('units' in self.ds[vv].attrs and
    #                 self.ds[vv].attrs['units'] in possible_units) or
    #                 ('units' in self.ds[vv].encoding and
    #                  self.ds[vv].encoding['units'].split(' ')[1] in possible_units)):

    #                 if varn is None:
    #                     varn = vv
    #                     logmsg = 'Detected {}'.format(varn)
    #                     logging.info(logmsg)
    #                 else:
    #                     errmsg = ' '.join([
    #                         '{0} already set! Two {0}',
    #                         'dimensions?']).format(varn)
    #                     raise ValueError(errmsg)
    #         if varn is None:
    #             errmsg = 'No dimension detected with units "{}"'.format(
    #                 ', '.join(possible_units))
    #             raise ValueError(errmsg)
    #     return varn

    def _set_grid(self, varn, allow_inverse=False):
        if varn == self._time_name:
            var = self.ds[varn].to_index().to_julian_date()
        else:
            var = self.ds[varn].data
        delta = np.unique(var[1:] - var[:-1])
        if len(delta) != 1:
            errmsg = 'No regular grid found for dimension {}'.format(varn)
            raise ValueError(errmsg)
        if delta[0] == 0:
            errmsg = 'Two equivalent values found for dimension {}.'.format(
                varn)
            raise ValueError(errmsg)
        elif delta[0] < 0:
            if allow_inverse:
                logmsg = '{} grid not increasing. Inverse dimension!'.format(varn)
                logging.warning(logmsg)
                self.ds = self.ds.reindex({varn: np.sort(self.ds[varn].data)})
                delta = -delta
            else:
                errmsg = ' '.join(['{} grid not increasing. This should',
                                   'not happen?!']).format(varn)
                raise ValueError(errmsg)

        self._grid[varn] = float(delta[0])
        logmsg = 'Set {} grid distance {}'.format(varn, delta[0])
        logging.info(logmsg)

    def _set_longitude_grid(self, lon_name):
        if lon_name is None:
            self._longitude_name = ut.get_longitude_name(self.ds)
        else:
            self._longitude_name = lon_name
        if self._longitude_name is None:
            errmsg = ' '.join(['Name of longitude dimension was not',
                               'given and could not be found'])
            raise ValueError(errmsg)
        self._set_grid(self._longitude_name)

    def _set_latitude_grid(self, lat_name):
        if lat_name is None:
            self._latitude_name = ut.get_latitude_name(self.ds)
        else:
            self._latitude_name = lat_name
        if self._latitude_name is None:
            errmsg = ' '.join(['Name of latitude dimension was not',
                               'given and could not be found'])
            raise ValueError(errmsg)
        self._set_grid(self._latitude_name, allow_inverse=True)

    def _set_pressure(self, varn, unit, var):
        if varn is None:
            logmsg = 'No pressure level selected'
            logging.debug(logmsg)
        elif varn not in self.ds:
            errmsg = '{} not found in the dataset'.format(varn)
            raise ValueError(errmsg)
        else:
            logmsg = 'Detected {}'.format(varn)
            logging.info(logmsg)

            # check pressure level and unit only if pressure was found
            if unit != self.ds[varn].attrs['units']:
                errmsg = 'Ambiguous unit for {}: {}!={}'.format(
                    varn, self.ds[varn].attrs['units'], unit)
                raise ValueError(errmsg)
            if var is None:
                logmsg = 'pressure_level=None, continue with all levels'
                logging.info(logmsg)
                raise NotImplementedError('pressure_level=None not yet implemented')
            else:
                if var in self.ds[varn].data:
                    logmsg = 'Value found in {}: {} {}'.format(
                        varn, var, unit)
                    logging.info(logmsg)
                else:
                    errmsg = 'Select {}: {} {}'.format(
                        varn, var, unit)
                    raise ValueError(errmsg)
        self._pressure_name = varn
        self._select_pressure(var)

    def _select_pressure(self, var):
        if self._pressure_name is not None and var is not None:
            self.ds = self.ds.sel(**{self._pressure_name:var})

    def _set_variable(self, varn, unit):
        if varn not in self.ds:
            errmsg = '{} not found in dataset'.format(varn)
            raise ValueError(errmsg)
        if self.ds[varn].attrs['units'] != unit:
            errmsg = 'Ambiguous unit for {}: {}!={}'.format(
                varn, self.ds[varn].attrs['units'], unit)
            raise ValueError(errmsg)
        if len(set(self.ds[varn].dims).difference([
                self._time_name,
                self._longitude_name,
                self._latitude_name])) != 0:
            errmsg = '{} has to have dimensions {}, {}, {}, not {}'.format(
                varn, self._time_name, self._longitude_name,
                self._latitude_name, self.ds[varn].dims)
            raise ValueError(errmsg)
        self._variable_name = varn

    def _swap_to_front(self, var, dimn):
        dims = list(var.dims)
        axis = dims.index(dimn)
        dims[0], dims[axis] = dims[axis], dims[0]
        return var.transpose(*dims), axis

    def _select_variables(self, keep_vars):
        if keep_vars is None:
            keep_vars = []
        else:
            if isinstance(keep_vars, str):
                keep_vars = [keep_vars]

        all_vars = self.ds.variables.keys()
        keep_vars += [self._time_name,
                      self._longitude_name,
                      self._latitude_name,
                      self._variable_name]
        drop_vars = set(all_vars).difference(keep_vars)
        self.ds = self.ds.drop(drop_vars)
