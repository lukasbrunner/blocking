#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author:
- Lukas Brunner || lukas.brunner@uni-graz.at

Abstract:

"""
from numpy.core import datetime64


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
    # NOTE: if there is more than one time dimension this will only find
    # the first one
    for varn in ds.variables.keys():
        if ('units' in ds[varn].encoding and
            ds[varn].encoding['units'].split(' ')[1] == 'since'):
            return varn
    for varn in ds.variables.keys():
        try:
            var = ds[varn].data[0]
        except IndexError:
            var = ds[varn].data
        if isinstance(var, datetime64):
            return varn
    return None


def get_season(ds, season):
    time_name = get_time_name(ds)
    ds = ds.sel(**{time_name: ds['{}.season'.format(time_name)]==season})
    return ds

def get_months(ds, months):
    time_name = get_time_name(ds)
    if isinstance(months, int):
        months = [months]
    ds = ds.sel(**{time_name: map(
        lambda x: x in months, ds['{}.month'.format(time_name)])})
    return ds

def get_period(ds, period):
    time_name = get_time_name(ds)
    if len(period) == 1:
        ds = ds.sel(**{time_name: period[0]})
    elif len(period) == 2:
        ds = ds.sel(**{time_name: slice(*period)})
    return ds


def get_time_subset(ds, period=None, months=None, season=None):
    if period is not None:
        ds = get_period(ds, period)
    if months is not None:
        ds = get_season(ds, season)
    if season is not None:
        ds =  get_months(ds, months)
    if None in [period, months, season]:
        ds.load()
    return ds
