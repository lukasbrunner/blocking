blocking v0.1
=============

A Python xarray-based blocking detection algorithm
--------------------------------------------------

by Lukas Brunner (Wegener Center/University of Graz,
lukas.brunner@live.at)

To acknowledge this code please cite Brunner et al. ([2017](https://doi.org/10.5194/amt-10-4727-2017)) and Brunner ([2018](http://iacweb.ethz.ch/staff/lukbrunn/welcome/files/Brunner2018_PhD.pdf)).

Please also note the [LICENSE](./LICENSE)


Download sample data set (for examples and tests)
-------------------------------------------------

Go to the ERA-Interim web page (an account my be required):
http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=pl/

Download the variable **Geopotential** at the **500hPa** pressure level as netCDF file. For a nice blocking pattern you should download at least one year (about 30MB). Save the data set in the _data_ folder and set the _INPATH_ variable in the scripts accordingly.
