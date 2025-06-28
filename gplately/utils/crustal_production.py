#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
import logging
import os
import sys

import pygmt
import xarray as xr

logger = logging.getLogger("gplately")


def compute_crustal_production_rate(agegrid_fn: str, age_threshold: int = 3):
    """mask out all values greater than the age threshold in the age grid;
    find the total area which is covered by the remaining age values;
    divide the area by the age threshold to obtain crustal production rate in km2/year

    Parameters
    ----------
    agegrid_fn: str
        the file path to the age grid file(NetCDF)
    age_threshold: int
        use the age threshold to mask age grid and calculate the crustal production rate

    Returns
    -------
    a floating point number(crustal production rate in km2/year)
    """

    if not os.path.isfile(agegrid_fn):
        logging.error(
            f"The given age grid file({agegrid_fn}) does not exist or is not a file."
        )

    if age_threshold < 3:
        logging.warn(
            f"The give age threshold({age_threshold}) is too small. It is highly recommended that the value >=3."
        )

    age_grid = xr.open_dataset(agegrid_fn)
    data = age_grid.to_dataarray()

    # just in case the age grid netcdf is 3-dimensional
    if len(data.dims) == 3:
        data = data[0]

    data = data.where(data < age_threshold + sys.float_info.epsilon)

    output_dataframe = pygmt.grdvolume(
        grid=data,
        contour=[0],
        output_type="pandas",
        region=[-180, 180, -90, 90],
        unit="k",  # specify the area unit, 'k' means kilometers
        f="g",  # To make sure your Cartesian grid is recognized as geographical, use the -fg option.
    )
    # print(output_dataframe)
    million_years = 1000000
    # the rate is in km2/year
    rate = output_dataframe.iloc[0, 1] / age_threshold / million_years
    return round(rate, 2)
