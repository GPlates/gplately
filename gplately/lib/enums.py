#
#    Copyright (C) 2024-2026 The University of Sydney, Australia
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

from enum import Enum


# The numbers representing the registration types are the same with PyGMT's grid registration enum values
class GridRegistration(Enum):
    # good for geoscience data because sampled at locations(points) on the grid
    Gridline = 0
    # good for image data, the data represents the average value of an area
    Pixel = 1


class LongitudeConvention(Enum):
    POSITIVE_360 = "0_360"
    SIGNED_180 = "-180_180"
    POSITIVE_180 = "0_180"  # regional grids within longitude range 0-180, not able to decide whether it is POSITIVE_360 or SIGNED_180
    OTHER = "other"  # for other longitude range not belonging to the above conventions, e.g. (360,720), (-540,-180), etc.
