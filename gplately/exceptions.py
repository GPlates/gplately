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


class ReconstructionTimeNotSet(Exception):
    """raise this exception when the reconstruction time is None."""

    def __init__(self):
        super().__init__(
            "The reconstruction time has not been set yet. Set `PlotTopologies.time` before calling plotting functions."
        )


class UnableToGetModelList(Exception):
    """raise this exception when failed to get a list of model names over Internet"""

    def __init__(self):
        super().__init__(
            "Unable to get a list of model names over Internet. Your network may be down or the servers may be offline."
        )
