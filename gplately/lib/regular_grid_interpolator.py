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

# pyright: reportMissingImports=false
# pyright: reportMissingModuleSource=false

from scipy.interpolate import RegularGridInterpolator as _RGI
import numpy as np


#
# Note: Michael Chin moved this class into this new file on 2026-06-25 and he didn't change code of this class.
#
class RegularGridInterpolator(_RGI):
    """A class to sample gridded data at a set of point coordinates using either linear or nearest-neighbour
    interpolation methods. It is a child class of `scipy 1.10`'s [`RegularGridInterpolator`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html) class.

    This will only work for scipy version 1.10 onwards.

    Attributes
    ----------
    points : tuple of ndarrays of float with shapes (m1, ), …, (mn, )
        Each array contains point coordinates that define the regular grid in n dimensions.
    values : ndarray
        The data on a regular grid. Note: the number of rows corresponds to the number of point latitudes, while the number
        of columns corresponds to the number of point longitudes.
    method : str, default=`linear`
        The method of interpolation to perform. Supported are "linear" and "nearest". Assumes “linear” by default.
    bounds_error : bool, default=false
        Choose whether to return a ValueError and terminate the interpolation if any provided sample points are out
        of grid bounds. By default, it is set to `False`. In this case, all out-of-bound point values are replaced
        with the `fill_value` (defined below) if supplied.
    fill_value : float, default=np.nan
        Used to replace point values that are out of grid bounds, provided that `bounds_error` is false.

    """

    def __init__(
        self, points, values, method="linear", bounds_error=False, fill_value=np.nan
    ):
        super(RegularGridInterpolator, self).__init__(
            points, values, method, bounds_error, fill_value
        )

    def __call__(self, xi, method=None, return_indices=False, return_distances=False):
        """Samples gridded data at a set of point coordinates. Uses either a linear or nearest-neighbour interpolation `method`.

        Uses the gridded data specified in the sample_grid method parameter. Note: if any provided sample points are out of
        grid bounds and a corresponding error message was suppressed (by specifying bounds_error=False), all out-of-bound
        point values are replaced with the self.fill_value attribute ascribed to the RegularGridInterpolator object (if it
        exists). Terminates otherwise.

        This is identical to scipy 1.10's RGI object.

        Parameters
        ----------
        xi : ndarray of shape (..., ndim)
            The coordinates of points to sample the gridded data at.

        method : str, default=None
            The method of interpolation to perform. Supported are "linear" and "Nearest". Assumes “linear” interpolation
            if None provided.

        return_indices : bool, default=False
            Choose whether to return indices of neighbouring sampling points.

        return_distances : bool, default=False
            Choose whether to return normal distances between interpolated points and neighbouring sampling points.

        Returns
        -------
        output_tuple : tuple of ndarrays
            The first ndarray in the output tuple holds the interpolated grid data. If sample point distances and indices are
            required, these are returned as subsequent tuple elements.

        Raises
        ------
        ValueError
            * Raised if the string method supplied is not “linear” or “nearest”.
            * Raised if the provided sample points for interpolation (xi) do not have the same dimensions as the supplied grid.
            * Raised if the provided sample points for interpolation include any point out of grid bounds. Alerts user which
            dimension (index) the point is located. Only raised if the RegularGridInterpolator attribute bounds_error is set
            to True. If suppressed, out-of-bound points are replaced with a set fill_value.
        """
        method = self.method if method is None else method
        if method not in ["linear", "nearest"]:
            raise ValueError("Method '%s' is not defined" % method)

        xi, xi_shape, ndim, nans, out_of_bounds = self._prepare_xi(xi)

        indices, norm_distances = self._find_indices(xi.T)

        if method == "linear":
            result = self._evaluate_linear(indices, norm_distances)
        elif method == "nearest":
            result = self._evaluate_nearest(indices, norm_distances)
        if not self.bounds_error and self.fill_value is not None:
            result[out_of_bounds] = self.fill_value  # type: ignore

        interp_output = result.reshape(xi_shape[:-1] + self.values.shape[ndim:])  # type: ignore
        output_tuple = [interp_output]

        if return_indices:
            output_tuple.append(indices)
        if return_distances:
            output_tuple.append(norm_distances)

        if return_distances or return_indices:
            return tuple(output_tuple)
        else:
            return output_tuple[0]

    def _prepare_xi(self, xi):
        try:
            from scipy.interpolate.interpnd import _ndim_coords_from_arrays
        except ImportError:
            # SciPy 1.15 renamed interpnd to _interpnd (see https://github.com/scipy/scipy/pull/21754).
            from scipy.interpolate._interpnd import _ndim_coords_from_arrays

        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError(
                "The requested sample points xi have dimension "
                f"{xi.shape[-1]} but this "
                f"RegularGridInterpolator has dimension {ndim}"
            )

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        # find nans in input
        nans = np.any(np.isnan(xi), axis=-1)

        if self.bounds_error:
            for i, p in enumerate(xi.T):
                if not np.logical_and(
                    np.all(self.grid[i][0] <= p), np.all(p <= self.grid[i][-1])
                ):
                    raise ValueError(
                        "One of the requested xi is out of bounds "
                        "in dimension %d" % i
                    )
            out_of_bounds = None
        else:
            out_of_bounds = self._find_out_of_bounds(xi.T)

        return xi, xi_shape, ndim, nans, out_of_bounds

    def _find_out_of_bounds(self, xi):
        # check for out of bounds xi
        out_of_bounds = np.zeros((xi.shape[1]), dtype=bool)
        # iterate through dimensions
        for x, grid in zip(xi, self.grid):
            out_of_bounds += x < grid[0]
            out_of_bounds += x > grid[-1]
        return out_of_bounds

    def _find_indices(self, xi):
        """Index identifier outsourced from scipy 1.9's
        RegularGridInterpolator to ensure stable
        operations with all versions of scipy >1.0.
        """
        # find relevant edges between which xi are situated
        indices = []
        # compute distance to lower edge in unity units
        norm_distances = []
        # iterate through dimensions
        for x, grid in zip(xi, self.grid):
            i = np.searchsorted(grid, x) - 1
            i[i < 0] = 0
            i[i > grid.size - 2] = grid.size - 2
            indices.append(i)

            # compute norm_distances, incl length-1 grids,
            # where `grid[i+1] == grid[i]`
            denom = grid[i + 1] - grid[i]
            with np.errstate(divide="ignore", invalid="ignore"):
                norm_dist = np.where(denom != 0, (x - grid[i]) / denom, 0)
            norm_distances.append(norm_dist)

        return indices, norm_distances

    def _evaluate_linear(self, indices, norm_distances):
        """Linear interpolator outsourced from scipy 1.9's
        RegularGridInterpolator to ensure stable
        operations with all versions of scipy >1.0.
        """
        import itertools

        # slice for broadcasting over trailing dimensions in self.values
        vslice = (slice(None),) + (None,) * (self.values.ndim - len(indices))

        # Compute shifting up front before zipping everything together
        shift_norm_distances = [1 - yi for yi in norm_distances]
        shift_indices = [i + 1 for i in indices]

        # The formula for linear interpolation in 2d takes the form:
        # values = self.values[(i0, i1)] * (1 - y0) * (1 - y1) + \
        #          self.values[(i0, i1 + 1)] * (1 - y0) * y1 + \
        #          self.values[(i0 + 1, i1)] * y0 * (1 - y1) + \
        #          self.values[(i0 + 1, i1 + 1)] * y0 * y1
        # We pair i with 1 - yi (zipped1) and i + 1 with yi (zipped2)
        zipped1 = zip(indices, shift_norm_distances)
        zipped2 = zip(shift_indices, norm_distances)

        # Take all products of zipped1 and zipped2 and iterate over them
        # to get the terms in the above formula. This corresponds to iterating
        # over the vertices of a hypercube.
        hypercube = itertools.product(*zip(zipped1, zipped2))
        values = 0.0
        for h in hypercube:
            edge_indices, weights = zip(*h)
            weight = 1.0
            for w in weights:
                weight *= w
            values += np.asarray(self.values[edge_indices]) * weight[vslice]  # type: ignore
        return values

    def _evaluate_nearest(self, indices, norm_distances):
        """Nearest neighbour interpolator outsourced from scipy 1.9's
        RegularGridInterpolator to ensure stable
        operations with all versions of scipy >1.0.
        """
        idx_res = [
            np.where(yi <= 0.5, i, i + 1) for i, yi in zip(indices, norm_distances)
        ]
        return self.values[tuple(idx_res)]
