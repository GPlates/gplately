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

import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def _as_2d_float_array(data, name):
    arr = numpy.ma.getdata(numpy.asanyarray(data)).astype(float, copy=True)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be a 2-D array, got shape {arr.shape}.")
    return arr


def _normalize_finite(arr, fill_value=0.0):
    normalized = numpy.full(arr.shape, fill_value, dtype=float)
    finite = numpy.isfinite(arr)
    if not numpy.any(finite):
        return normalized

    arr_min = numpy.nanmin(arr)
    arr_max = numpy.nanmax(arr)
    if arr_max <= arr_min:
        return normalized

    normalized[finite] = (arr[finite] - arr_min) / (arr_max - arr_min)
    return normalized


def set_shade(
    a,
    intensity=None,
    cmap=None,
    vmin=None,
    vmax=None,
    scale=10.0,
    azdeg=165.0,
    altdeg=45.0,
):
    """Apply soft-light shading to a 2-D grid.

    Parameters
    ----------
    a : array-like
        A 2-D data array.
    intensity : array-like, optional
        A 2-D intensity source used to derive hillshade. If omitted, ``a`` is
        used.
    cmap : matplotlib colormap, optional
        Colormap used to map ``a`` into RGB before blending.
    vmin, vmax : float, optional
        Data bounds for colormap normalization. If omitted, finite min/max of
        ``a`` are used.
    scale, azdeg, altdeg : float, optional
        Hillshade parameters passed to :func:`hillshade`.

    Returns
    -------
    numpy.ndarray
        RGB image array with shape ``(ny, nx, 3)``.
    """
    if cmap is None:
        cmap = getattr(plt, "get_cmap")("jet")

    data = _as_2d_float_array(a, "a")
    finite_data = numpy.isfinite(data)

    if intensity is not None:
        intensity = _as_2d_float_array(intensity, "intensity")
        if intensity.shape != data.shape:
            raise ValueError(
                "intensity must have the same shape as a: "
                f"{intensity.shape} != {data.shape}."
            )

    if intensity is None:
        intensity_img = hillshade(data, scale=scale, azdeg=azdeg, altdeg=altdeg)
    else:
        intensity_img = hillshade(intensity, scale=scale, azdeg=azdeg, altdeg=altdeg)

    if not numpy.any(finite_data):
        data_norm = numpy.zeros_like(data)
    else:
        if vmin is None:
            vmin = numpy.nanmin(data)
        if vmax is None:
            vmax = numpy.nanmax(data)

        if vmax <= vmin:
            data_norm = numpy.zeros_like(data)
        else:
            data_norm = numpy.zeros_like(data)
            data_norm[finite_data] = (data[finite_data] - vmin) / float(vmax - vmin)

    rgb = cmap(data_norm)[:, :, :3]

    # Use neutral illumination on invalid cells so their base color is preserved.
    intensity_img[~finite_data] = 0.5

    # Form an RGB equivalent of intensity.
    d = intensity_img.repeat(3).reshape(rgb.shape)
    # simulate illumination based on pegtop algorithm.
    rgb = 2 * d * rgb + (rgb**2) * (1 - 2 * d)
    return rgb


def hillshade(data, scale=10.0, azdeg=165.0, altdeg=45.0):
    """Convert a 2-D array to normalized hillshade intensity.

    Returns values in the ``[0, 1]`` range, with ``0.5`` as neutral intensity
    when variation is insufficient to estimate relief.
    """
    data = _as_2d_float_array(data, "data")
    if scale <= 0:
        raise ValueError("scale must be > 0.")

    finite = numpy.isfinite(data)
    if not numpy.any(finite):
        return numpy.full(data.shape, 0.5, dtype=float)

    # Fill invalid values with the finite mean to reduce edge artifacts in gradients.
    fill_value = numpy.nanmean(data[finite])
    filled = numpy.where(finite, data, fill_value)

    # convert alt, az to radians
    az = azdeg * numpy.pi / 180.0
    alt = altdeg * numpy.pi / 180.0
    # gradient in x and y directions
    dx, dy = numpy.gradient(filled / float(scale))
    slope = 0.5 * numpy.pi - numpy.arctan(numpy.hypot(dx, dy))
    aspect = numpy.arctan2(dx, dy)
    intensity = numpy.sin(alt) * numpy.sin(slope) + numpy.cos(alt) * numpy.cos(
        slope
    ) * numpy.cos(-az - aspect - 0.5 * numpy.pi)

    intensity = _normalize_finite(intensity, fill_value=0.5)
    intensity[~finite] = 0.5
    return intensity


def get_topo_cmap():
    """Return a topography/bathymetry colormap used by gplately plots."""
    colors = [
        (-10927, [10, 0, 121]),
        (-10500, [26, 0, 137]),
        (-10000, [38, 0, 152]),
        (-9500, [27, 3, 166]),
        (-9000, [16, 6, 180]),
        (-8500, [5, 9, 193]),
        (-8000, [0, 14, 203]),
        (-7500, [0, 22, 210]),
        (-7000, [0, 30, 216]),
        (-6500, [0, 39, 223]),
        (-6000, [12, 68, 231]),
        (-5500, [26, 102, 240]),
        (-5000, [19, 117, 244]),
        (-4500, [14, 133, 249]),
        (-4000, [21, 158, 252]),
        (-3500, [30, 178, 255]),
        (-3000, [43, 186, 255]),
        (-2500, [55, 193, 255]),
        (-2000, [65, 200, 255]),
        (-1500, [79, 210, 255]),
        (-1000, [94, 223, 255]),
        (-500, [138, 227, 255]),
        (-0.001, [188, 230, 255]),
        (-0.0005, [51, 102, 0]),
        (100, [51, 204, 102]),
        (200, [187, 228, 146]),
        (500, [255, 220, 185]),
        (1000, [243, 202, 137]),
        (1500, [230, 184, 88]),
        (2000, [217, 166, 39]),
        (2500, [168, 154, 31]),
        (3000, [164, 144, 25]),
        (3500, [162, 134, 19]),
        (4000, [159, 123, 13]),
        (4500, [156, 113, 7]),
        (5000, [153, 102, 0]),
        (5500, [162, 89, 89]),
        (6000, [178, 118, 118]),
        (6500, [183, 147, 147]),
        (7000, [194, 176, 176]),
        (7500, [204, 204, 204]),
        (8726, [229, 229, 229]),
    ]

    color_list = [
        (
            float(color[0] - colors[0][0]) / (colors[-1][0] - colors[0][0]),
            [x / 255.0 for x in color[1]],
        )
        for color in colors
    ]

    return LinearSegmentedColormap.from_list("topo", color_list, N=1024)
