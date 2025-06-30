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
import math
from typing import Union

import pygplates
from plate_model_manager import PlateModelManager


def reconstruct_points(
    lons: list[float],
    lats: list[float],
    model_name: str,
    times: Union[float, list[float]],
    pids: Union[int, list[int], None] = None,
    valid_time: Union[tuple[float, float], list[tuple[float, float]], None] = None,
    anchor_plate_id: int = 0,
    ignore_valid_time: bool = False,
    reverse: bool = False,
):
    """Reconstruct points to the given time(s).

    Parameters
    ----------
    lons: a list of float
        The longitudes of points.
    lats: a list of float
        The latitudes of points.
    model_name: str
        The plate model name.
    times: float or a list of float
        The reconstruction time(s)
    pids: a list of int, int or None, optional, default=None
        Allow user to provide plate IDs to speed up the reconstruction.
        If user provided ``pids`` without ``valid_time``.
        The reconstruction will not consider the valid time range, same as ``ignore_valid_time``.
    valid_time: tuple of 2, a list of tuple of 2, or None, optional, default=None
        Allow user to provide valid times for points.
    anchor_plate_id: int, optional, default=0
        The anchor plate ID.
    ignore_valid_time: bool, optional default=False
        If True, reconstruct the points without considering the valid time range.
    reverse: bool, optional, default=False
        If True, do reverse reconstruction to present-day (0Ma).
        The reverse reconstruction works only with a single reconstruction time.

    Returns
    -------
    A list of :class:`dict`.
        The returned :class:`dict` contains:

        * lons
        * lats
        * pids (plate IDs)
        * begin_times
        * end_times
        * time (the reconstruction time)

    """
    _model = PlateModelManager().get_model(model_name)
    # user must provide a valid plate model name. See https://gplates.github.io/gplately/sphinx-latest/html/use_cases.html#id1
    assert _model

    static_polygon_fc = pygplates.FeatureCollection()
    static_polygon_files = _model.get_static_polygons()
    assert static_polygon_files
    for f in static_polygon_files:
        static_polygon_fc.add(pygplates.FeatureCollection(f))

    return reconstruct_points_impl(
        lons,
        lats,
        pygplates.RotationModel(_model.get_rotation_model()),
        pygplates.FeatureCollection(static_polygon_fc),
        times,
        pids,
        valid_time,
        anchor_plate_id,
        ignore_valid_time,
        reverse,
    )


def reconstruct_points_impl(
    lons: list[float],
    lats: list[float],
    rotation_model: pygplates.RotationModel,
    static_polygons: pygplates.FeatureCollection,
    times: Union[float, list[float]],
    pids: Union[int, list[int], None] = None,
    valid_time: Union[tuple[float, float], list[tuple[float, float]], None] = None,
    anchor_plate_id: int = 0,
    ignore_valid_time: bool = False,
    reverse: bool = False,
):
    """Similar to :func:`gplately.reconstruct_points`. The only difference is that this function allows
    user to provide rotation model and static polygons feature collection instead of a model name.

    Parameters
    ----------
    rotation_model: pygplates.RotationModel
        The rotation model as a ``pygplates.RotationModel`` object.
    static_polygons: pygplates.FeatureCollection
        The static polygons as a ``pygplates.FeatureCollection`` object.


    .. seealso::

        See :func:`gplately.reconstruct_points` for other parameter details.

    """
    # the length of lons and lats must be the same.
    assert len(lons) == len(lats)

    if isinstance(pids, list):
        # if user has provided a list of pids for the points, the length must be the same as lons and lats.
        assert len(lons) == len(pids)
    if pids is None or isinstance(pids, int):
        _pids = [pids] * len(lons)
    else:
        _pids = pids

    if isinstance(valid_time, list):
        # if user has provided a list of valid_time for the points, the length must be the same as lons and lats.
        assert len(lons) == len(valid_time)
    if valid_time is None or isinstance(valid_time, tuple):
        _valid_time = [valid_time] * len(lons)
    else:
        _valid_time = valid_time

    if isinstance(times, list):
        _times: list[float] = times  # type: ignore
    else:
        _times: list[float] = [float(times)]

    # create point features
    p_index = 0
    point_features = []
    for lat, lon, pid, v_time in zip(lats, lons, _pids, _valid_time):
        point_feature = pygplates.Feature()
        point_feature.set_geometry(pygplates.PointOnSphere(lat, lon))
        point_feature.set_name(str(p_index))  # type: ignore
        if pid is not None:
            point_feature.set_reconstruction_plate_id(pid)  # type: ignore
        if v_time is not None:
            # the valid time must be a tuple of 2 (begin_time, end_time)
            assert len(v_time) == 2
            point_feature.set_valid_time(v_time[0], v_time[1])  # type: ignore
        point_features.append(point_feature)
        p_index += 1

    # if user has provided plate id(s), do not partition(slow)
    if pids is None:
        properties_to_copy = [pygplates.PartitionProperty.reconstruction_plate_id]
        if not ignore_valid_time:
            properties_to_copy.append(pygplates.PartitionProperty.valid_time_period)

        if reverse:
            # if it is reverse reconstruction, the "times" must be a single time
            assert len(_times) == 1
            partition_time = _times[0]
        else:
            partition_time = 0.0

        # LOOK HERE !!!!
        # it seems when the partition_time is not 0
        # the returned assigned_point_features contains reverse reconstructed present-day geometries.
        # so, when doing reverse reconstruct, do not reverse reconstruct again later.
        assigned_point_features = pygplates.partition_into_plates(
            static_polygons,
            rotation_model,
            point_features,
            properties_to_copy=properties_to_copy,
            reconstruction_time=partition_time,  # type: ignore
        )
    else:
        assigned_point_features = point_features

    if reverse:
        # we still need reverse reconstruct if the points had not been partitioned.
        # if user had provided the pids, we did not do partition. So, we need to call pygplates.reverse_reconstruct in that case.
        if pids is not None:
            need_reverse_reconstruct = True
        else:
            need_reverse_reconstruct = False
        ret = _get_reverse_reconstructed_data(
            assigned_point_features,
            rotation_model,
            _times[0],  # can only do single time for reverse reconstruction
            anchor_plate_id,
            need_reverse_reconstruct,
        )

        return [ret]
    else:
        # handle normal reconstruction
        ret = []
        for t in _times:
            ret.append(
                _get_reconstructed_data(
                    assigned_point_features, rotation_model, t, anchor_plate_id
                )
            )
        return ret


def reverse_reconstruct_points(
    lons: list[float],
    lats: list[float],
    model_name: str,
    time: float,
    pids: Union[int, list[int], None] = None,
    valid_time: Union[tuple[float, float], list[tuple[float, float]], None] = None,
    anchor_plate_id: int = 0,
    ignore_valid_time: bool = False,
):
    """Wrapper function to reverse reconstruct points with :func:`gplately.reconstruct_points`.

    .. seealso::

        See :func:`gplately.reconstruct_points` for parameter details.
    """
    return reconstruct_points(
        lons,
        lats,
        model_name,
        times=time,
        pids=pids,
        valid_time=valid_time,
        anchor_plate_id=anchor_plate_id,
        ignore_valid_time=ignore_valid_time,
        reverse=True,
    )


def reverse_reconstruct_points_impl(
    lons: list[float],
    lats: list[float],
    rotation_model: pygplates.RotationModel,
    static_polygons: pygplates.FeatureCollection,
    time: float,
    pids: Union[int, list[int], None] = None,
    valid_time: Union[tuple[float, float], list[tuple[float, float]], None] = None,
    anchor_plate_id: int = 0,
    ignore_valid_time: bool = False,
):
    """Wrapper function to reverse reconstruct points with :func:`gplately.reconstruct_points_impl`.

    .. seealso::

        See :func:`gplately.reconstruct_points` for parameter details.
    """
    return reconstruct_points_impl(
        lons=lons,
        lats=lats,
        rotation_model=rotation_model,
        static_polygons=static_polygons,
        times=time,
        pids=pids,
        valid_time=valid_time,
        anchor_plate_id=anchor_plate_id,
        ignore_valid_time=ignore_valid_time,
        reverse=True,
    )


def _get_reconstructed_data(
    assigned_point_features, rotation_model, time, anchor_plate_id
):
    reconstructed_feature_geometries = []
    pygplates.reconstruct(  # type: ignore
        assigned_point_features,
        rotation_model,
        reconstructed_feature_geometries,
        time,
        anchor_plate_id=anchor_plate_id,
    )
    point_count = len(assigned_point_features)
    r_lons = point_count * [math.nan]
    r_lats = point_count * [math.nan]
    r_pids = point_count * [math.nan]
    r_btimes = point_count * [math.nan]
    r_etimes = point_count * [math.nan]
    for rfg in reconstructed_feature_geometries:
        lat, lon = rfg.get_reconstructed_geometry().to_lat_lon()
        feature = rfg.get_feature()
        idx = int(feature.get_name())
        r_lons[idx] = lon
        r_lats[idx] = lat
    for idx in range(len(assigned_point_features)):
        feature = assigned_point_features[idx]
        r_pids[idx] = feature.get_reconstruction_plate_id()
        btime, etime = feature.get_valid_time()
        r_btimes[idx] = btime
        r_etimes[idx] = etime

    return {
        "lons": r_lons,
        "lats": r_lats,
        "pids": r_pids,
        "begin_times": r_btimes,
        "end_times": r_etimes,
        "time": time,
    }


def _get_reverse_reconstructed_data(
    assigned_point_features,
    rotation_model,
    time,
    anchor_plate_id,
    need_reverse_reconstruct,
):
    if need_reverse_reconstruct:
        pygplates.reverse_reconstruct(  # type: ignore
            assigned_point_features,
            rotation_model,
            time,
            anchor_plate_id=anchor_plate_id,
        )
    point_count = len(assigned_point_features)
    r_lons = point_count * [math.nan]
    r_lats = point_count * [math.nan]
    r_pids = point_count * [math.nan]
    r_btimes = point_count * [math.nan]
    r_etimes = point_count * [math.nan]
    for feature in assigned_point_features:
        lat, lon = feature.get_geometry().to_lat_lon()
        idx = int(feature.get_name())
        r_lons[idx] = lon
        r_lats[idx] = lat
        r_pids[idx] = feature.get_reconstruction_plate_id()
        btime, etime = feature.get_valid_time()
        r_btimes[idx] = btime
        r_etimes[idx] = etime
    return {
        "lons": r_lons,
        "lats": r_lats,
        "pids": r_pids,
        "begin_times": r_btimes,
        "end_times": r_etimes,
        "time": time,
    }
