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
from matplotlib.colors import LinearSegmentedColormap, hsv_to_rgb


def parse_old_cpt_row(row, color_model="RGB"):
    """parse one row of text from a cpt file in old format, such as the four rows below from a cpt file

    0	210	0	0	10	210	0	0
    10	230	40	0	20	230	40	0
    20	245	60	0	30	245	60	0
    30	255	98	0	40	255	98	0

    Parameters
    ----------
    row : a list of strings
        A list of strings representing a line from a cpt file in old format.
    color_model : str, optional
        The color model to use, either "RGB" or "HSV". Default is "RGB".

    Returns
    -------
    tuple of list of float and list of list of float
        A tuple containing the values and colors(RGB values in a closed interval [0, 1]) parsed from the row.
    """

    values = []
    colors = []
    values.append(float(row[0]))
    values.append(float(row[4]))
    if color_model == "RGB":
        colors.append(
            [
                float(row[1]) / 255.0,
                float(row[2]) / 255.0,
                float(row[3]) / 255.0,
            ]
        )
        colors.append(
            [
                float(row[5]) / 255.0,
                float(row[6]) / 255.0,
                float(row[7]) / 255.0,
            ]
        )
    else:
        colors.append(
            list(hsv_to_rgb([float(row[1]) / 360.0, float(row[2]), float(row[3])]))
        )
        colors.append(
            list(hsv_to_rgb([float(row[5]) / 360.0, float(row[6]), float(row[7])]))
        )
    return values, colors


def parse_new_cpt_row(row, color_model="RGB"):
    """parse one row of text from a cpt file in new format, such as the four rows below from a cpt file

    0               64/0/64         0.052632        64/0/192        L
    0.052632        64/0/192        0.105263        0/64/255        L
    0.105263        0/64/255        0.157895        0/128/255       L
    0.157895        0/128/255       0.210526        0/160/255       L

    Parameters
    ----------
    row : a list of strings
        A list of strings representing a line from a cpt file in new format.
    color_model : str, optional
        The color model to use, either "RGB" or "HSV". Default is "RGB

    Returns
    -------
    tuple of list of float and list of list of float
        A tuple containing the values and colors(RGB values in a closed interval [0, 1]) parsed from the row.
    """
    # convert the new format into old format
    color_1 = row[1].split("/")
    color_2 = row[3].split("/")
    return parse_old_cpt_row(
        [row[0]] + color_1 + [row[2]] + color_2, color_model=color_model
    )


def get_cmap_from_gmt_cpt(cpt_file):
    """given a gmt cpt file path, return a `matplotlib.colors.LinearSegmentedColormap` object.

    Parameters
    ----------
    cpt_file : str
        The path to the gmt cpt file.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        A matplotlib colormap object.
    """
    values = []
    colors = []
    color_model = "RGB"
    with open(cpt_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if "COLOR_MODEL" in line and "HSV" in line:
                color_model = "HSV"
            if line[0] in ["#", "B", "F", "N"]:
                continue
            vals = line.split()
            if len(vals) == 8:
                vs, cs = parse_old_cpt_row(vals, color_model=color_model)
                values.extend(vs)
                colors.extend(cs)
            elif len(vals) == 5 or len(vals) == 4:
                vs, cs = parse_new_cpt_row(vals, color_model=color_model)
                values.extend(vs)
                colors.extend(cs)

    if values is None or colors is None:
        raise ValueError(f"No valid colour data found in cpt file: {cpt_file}")

    colour_list = []
    for i in range(len(values)):
        colour_list.append(
            ((values[i] - values[0]) / (values[-1] - values[0]), [x for x in colors[i]])
        )

    return LinearSegmentedColormap.from_list("cpt-cmap", colour_list)
