#!/usr/bin/env python3
import os

from matplotlib.colors import LinearSegmentedColormap, hsv_to_rgb


def parse_old_cpt(row, values, colors, color_model="RGB"):
    """parse one row of the cpt file

    0	210	0	0	10	210	0	0
    10	230	40	0	20	230	40	0
    20	245	60	0	30	245	60	0
    30	255	98	0	40	255	98	0

    """
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


def parse_new_cpt(row, values, colors, color_model="RGB"):
    """parse one row of the new cpt file

    0               64/0/64         0.052632        64/0/192        L
    0.052632        64/0/192        0.105263        0/64/255        L
    0.105263        0/64/255        0.157895        0/128/255       L
    0.157895        0/128/255       0.210526        0/160/255       L
    0.210526        0/160/255       0.263158        64/192/255      L
    """
    # convert the new format into old format
    color_1 = row[1].split("/")
    color_2 = row[3].split("/")
    parse_old_cpt(
        [row[0]] + color_1 + [row[2]] + color_2, values, colors, color_model=color_model
    )


def get_cm_from_gmt_cpt(cpt_file):
    """given a gmt cpt file path, return a matplotlib matplotlib.colors.LinearSegmentedColormap object"""
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
                parse_old_cpt(vals, values, colors, color_model=color_model)
            elif len(vals) == 5 or len(vals) == 4:
                parse_new_cpt(vals, values, colors, color_model=color_model)

    colour_list = []
    # print(colors)
    for i in range(len(values)):
        colour_list.append(
            ((values[i] - values[0]) / (values[-1] - values[0]), [x for x in colors[i]])
        )
    return LinearSegmentedColormap.from_list("cpt-cmap", colour_list)


if __name__ == "__main__":
    # gmt makecpt -Chaxby > haxby.cpt
    file_path = os.path.dirname(os.path.realpath(__file__))
    print(get_cm_from_gmt_cpt(f"{file_path}/../data/GMT_haxby_old.cpt"))
    print(get_cm_from_gmt_cpt(f"{file_path}/../data/GMT_haxby_new.cpt"))
