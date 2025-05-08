#!/usr/bin/env python3

# You may use the auxiliary functions as a shortcut to create the `PlateReconstruction` and `PlotTopologies` objects.

from gplately.auxiliary import get_gplot, get_plate_reconstruction

# use the auxiliary function to create a PlateReconstruction object
plate_reconstruction_obj = get_plate_reconstruction("Muller2019")
print(plate_reconstruction_obj)

# use the auxiliary function to create a PlotTopologies object
plot_topologies_obj = get_gplot("Muller2019", time=140)
print(plot_topologies_obj)

# there is a PlateReconstruction object inside the PlotTopologies object.
# so, in most cases, a single get_gplot() call is enough.
# you can retrieve the PlateReconstruction object from a PlotTopologies object later,
# for example
plate_reconstruction_obj_within_plot_topologies_obj = (
    plot_topologies_obj.plate_reconstruction
)
print(plate_reconstruction_obj_within_plot_topologies_obj)
