import pygplates


class SeafloorDeactivatePoints(
    pygplates.ReconstructedGeometryTimeSpan.DeactivatePoints
):
    def __init__(self):
        super(SeafloorDeactivatePoints, self).__init__()
        # Other initialisation you may want...
        ...

    def deactivate(
        self,
        prev_point,
        prev_location,
        prev_time,
        current_point,
        current_location,
        current_time,
    ):
        for (
            prev_boundary_sub_segment
        ) in prev_location.located_in_resolved_boundary().get_boundary_sub_segments():
            # Use feature-specific collision parameters if found (falling back to global collision parameters).
            if (
                prev_boundary_sub_segment.get_feature().get_feature_type()
                == pygplates.FeatureType.gpml_subduction_zone
            ):
                # parameters for subduction zone
                deactivator = (
                    pygplates.ReconstructedGeometryTimeSpan.DefaultDeactivatePoints(
                        threshold_velocity_delta=0.5,  # cms/yr
                        threshold_distance_to_boundary=10,  # kms/myr
                        deactivate_points_that_fall_outside_a_network=True,
                    )
                )
                # ??????
                # I could not get my head around here...
                # How DefaultDeactivatePoints knows that it should use the parameters against gpml_subduction_zone?
                # DefaultDeactivatePoints could have used the parameters to check against other sub segments and return True?
                # sorry, I might be a bit incoherence here. I just could not understand. Maybe it is because I don't really know how DefaultDeactivatePoints works.
                if deactivator.deactivate(
                    prev_point,
                    prev_location,
                    prev_time,
                    current_point,
                    current_location,
                    current_time,
                ):
                    return True
            else:
                # parameters for other feature types
                deactivator = (
                    pygplates.ReconstructedGeometryTimeSpan.DefaultDeactivatePoints(
                        threshold_velocity_delta=0.789,  # cms/yr
                        threshold_distance_to_boundary=12.345,  # kms/myr
                        deactivate_points_that_fall_outside_a_network=True,
                    )
                )
                if deactivator.deactivate(
                    prev_point,
                    prev_location,
                    prev_time,
                    current_point,
                    current_location,
                    current_time,
                ):
                    return True

        return False
