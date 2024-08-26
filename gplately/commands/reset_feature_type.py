import pygplates

new_fc = pygplates.FeatureCollection()
for f in pygplates.FeatureCollection("shapes_continents.gpml"):
    # print(f.get_feature_type())
    new_f = pygplates.Feature(
        pygplates.FeatureType.gpml_unclassified_feature,
        f.get_feature_id(),
        pygplates.VerifyInformationModel.yes,  # in case the user provided a wrong feature type
    )
    for p in f:
        new_f.add(
            p.get_name(),
            p.get_time_dependent_value(),
            pygplates.VerifyInformationModel.no,
        )
    new_fc.add(new_f)

new_fc.write("new_shapes_continents.gpmlz")
