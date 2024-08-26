import pygplates

new_fc = pygplates.FeatureCollection()
for f in pygplates.FeatureCollection("shapes_continents.gpml"):
    # print(f.get_feature_type())
    new_f = pygplates.Feature(
        pygplates.FeatureType.gpml_unclassified_feature,
        f.get_feature_id(),
        pygplates.VerifyInformationModel.yes,
    )
    for p in f:
        new_f.add(p.get_name(), p.get_value(), pygplates.VerifyInformationModel.yes)
    new_fc.add(new_f)

new_fc.write("new_shapes_continents.gpmlz")
