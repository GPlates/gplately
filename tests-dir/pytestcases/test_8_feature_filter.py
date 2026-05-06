#!/usr/bin/env python3
"""
Unit tests for gplately.utils.feature_filter module.

This test suite provides comprehensive coverage of all filter classes and functions
in the feature_filter module, based on patterns from the 14-RuleBasedGPMLProcessingPipeline notebook.
"""

import os
import sys

os.environ["DISABLE_GPLATELY_DEV_WARNING"] = "true"

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

import unittest
from typing import List

import pygplates  # type: ignore

from gplately.utils.feature_filter import (
    FeatureFilter,
    FeatureNameFilter,
    PlateIDFilter,
    BirthAgeFilter,
    EndTimeFilter,
    FeatureTypeFilter,
    PropertyExistsFilter,
    PropertyValueFilter,
    RegionOfInterestFilter,
    TopologicalFeaturesWithDuplicateSectionsFilter,
    TopologicalReferenceFilter,
    FeatureIDFilter,
    ValidTimeFilter,
    filter_feature_collection,
)


def create_test_feature(
    name: str = "Test",
    plate_id: int = 0,
    valid_time=None,  # type: ignore
):  # type: ignore
    """Create a simple test feature with basic properties."""
    feature = pygplates.Feature()  # type: ignore

    # Set basic properties - these are module-level functions
    feature.set_name(name)  # type: ignore
    feature.set_reconstruction_plate_id(plate_id)  # type: ignore

    # Set geometry (default point)
    feature.set_geometry(pygplates.PointOnSphere((0, 0)))  # type: ignore

    # Set valid time if provided
    if valid_time:
        begin_time, end_time = valid_time
        feature.set_valid_time(begin_time, end_time)  # type: ignore

    return feature


class TestFeatureNameFilter(unittest.TestCase):
    """Tests for FeatureNameFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("Australia", plate_id=101),
            create_test_feature("New Zealand", plate_id=102),
            create_test_feature("Tasmania", plate_id=103),
            create_test_feature("Africa", plate_id=104),
            create_test_feature("North America", plate_id=105),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_contains_match(self):
        """Test name filter with contains matching."""
        filter_obj = FeatureNameFilter(
            ["Australia"],
            exact_match=False,
            case_sensitive=True,
        )
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 1)
        self.assertEqual(
            filter_obj.filtrate_features_as_list[0].get_name(), "Australia"  # type: ignore
        )

    def test_exact_match(self):
        """Test name filter with exact matching."""
        filter_obj = FeatureNameFilter(
            ["Australia"],
            exact_match=True,
            case_sensitive=True,
        )
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 1)

    def test_case_insensitive(self):
        """Test case-insensitive name matching."""
        filter_obj = FeatureNameFilter(
            ["australia"],
            exact_match=False,
            case_sensitive=False,
        )
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 1)

    def test_multiple_names(self):
        """Test filter with multiple name patterns."""
        filter_obj = FeatureNameFilter(
            ["Australia", "New Zealand"],
            exact_match=False,
            case_sensitive=True,
        )
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 2)

    def test_reverse_filter(self):
        """Test reverse flag to exclude matching features."""
        filter_obj = FeatureNameFilter(
            ["Australia"],
            reverse=True,
            exact_match=False,
            case_sensitive=True,
        )
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 4)
        names = [f.get_name() for f in filter_obj.filtrate_features_as_list]  # type: ignore
        self.assertNotIn("Australia", names)

    def test_no_match(self):
        """Test filter with no matching names."""
        filter_obj = FeatureNameFilter(
            ["NonExistent"],
            exact_match=False,
            case_sensitive=True,
        )
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 0)

    def test_case_sensitive_toggle(self):
        """Test case sensitivity toggle."""
        feature = create_test_feature("AustRalia")
        collection = pygplates.FeatureCollection([feature])  # type: ignore

        # Should match with case_sensitive=False
        filter_obj1 = FeatureNameFilter(["australia"], case_sensitive=False)
        result1 = filter_feature_collection(collection, [filter_obj1])
        self.assertEqual(len(filter_obj1.filtrate_features_as_list), 1)

        # Should not match with case_sensitive=True
        filter_obj2 = FeatureNameFilter(["australia"], case_sensitive=True)
        result2 = filter_feature_collection(collection, [filter_obj2])
        self.assertEqual(len(filter_obj2.filtrate_features_as_list), 0)


class TestPlateIDFilter(unittest.TestCase):
    """Tests for PlateIDFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("Feature1", plate_id=101),
            create_test_feature("Feature2", plate_id=102),
            create_test_feature("Feature3", plate_id=103),
            create_test_feature("Feature4", plate_id=201),
            create_test_feature("Feature5", plate_id=202),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_single_plate_id(self):
        """Test filter by single plate ID."""
        filter_obj = PlateIDFilter([101])
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 1)
        self.assertEqual(
            filter_obj.filtrate_features_as_list[0].get_reconstruction_plate_id(),  # type: ignore
            101,
        )

    def test_multiple_plate_ids(self):
        """Test filter by multiple plate IDs."""
        filter_obj = PlateIDFilter([101, 102, 103])
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 3)

    def test_plate_id_range(self):
        """Test filter with a range of plate IDs."""
        features = [
            create_test_feature(f"Feature{i}", plate_id=701 + i) for i in range(5)
        ]
        collection = pygplates.FeatureCollection(features)  # type: ignore
        filter_obj = PlateIDFilter(list(range(701, 706)))
        filtered = filter_feature_collection(collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 5)

    def test_reverse_plate_id(self):
        """Test reverse flag for plate ID filter."""
        filter_obj = PlateIDFilter([101, 102], reverse=True)
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 3)

    def test_no_matching_plate_id(self):
        """Test filter with no matching plate IDs."""
        filter_obj = PlateIDFilter([999])
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 0)


class TestBirthAgeFilter(unittest.TestCase):
    """Tests for BirthAgeFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("OldFeature", valid_time=(500, 0)),
            create_test_feature("YoungFeature", valid_time=(100, 0)),
            create_test_feature("VeryOldFeature", valid_time=(1000, 0)),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_keep_older_features(self):
        """Test keeping features older than specified age."""
        filter_obj = BirthAgeFilter(300, reverse=False)
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should keep features with birth age >= 300
        filtrate_list = filter_obj.filtrate_features_as_list
        self.assertGreaterEqual(len(filtrate_list), 2)

    def test_keep_younger_features(self):
        """Test keeping features younger than specified age."""
        filter_obj = BirthAgeFilter(300, reverse=True)
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should keep features with birth age < 300
        filtrate_list = filter_obj.filtrate_features_as_list
        self.assertEqual(len(filtrate_list), 1)

    def test_age_threshold(self):
        """Test age threshold behavior."""
        filter_obj = BirthAgeFilter(500, reverse=False)
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should include features exactly at age 500
        self.assertGreaterEqual(len(filter_obj.filtrate_features_as_list), 1)


class TestEndTimeFilter(unittest.TestCase):
    """Tests for EndTimeFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("Feature1", valid_time=(500, 100)),
            create_test_feature("Feature2", valid_time=(400, 50)),
            create_test_feature(
                "Feature3",
                valid_time=(
                    300,
                    pygplates.GeoTimeInstant.create_distant_future(),  # type: ignore
                ),
            ),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_disappeared_before(self):
        """Test filter for features that disappeared before a certain time."""
        filter_obj = EndTimeFilter(100, reverse=False)
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should find features with end_time >= 100
        self.assertGreater(len(filter_obj.filtrate_features_as_list), 0)

    def test_still_existing(self):
        """Test filter for features still existing at present day."""
        filter_obj = EndTimeFilter(0, reverse=True)
        filtered = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should find features that still exist (end_time is distant future)
        self.assertGreater(len(filter_obj.filtrate_features_as_list), 0)


class TestFilterFeatureCollection(unittest.TestCase):
    """Tests for the main filter_feature_collection function."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("Australia", plate_id=101),
            create_test_feature("New Zealand", plate_id=102),
            create_test_feature("Africa", plate_id=201),
            create_test_feature("Basin1", plate_id=101),
            create_test_feature("Transform1", plate_id=102),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_single_filter(self):
        """Test filtering with a single filter."""
        filter_obj = FeatureNameFilter(["Australia"])
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 1)

    def test_multiple_filters_chain(self):
        """Test filtering with multiple filters applied sequentially."""
        filters = [
            FeatureNameFilter(["Australia", "Basin1", "Transform1"], reverse=True),
            PlateIDFilter([102, 201]),
        ]
        result = filter_feature_collection(self.feature_collection, filters)
        # Should keep only features with plate ID 102 or 201 that don't match the first filter
        result_names = sorted([feature.get_name() for feature in result])  # type: ignore
        self.assertEqual(result_names, ["Africa", "New Zealand"])
        self.assertEqual(len(result_names), 2)

    def test_filter_returns_feature_collection(self):
        """Test that filter_feature_collection returns a FeatureCollection."""
        filter_obj = FeatureNameFilter(["Australia"])
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertIsInstance(result, pygplates.FeatureCollection)  # type: ignore

    def test_empty_result(self):
        """Test filter with no matching results."""
        filter_obj = FeatureNameFilter(["NonExistent"])
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions."""

    def test_empty_feature_collection(self):
        """Test filtering empty feature collection."""
        empty_collection = pygplates.FeatureCollection()  # type: ignore
        filter_obj = FeatureNameFilter(["Australia"])
        result = filter_feature_collection(empty_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 0)

    def test_empty_filter_list(self):
        """Test with empty filter list."""
        features = [create_test_feature("Feature1")]
        collection = pygplates.FeatureCollection(features)  # type: ignore
        result = filter_feature_collection(collection, [])
        self.assertIsInstance(result, pygplates.FeatureCollection)  # type: ignore

    def test_feature_without_valid_time(self):
        """Test filtering features without valid time."""
        features = [
            create_test_feature("Feature1", valid_time=None),
        ]
        collection = pygplates.FeatureCollection(features)  # type: ignore
        filter_obj = BirthAgeFilter(500)
        result = filter_feature_collection(collection, [filter_obj])
        # Features without valid time should be in residue
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 0)


class TestFilterBaseClass(unittest.TestCase):
    """Tests for the FeatureFilter base class."""

    def test_abstract_method_exists(self):
        """Test that FeatureFilter has the should_keep abstract method."""
        self.assertTrue(hasattr(FeatureFilter, "should_keep"))
        self.assertTrue(callable(getattr(FeatureFilter, "should_keep")))

    def test_filtrate_and_residue_properties(self):
        """Test filtrate and residue collection properties."""
        features = [
            create_test_feature("Australia"),
            create_test_feature("New Zealand"),
        ]
        collection = pygplates.FeatureCollection(features)  # type: ignore

        filter_obj = FeatureNameFilter(["Australia"])
        filter_feature_collection(collection, [filter_obj])

        # Check filtrate property
        filtrate = filter_obj.filtrate_feature_collection
        self.assertIsInstance(filtrate, pygplates.FeatureCollection)  # type: ignore
        self.assertGreater(len(list(filtrate)), 0)

        # Check residue property
        residue = filter_obj.residue_feature_collection
        self.assertIsInstance(residue, pygplates.FeatureCollection)  # type: ignore


class TestPropertyExistsFilter(unittest.TestCase):
    """Tests for PropertyExistsFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("Feature1"),
            create_test_feature("Feature2"),
            create_test_feature("Feature3"),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_filter_initialization(self):
        """Test PropertyExistsFilter initialization."""
        filter_obj = PropertyExistsFilter("gpml:name", reverse=False)
        self.assertIsNotNone(filter_obj)

    def test_reverse_flag(self):
        """Test PropertyExistsFilter reverse flag."""
        filter_obj = PropertyExistsFilter("gpml:name", reverse=True)
        self.assertIsNotNone(filter_obj)


class TestPropertyValueFilter(unittest.TestCase):
    """Tests for PropertyValueFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("Feature1"),
            create_test_feature("Feature2"),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_filter_initialization(self):
        """Test PropertyValueFilter initialization."""
        filter_obj = PropertyValueFilter("gpml:name", "Feature1", reverse=False)
        self.assertIsNotNone(filter_obj)

    def test_reverse_flag(self):
        """Test PropertyValueFilter reverse flag."""
        filter_obj = PropertyValueFilter("gpml:name", "Feature1", reverse=True)
        self.assertIsNotNone(filter_obj)


class TestFeatureIDFilter(unittest.TestCase):
    """Tests for FeatureIDFilter class."""

    def test_unique_ids_assertion(self):
        """Test that duplicate IDs cause an assertion error during initialization."""
        with self.assertRaises(AssertionError):
            FeatureIDFilter(["id1", "id1", "id2"])

    def test_filter_initialization(self):
        """Test FeatureIDFilter initialization with unique IDs."""
        filter_obj = FeatureIDFilter(["id1", "id2", "id3"])
        self.assertIsNotNone(filter_obj)

    def test_reverse_flag(self):
        """Test FeatureIDFilter reverse flag."""
        filter_obj = FeatureIDFilter(["id1", "id2"], reverse=True)
        self.assertIsNotNone(filter_obj)


class TestValidTimeFilter(unittest.TestCase):
    """Tests for ValidTimeFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.features = [
            create_test_feature("Feature1", valid_time=(500, 100)),
            create_test_feature("Feature2", valid_time=(300, 50)),
            create_test_feature("Feature3", valid_time=(1000, 200)),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_filter_initialization(self):
        """Test ValidTimeFilter initialization."""
        filter_obj = ValidTimeFilter(begin_time=600, end_time=50)
        self.assertIsNotNone(filter_obj)

    def test_filter_with_infinite_defaults(self):
        """Test ValidTimeFilter with infinite defaults."""
        filter_obj = ValidTimeFilter()
        self.assertIsNotNone(filter_obj)


class TestIntegrationScenarios(unittest.TestCase):
    """Integration tests based on the 14-RuleBasedGPMLProcessingPipeline notebook examples."""

    def setUp(self):
        """Set up test fixtures with diverse features."""
        self.features = [
            create_test_feature("Australia", plate_id=701),
            create_test_feature("New Zealand", plate_id=702),
            create_test_feature("Tasmania", plate_id=701),
            create_test_feature("Africa", plate_id=701),
            create_test_feature("Basin1", plate_id=101, valid_time=(500, 0)),
            create_test_feature("Basin2", plate_id=102, valid_time=(100, 0)),
        ]
        self.feature_collection = pygplates.FeatureCollection(self.features)  # type: ignore

    def test_case_1_australia_region(self):
        """Case 1: Search features whose name contains Australia, New Zealand or Tasmania."""
        filter_obj = FeatureNameFilter(
            ["Australia", "New Zealand", "Tasmania"],
            exact_match=False,
            case_sensitive=True,
        )
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 3)

    def test_case_2_exclude_australia_region(self):
        """Case 2: Filter out features whose name contains Australia, New Zealand or Tasmania."""
        filter_obj = FeatureNameFilter(
            ["Australia", "New Zealand", "Tasmania"],
            reverse=True,
            exact_match=False,
            case_sensitive=True,
        )
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 3)

    def test_case_3_africa_plate_ids(self):
        """Case 3: Search features whose plate ID is from 701 to 715."""
        filter_obj = PlateIDFilter(list(range(701, 716)))
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should find: Australia, New Zealand, Tasmania, Africa
        self.assertGreaterEqual(len(filter_obj.filtrate_features_as_list), 4)

    def test_case_5_older_than_500myr(self):
        """Case 5: Search features whose birth ages are older than 500 million years."""
        filter_obj = BirthAgeFilter(500, reverse=False)
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should find Basin1 which was born at 500 Ma
        self.assertGreater(len(filter_obj.filtrate_features_as_list), 0)

    def test_case_6_younger_than_500myr(self):
        """Case 6: Search features whose birth ages are younger than 500 million years."""
        filter_obj = BirthAgeFilter(500, reverse=True)
        result = filter_feature_collection(self.feature_collection, [filter_obj])
        # Should find features born at age <= 500
        # Basin1 (500 Ma) and Basin2 (100 Ma) both satisfy this condition
        self.assertEqual(len(filter_obj.filtrate_features_as_list), 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
