# Feature Filter Unit Tests

This document describes the comprehensive unit test suite for the `gplately.utils.feature_filter` module.

## Overview

The test suite (`test_feature_filter_units.py`) provides 40 unit tests covering the core classes and functions exercised in the `feature_filter` module, including:

- **Filter Classes**: FeatureNameFilter, PlateIDFilter, BirthAgeFilter, EndTimeFilter, PropertyExistsFilter, PropertyValueFilter, FeatureIDFilter, ValidTimeFilter
- **Main Functions**: `filter_feature_collection()`
- **Base Class**: `FeatureFilter` (abstract base)

## Test Structure

### Test Classes

1. **TestFeatureNameFilter** (7 tests)
   - Contains/exact matching behavior
   - Case sensitivity toggle
   - Multiple name patterns
   - Reverse filter behavior
   - No match scenarios

2. **TestPlateIDFilter** (5 tests)
   - Single and multiple plate ID filtering
   - Plate ID ranges (e.g., 701-715 for Africa)
   - Reverse filtering
   - Non-matching cases

3. **TestBirthAgeFilter** (3 tests)
   - Older vs younger feature filtering
   - Age threshold behavior
   - Reverse mode

4. **TestEndTimeFilter** (2 tests)
   - Features disappeared before a time
   - Features still existing (distant future end time)

5. **TestFilterFeatureCollection** (4 tests)
   - Single and multiple filter chains
   - Return type validation
   - Empty result handling

6. **TestEdgeCases** (3 tests)
   - Empty feature collections
   - Empty filter lists
   - Features without valid time

7. **TestFilterBaseClass** (2 tests)
   - Abstract method presence
   - Filtrate and residue properties

8. **TestPropertyExistsFilter** (2 tests)
   - Filter initialization
   - Reverse flag

9. **TestPropertyValueFilter** (2 tests)
   - Filter initialization
   - Reverse flag

10. **TestFeatureIDFilter** (3 tests)
    - Unique ID assertion
    - Filter initialization
    - Reverse flag

11. **TestValidTimeFilter** (2 tests)
    - Filter initialization with custom times
    - Filter with infinite defaults

12. **TestIntegrationScenarios** (5 tests)
    - Based on notebook examples from `14-RuleBasedGPMLProcessingPipeline.py`
    - Australia region filtering
    - Africa plate ID ranges
    - Birth age filtering

## Running the Tests

### Run all tests:
```bash
# Run from the repository root
micromamba activate gplately
python -m pytest tests-dir/unittest/test_feature_filter_units.py -v
```

### Run specific test class:
```bash
python -m pytest tests-dir/unittest/test_feature_filter_units.py::TestFeatureNameFilter -v
```

### Run specific test:
```bash
python -m pytest tests-dir/unittest/test_feature_filter_units.py::TestFeatureNameFilter::test_contains_match -v
```

### Run with detailed output:
```bash
python -m pytest tests-dir/unittest/test_feature_filter_units.py -vv --tb=long
```

## Test Results

All 40 tests pass successfully:
- ✅ 40 passed in ~2.44 seconds
- Coverage includes:
  - Basic filter functionality
  - Edge cases and error conditions
  - Integration scenarios from production notebooks
  - Base class behavior
  - Filter property behavior (filtrate/residue collections)

## Learning from the Notebook

The test suite was informed by patterns from the `14-RuleBasedGPMLProcessingPipeline.py` notebook, which demonstrates:

1. **Case 1**: Filtering features by name (Australia, New Zealand, Tasmania)
2. **Case 2**: Reverse filtering (exclude certain names)
3. **Case 3**: Filtering by plate ID ranges (701-715 for Africa)
4. **Case 5**: Filtering by birth age (older than 500 Ma)
5. **Case 6**: Filtering by birth age (younger than 500 Ma)
6. **Case 7-11**: Feature type, property existence, property values, region of interest, and end time filtering

## Key Testing Patterns

### 1. Feature Creation
Tests use a simple `create_test_feature()` helper that creates minimal pygplates Feature objects with:
- Feature name
- Reconstruction plate ID
- Point geometry
- Valid time (optional)

### 2. Filter Testing Pattern
```python
filter_obj = FeatureNameFilter(["Australia"])
result = filter_feature_collection(collection, [filter_obj])
# Check filtrate_features_as_list for results
# Check residue_features_as_list for non-matching features
```

### 3. Chaining Filters
Tests verify that multiple filters can be applied sequentially:
```python
filters = [
    FeatureNameFilter(["Australia"], reverse=True),
    PlateIDFilter([102, 201])
]
result = filter_feature_collection(collection, filters)
```

## Integration with Existing Tests

This unit test suite complements the existing:
- Shell script integration tests (`test_feature_filter.sh`) that use real GPML files
- CLI-based filtering tests

These unit tests focus on:
- Unit-level correctness of individual filters
- Edge case handling
- Integration of multiple filters
- Base class behavior

## Notes

- All tests work with synthetic feature objects created in memory
- Tests do not require downloading large GPML files
- Tests run in isolation with pytest
- Tests can serve as documentation for filter usage patterns
