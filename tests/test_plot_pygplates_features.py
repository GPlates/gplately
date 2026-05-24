import importlib.util
import sys
import types
from pathlib import Path


class FakeGeoDataFrame:
    def __init__(self, data, geometry=None, crs=None):
        self.geometry = data.get("geometry", [])
        self.crs = crs

    def __len__(self):
        return len(self.geometry)

    def plot(self, ax=None, **kwargs):
        return ax

    def to_crs(self, projection):
        self.crs = projection
        return self


class FakeFeature:
    def __init__(self, valid_end, geometries):
        self._valid_end = valid_end
        self._geometries = geometries

    def get_valid_time(self, _):
        return (None, self._valid_end)

    def get_geometries(self):
        return self._geometries


class FakeAxis:
    def __init__(self):
        self.calls = []

    def add_geometries(self, geometries, crs=None, **kwargs):
        self.calls.append({"geometries": geometries, "crs": crs, "kwargs": kwargs})


def _load_modules(monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]

    gplately_mod = types.ModuleType("gplately")
    gplately_mod.__path__ = []
    mapping_mod = types.ModuleType("gplately.mapping")
    mapping_mod.__path__ = []

    utils_mod = types.ModuleType("gplately.utils")
    utils_mod.__path__ = []
    plot_utils_mod = types.ModuleType("gplately.utils.plot_utils")
    plot_utils_mod.plot_subduction_teeth = lambda *args, **kwargs: None

    tools_mod = types.ModuleType("gplately.tools")
    tools_mod.EARTH_RADIUS = 6371

    grids_mod = types.ModuleType("gplately.grids")
    grids_mod.Raster = type("Raster", (), {})

    geometry_mod = types.ModuleType("gplately.geometry")
    geometry_mod.pygplates_to_shapely = lambda geometry: f"shapely({geometry})"

    geopandas_mod = types.ModuleType("geopandas")
    geodataframe_mod = types.ModuleType("geopandas.geodataframe")
    geodataframe_mod.GeoDataFrame = FakeGeoDataFrame
    geopandas_mod.geodataframe = geodataframe_mod

    class PlateCarree:
        pass

    class Robinson:
        def __init__(self, central_longitude=None):
            self.central_longitude = central_longitude

    ccrs_mod = types.ModuleType("cartopy.crs")
    ccrs_mod.PlateCarree = PlateCarree
    ccrs_mod.Robinson = Robinson
    cartopy_mod = types.ModuleType("cartopy")
    cartopy_mod.crs = ccrs_mod

    distant_future = object()

    class GeoTimeInstant:
        @staticmethod
        def create_distant_future():
            return distant_future

    pygplates_mod = types.ModuleType("pygplates")
    pygplates_mod.Feature = FakeFeature
    pygplates_mod.FeatureCollection = list
    pygplates_mod.GeoTimeInstant = GeoTimeInstant

    pygmt_mod = types.ModuleType("pygmt")
    pygmt_mod.config = lambda **kwargs: None

    monkeypatch.setitem(sys.modules, "gplately", gplately_mod)
    monkeypatch.setitem(sys.modules, "gplately.mapping", mapping_mod)
    monkeypatch.setitem(sys.modules, "gplately.utils", utils_mod)
    monkeypatch.setitem(sys.modules, "gplately.utils.plot_utils", plot_utils_mod)
    monkeypatch.setitem(sys.modules, "gplately.tools", tools_mod)
    monkeypatch.setitem(sys.modules, "gplately.grids", grids_mod)
    monkeypatch.setitem(sys.modules, "gplately.geometry", geometry_mod)
    monkeypatch.setitem(sys.modules, "geopandas", geopandas_mod)
    monkeypatch.setitem(sys.modules, "geopandas.geodataframe", geodataframe_mod)
    monkeypatch.setitem(sys.modules, "cartopy", cartopy_mod)
    monkeypatch.setitem(sys.modules, "cartopy.crs", ccrs_mod)
    monkeypatch.setitem(sys.modules, "pygplates", pygplates_mod)
    monkeypatch.setitem(sys.modules, "pygmt", pygmt_mod)

    modules = {}
    for name, relative_path in (
        ("gplately.mapping.plot_engine", "gplately/mapping/plot_engine.py"),
        ("gplately.mapping.cartopy_plot", "gplately/mapping/cartopy_plot.py"),
        ("gplately.mapping.pygmt_plot", "gplately/mapping/pygmt_plot.py"),
    ):
        spec = importlib.util.spec_from_file_location(name, repo_root / relative_path)
        module = importlib.util.module_from_spec(spec)
        monkeypatch.setitem(sys.modules, name, module)
        spec.loader.exec_module(module)
        modules[name] = module

    return modules, distant_future, PlateCarree


def test_cartopy_plot_pygplates_features_filters_and_plots(monkeypatch):
    modules, distant_future, plate_carree_cls = _load_modules(monkeypatch)

    engine = modules["gplately.mapping.cartopy_plot"].CartopyPlotEngine()
    ax = FakeAxis()

    valid_feature = FakeFeature(0, ["geometry-a"])
    future_feature = FakeFeature(10, ["geometry-b"])
    distant_feature = FakeFeature(distant_future, ["geometry-c"])

    result = engine.plot_pygplates_features(ax, [valid_feature, future_feature, distant_feature])

    assert result is ax
    assert len(ax.calls) == 2
    assert ax.calls[0]["geometries"] == ["shapely(geometry-a)"]
    assert ax.calls[1]["geometries"] == ["shapely(geometry-c)"]
    assert isinstance(ax.calls[0]["crs"], plate_carree_cls)
    assert ax.calls[0]["kwargs"]["edgecolor"] == "blue"
    assert ax.calls[0]["kwargs"]["facecolor"] == "none"


def test_pygmt_plot_pygplates_features_uses_geodataframe_path(monkeypatch):
    modules, _, _ = _load_modules(monkeypatch)

    engine = modules["gplately.mapping.pygmt_plot"].PygmtPlotEngine()
    captured = {}

    def fake_plot_geo_data_frame(ax_or_fig, gdf, **kwargs):
        captured["ax_or_fig"] = ax_or_fig
        captured["gdf"] = gdf
        captured["kwargs"] = kwargs
        return "plotted"

    engine.plot_geo_data_frame = fake_plot_geo_data_frame

    valid_feature = FakeFeature(0, ["geometry-a"])
    invalid_feature = FakeFeature(5, ["geometry-b"])

    result = engine.plot_pygplates_features("figure", [valid_feature, invalid_feature], edgecolor="red")

    assert result == "plotted"
    assert captured["ax_or_fig"] == "figure"
    assert len(captured["gdf"]) == 1
    assert captured["gdf"].geometry == ["shapely(geometry-a)"]
    assert captured["kwargs"]["edgecolor"] == "red"


def test_pygmt_plot_pygplates_features_returns_input_when_no_geometry(monkeypatch):
    modules, _, _ = _load_modules(monkeypatch)

    engine = modules["gplately.mapping.pygmt_plot"].PygmtPlotEngine()
    invalid_feature = FakeFeature(5, ["geometry-b"])

    result = engine.plot_pygplates_features("figure", [invalid_feature])

    assert result == "figure"
