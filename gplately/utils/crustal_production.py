import pygmt
import xarray as xr

# age threshold
x = 14

# https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc
age_grid = xr.open_dataset("Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc")
data = age_grid.to_dataarray()

data = data.where(data < x)

output_dataframe = pygmt.grdvolume(grid=data, contour=[0], output_type="pandas")
print(output_dataframe)
print(output_dataframe.iloc[0, 1])
