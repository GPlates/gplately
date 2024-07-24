#!/usr/bin/env python3


from common import *

import gplately

print(gplately.__file__)
from gplately.utils import crustal_production

if __name__ == "__main__":
    # https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc
    for i in range(3, 10):
        crustal_production_rate = crustal_production.compute_crustal_production(
            "Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc", i
        )
        print(f"{i}: {crustal_production_rate}")
