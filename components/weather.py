"""Weather Module"""

import numpy as np
import pickle


class WeatherData:
    """Data class that loads weather data from specific city"""

    def __init__(self, df):
        self.ambient_air_temps: None
        self.irradiance_direct_norm: None
        self.irradiance_diffuse_horiz: None
        self.solar_latitude = df.solar_latitude[0]
        self.city = df.city[0]

        with open("data/climate_data.pickle", "rb") as file:
            climate_data = pickle.load(file)

        self.ambient_air_temps = np.array(climate_data[self.city].TAir)
        self.irradiance_direct_norm = np.array(climate_data[self.city].IDirNorm)
        self.irradiance_diffuse_horiz = np.array(climate_data[self.city].IDiffHor)
