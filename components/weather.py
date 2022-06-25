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

    def calc_average_temp_winter(self):
        """Returns average air temperature in winter: average of Q1 and Q4"""
        return np.array([self.ambient_air_temps[6570:8760], self.ambient_air_temps[0:2190]]).mean()

    def calc_average_temp_summer(self):
        """Returns average air temperature in summer: average of Q2 and Q3"""
        return self.ambient_air_temps[2190:6570].mean()
