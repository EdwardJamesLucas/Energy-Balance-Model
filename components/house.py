"""Component Module: House"""

import numpy as np
from components.physical_constants import WATER_HEAT_CAPACITY


class HouseData:
    """Data class for storing house data from GUI."""

    def __init__(self, df):
        self.heating_demand = df.heating_demand[0]
        self.footprint = df.footprint[0]
        self.target_temp = df.target_temp[0]
        self.massflow = df.massflow[0]
        self.temp_drop = df.temp_drop[0]
        self.autarky_thermal = df.autarky_thermal[0]


class HouseCalculator:
    """Calculator class for performing calculations on house data."""

    def __init__(self, house_data: HouseData, simulator):
        self.hd = house_data
        self.hd.time_step = simulator.time_step

    def calc_heating(self, ambient_air_temps):
        """Calculator method for determining required heating enthalpy based on thermal autarky."""
        self.hd.temp_grad = np.array([max(0, x) for x in self.hd.target_temp - ambient_air_temps])
        self.hd.summed_temp_grad = self.hd.temp_grad.sum()
        self.hd.heating = (
            self.hd.heating_demand * 1000 * self.hd.footprint * (self.hd.temp_grad / self.hd.summed_temp_grad)
        )
        self.hd.heating_enthalpy = (
            self.hd.heating * self.hd.time_step * self.hd.autarky_thermal
        )  # here should be *1000 to get into joules

    def calc_heating_mass_flow(self, time_step):
        """Calculator method for determining massflow of circulated HTF."""
        self.hd.massflow = self.hd.heating_enthalpy / (self.hd.temp_drop * time_step * WATER_HEAT_CAPACITY)
