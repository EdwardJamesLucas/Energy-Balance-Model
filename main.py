# Energy Balance Model
# Edward James Lucas 2022
# .\.venv\Scripts\activate

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle

from physics.PhysicalConstants import (
    WATER_DENSITY,
    WATER_HEAT_CAPACITY,
    WATER_CONDUCTION_COEFF,
)


def read_variables(excel_file_name: str, excel_sheet_name) -> pd.DataFrame:
    df = pd.read_excel(excel_file_name, excel_sheet_name)
    return df


def celsius_to_kelvin(temp_in_celsius):
    return temp_in_celsius + 273.0


def kelvin_to_celsius(temp_in_kelvin):
    return temp_in_kelvin - 273.0


class StorageData:
    """Storage class for storing storage data from GUI."""

    def __init__(self, df):
        self.tank_height = df.tank_height[0]
        self.tank_footprint = df.tank_footprint[0]
        self.no_of_layers = df.no_of_layers[0]
        self.time_step = df.time_step[0]
        self.no_of_time_steps = df.no_of_time_steps[0]
        self.tank_u_value = df.tank_u_value[0]
        self.charging_massflow = df.charging_massflow[0]
        self.discharging_massflow = df.discharging_massflow[0]
        self.charging_temp = df.charging_temp[0]
        self.discharging_temp_reduction = df.discharging_temp_reduction[0]
        self.discharging_min_return_water_temp = df.discharging_min_return_water_temp[0]
        self.initially_charged = df.initially_charged[0]
        self.init_charging_temp = df.init_charging_temp[0]
        self.ambient_temp = None


class StorageCalculator:
    """Calculator class for performing calculations on storage data."""

    def __init__(self, StorageData):
        self.sd = StorageData
        self.calc_init_dimensions()
        self.calc_init_masses()
        self.calc_init_temps()
        self.calc_init_enthalpies()
        self.calc_losses_area()

    def calc_init_dimensions(self):
        self.sd.sl_height = self.sd.tank_height / self.sd.no_of_layers
        self.sd.sl_diameter = math.sqrt(4 * self.sd.tank_footprint / math.pi)

    def calc_init_masses(self):
        self.sd.sl_mass = np.array(
            [self.sd.sl_height * self.sd.tank_footprint * WATER_DENSITY for x in range(self.sd.no_of_layers)]
        )
        self.sd.charging_mass = self.sd.charging_massflow * self.sd.time_step
        self.sd.discharging_mass = self.sd.discharging_massflow * self.sd.time_step

        if self.sd.charging_mass > self.sd.sl_mass.max():
            print("Charging fills entire storage layer within a single time step, could be problematic")
        if self.sd.discharging_mass > self.sd.sl_mass.max():
            print("Discharging empties entire storage layer within a single time step, could be problematic")

    def calc_init_temps(self):
        self.sd.sl_temps_across_time = np.ones((self.sd.no_of_time_steps, self.sd.no_of_layers))
        if self.sd.initially_charged:
            self.sd.sl_temps = np.array(
                [celsius_to_kelvin(self.sd.init_charging_temp) for x in range(self.sd.no_of_layers)]
            )
        else:
            self.sd.sl_temps = np.array(
                [celsius_to_kelvin(self.sd.discharging_min_return_water_temp) for x in range(self.sd.no_of_layers)]
            )

    def calc_init_enthalpies(self):
        self.sd.sl_enthalpy = (
            np.array([(self.sd.sl_mass[x] * WATER_HEAT_CAPACITY) for x in range(self.sd.no_of_layers)])
            * self.sd.sl_temps
        )
        self.sd.sl_enthalpy_added = np.array([x for x in range(self.sd.no_of_layers)], dtype=np.float64)
        self.sd.sl_enthalpy_rest = np.array([x for x in range(self.sd.no_of_layers)], dtype=np.float64)

    def calc_losses_area(self):
        self.sd.sl_enthalpy_losses_area = np.array(
            [(math.pi * self.sd.sl_diameter * self.sd.sl_height) for x in range(self.sd.no_of_layers)]
        )
        self.sd.sl_enthalpy_losses_area[[0, -1]] += self.sd.tank_footprint

    def is_charging(self, received_enthalpy):
        if received_enthalpy > 0:
            return True

    def is_discharging(self, discharging_mass):
        if discharging_mass > 0:
            return True

    def calc_charging_enthalpy(self, charging_temp_lift: float, charged_by_pvt: bool):
        if charged_by_pvt:
            self.sd.charging_enthalpy = np.array(
                [
                    self.sd.charging_mass
                    * WATER_HEAT_CAPACITY
                    * celsius_to_kelvin(
                        min(
                            kelvin_to_celsius(self.sd.sl_temps[-1]) + charging_temp_lift,
                            self.sd.charging_temp,
                        )
                    )
                ]
            )
        else:
            self.sd.charging_enthalpy = np.array(
                [self.sd.charging_mass * WATER_HEAT_CAPACITY * celsius_to_kelvin(self.sd.charging_temp)]
            )

    def calc_discharging_enthalpy(self, discharging_temp_reduction: float, discharged_by_house: bool):
        if discharged_by_house:
            pass
        else:
            self.sd.discharging_enthalpy = np.array(
                [
                    self.sd.discharging_mass
                    * WATER_HEAT_CAPACITY
                    * max(
                        celsius_to_kelvin(self.sd.discharging_min_return_water_temp),
                        (self.sd.sl_temps[0] - discharging_temp_reduction),
                    )
                ]
            )

    def record_layer_temps(self, tstep):
        self.sd.sl_temps_across_time[tstep, :] = kelvin_to_celsius(self.sd.sl_temps)

    def calc_layer_enthalpy_change_charging(self):
        for i, (av, bv) in enumerate(zip(self.sd.sl_mass, self.sd.sl_temps)):
            self.sd.sl_enthalpy_added[i] = self.sd.charging_mass * bv * WATER_HEAT_CAPACITY
            self.sd.sl_enthalpy_rest[i] = (av - self.sd.charging_mass) * bv * WATER_HEAT_CAPACITY

        self.sd.sl_enthalpy_change_charging = (
            (np.concatenate([self.sd.charging_enthalpy, self.sd.sl_enthalpy_added[:-1]]))
            + self.sd.sl_enthalpy_rest
            - self.sd.sl_enthalpy
        )

    def calc_layer_enthalpy_change_discharging(self):
        for i, (av, bv) in enumerate(zip(self.sd.sl_mass, self.sd.sl_temps)):
            self.sd.sl_enthalpy_added[i] = self.sd.discharging_mass * bv * WATER_HEAT_CAPACITY
            self.sd.sl_enthalpy_rest[i] = (av - self.sd.discharging_mass) * bv * WATER_HEAT_CAPACITY

        self.sd.sl_enthalpy_change_discharging = (
            (np.concatenate([self.sd.sl_enthalpy_added[1:], self.sd.discharging_enthalpy]))
            + self.sd.sl_enthalpy_rest
            - self.sd.sl_enthalpy
        )

    def calc_layer_losses(self, ambient_air_temp):
        self.sd.sl_enthalpy_losses = (
            self.sd.sl_enthalpy_losses_area
            * self.sd.tank_u_value
            * (celsius_to_kelvin(ambient_air_temp) - self.sd.sl_temps)
            * self.sd.time_step
        )

    def calc_layer_conduction(self):
        self.sd.sl_enthalpy_conduction = (
            self.sd.tank_footprint
            * WATER_CONDUCTION_COEFF
            / self.sd.sl_height
            * (
                np.concatenate(
                    [
                        np.diff(self.sd.sl_temps[:2]),
                        self.sd.sl_temps[:-2] + self.sd.sl_temps[2:] - 2 * self.sd.sl_temps[1:-1],
                        np.diff(self.sd.sl_temps[-1:-3:-1]),
                    ]
                )
            )
            * self.sd.time_step
        )

    def calc_layer_enthalpy(self):
        self.sd.sl_enthalpy = (
            self.sd.sl_enthalpy
            + self.sd.sl_enthalpy_change_charging
            + self.sd.sl_enthalpy_change_discharging
            + self.sd.sl_enthalpy_losses
            + self.sd.sl_enthalpy_conduction
        )

    def calc_layer_temps(self):
        self.sd.sl_temps = (self.sd.sl_enthalpy) / (self.sd.sl_mass * WATER_HEAT_CAPACITY)

    def reorder_layer_temps(self):
        self.sd.sl_temps = np.sort(self.sd.sl_temps)[::-1]


class PVTData:
    """Data class for storing PVT data from GUI."""

    def __init__(self, df):
        self.max_efficiency = df.max_efficiency[0]
        self.loss_coeff_1 = df.loss_coeff_1[0]
        self.loss_coeff_2 = df.loss_coeff_2[0]
        self.massflow = df.massflow[0]
        self.time_step = df.time_step[0]
        self.no_of_time_steps = df.no_of_time_steps[0]
        self.inclination_angle_south = df.inclination_angle_south[0]
        self.inclination_angle_west = df.inclination_angle_west[0]
        self.inclination_angle_east = df.inclination_angle_east[0]
        self.surface_area_horiz = df.surface_area_horiz[0]
        self.surface_area_south = df.surface_area_south[0]
        self.surface_area_west = df.surface_area_west[0]
        self.surface_area_east = df.surface_area_east[0]
        self.max_return_temp = df.max_return_temp[0]
        self.collector_temp = df.collector_temp[0]
        self.solar_latitute = df.solar_latitute[0]  # luzern


class PVTCalculator:
    """Calculator class for performing calculations on PVT data."""

    def __init__(self, PVTData) -> PVTData:
        self.pvt = PVTData
        self.calc_total_area()

    def calc_total_area(self):
        self.pvt.surface_area_total = (
            self.pvt.surface_area_south + self.pvt.surface_area_west + self.pvt.surface_area_east
        )

    def calc_solar_declination(self, tstep):
        self.pvt.solar_declination = 23.45 * np.sin(2 * math.pi * (284.5 + np.floor(tstep / 24)) / 365)

    def calc_solar_time(self, tstep):
        self.pvt.solar_time = (tstep - (np.floor(tstep / 24) * 24) - 12.5) * 15

    def calc_solar_elevation(self):
        self.pvt.solar_elevation = np.degrees(
            np.array(
                [
                    max(0, x)
                    for x in (
                        np.arcsin(
                            (
                                np.sin(np.radians(self.pvt.solar_latitute))
                                * np.sin(np.radians(self.pvt.solar_declination))
                            )
                            + (
                                np.cos(np.radians(self.pvt.solar_latitute))
                                * np.cos(np.radians(self.pvt.solar_declination))
                                * np.cos(np.radians(self.pvt.solar_time))
                            )
                        )
                    )
                ]
            )
        )

    def calc_solar_azimuth(self):
        self.pvt.solar_azimuth = np.degrees(
            np.arccos(
                np.array(
                    [
                        min(1, x)
                        for x in (
                            (
                                np.sin(np.radians(self.pvt.solar_elevation))
                                * np.sin(np.radians(self.pvt.solar_latitute))
                                - np.sin(np.radians(self.pvt.solar_declination))
                            )
                            / (
                                np.cos(np.radians(self.pvt.solar_elevation))
                                * np.cos(np.radians(self.pvt.solar_latitute))
                            )
                        )
                    ]
                )
            )
        )

    def calc_solar_geometry(self, timesteps):
        self.calc_solar_declination(timesteps)
        self.calc_solar_time(timesteps)
        self.calc_solar_elevation()
        self.calc_solar_azimuth()

    def calc_solar_irrad(self, IDirNorm, IDiffHorz):
        self.pvt.solar_irrad = (
            self.pvt.surface_area_horiz * (IDiffHorz + IDirNorm * np.sin(np.radians(self.pvt.solar_elevation)))
            + self.pvt.surface_area_south
            * (
                IDiffHorz * 0.5 * (1 + np.cos(np.radians(self.pvt.inclination_angle_south)))
                + IDirNorm
                * np.array(
                    [
                        max(0, x)
                        for x in (
                            np.sin(np.radians(self.pvt.inclination_angle_south))
                            * np.cos(np.radians(self.pvt.solar_elevation))
                            * np.cos(np.radians(self.pvt.solar_azimuth))
                            + np.cos(np.radians(self.pvt.inclination_angle_south))
                            * np.sin(np.radians(self.pvt.solar_elevation))
                        )
                    ]
                )
            )
            + self.pvt.surface_area_west
            * (
                IDiffHorz * 0.5 * (1 + np.cos(np.radians(self.pvt.inclination_angle_west)))
                + IDirNorm
                * np.array(
                    [
                        max(0, x)
                        for x in (
                            np.sin(np.radians(self.pvt.inclination_angle_west))
                            * np.cos(np.radians(self.pvt.solar_elevation))
                            * np.cos(np.radians(self.pvt.solar_azimuth))
                            + np.cos(np.radians(self.pvt.inclination_angle_west))
                            * np.sin(np.radians(self.pvt.solar_elevation))
                        )
                    ]
                )
            )
            + self.pvt.surface_area_east
            * (
                IDiffHorz * 0.5 * (1 + np.cos(np.radians(self.pvt.inclination_angle_east)))
                + IDirNorm
                * np.array(
                    [
                        max(0, x)
                        for x in (
                            np.sin(np.radians(self.pvt.inclination_angle_east))
                            * np.cos(np.radians(self.pvt.solar_elevation))
                            * np.cos(np.radians(self.pvt.solar_azimuth))
                            + np.cos(np.radians(self.pvt.inclination_angle_east))
                            * np.sin(np.radians(self.pvt.solar_elevation))
                        )
                    ]
                )
            )
        )

    def calc_solar_enthalpy(self):
        self.pvt.solar_enthalpy = self.pvt.solar_irrad * self.pvt.time_step  # * self.pvt.max_efficiency

    def calc_solar_temp_lift(self):
        self.pvt.solar_temp_lift = self.pvt.solar_enthalpy / (
            self.pvt.massflow * self.pvt.time_step * WATER_HEAT_CAPACITY
        )
        self.pvt.solar_temp_lift_hourly_avg = np.repeat(
            (np.sum(self.pvt.solar_enthalpy.reshape(365, 24), 1) / 24)
            / (self.pvt.massflow * self.pvt.time_step * WATER_HEAT_CAPACITY),
            24,
        )


class HouseData:
    """Data class for storing house data from GUI."""

    def __init__(self, df):
        self.is_included = df.is_included[0]
        self.heating_demand = df.heating_demand[0]
        self.footprint = df.footprint[0]
        self.target_temp = df.target_temp[0]
        self.massflow = df.massflow[0]


class HouseCalculator:
    """Calculator class for performing calculations on house data."""

    def __init__(self, HouseData):
        self.hd = HouseData

    def calc_heating(self, ambient_air_temps):
        self.hd.temp_grad = np.array([max(0, x) for x in self.hd.target_temp - ambient_air_temps])
        self.hd.summed_temp_grad = self.hd.temp_grad.sum()
        self.hd.heating = (
            self.hd.heating_demand * 1000 * self.hd.footprint * (self.hd.temp_grad / self.hd.summed_temp_grad)
        )
        self.hd.heating_enthalpy = self.hd.heating * 3600

    def calc_heating_temp_drop(self, time_step):
        self.hd.heating_temp_drop = self.hd.heating_enthalpy / (self.hd.massflow * time_step * WATER_HEAT_CAPACITY)


class SysMediator:
    """Overarching class that operates the individual calculator classes in accordance with the system layout."""

    def __init__(self, df, StorageCalculator, PVTCalculator, HouseCalculator):
        self.storage_is_included = df.storage_is_included[0]
        self.pvt_is_included = df.pvt_is_included[0]
        self.house_is_included = df.house_is_included[0]
        self.pvt_charges_storage = df.pvt_charges_storage[0]
        self.house_discharges_storage = df.house_discharges_storage[0]
        self.city = df.city[0]
        self.no_of_time_steps = df.no_of_time_steps[0]
        self.time_step = df.time_step[0]
        self.time_steps = np.array([x + 1 for x in range(self.no_of_time_steps)])

        self.sc = StorageCalculator
        self.pvtc = PVTCalculator
        self.hc = HouseCalculator

        self.print_system_layout()

    def print_system_layout(self):
        if self.pvt_charges_storage:
            print("PVT is charging storage when able")
        if self.house_discharges_storage:
            print("House is discharging storage when able")

    def storage_calculations(self, solar_lift, heating_drop, ambient_temps, tstep):
        self.sc.calc_charging_enthalpy(solar_lift[tstep], self.pvt_charges_storage)
        self.sc.calc_discharging_enthalpy(heating_drop[tstep], self.house_discharges_storage)
        self.sc.calc_layer_enthalpy_change_charging()
        self.sc.calc_layer_enthalpy_change_discharging()
        self.sc.calc_layer_losses(ambient_temps[tstep])
        self.sc.calc_layer_conduction()
        self.sc.calc_layer_enthalpy()
        self.sc.calc_layer_temps()
        self.sc.reorder_layer_temps()
        self.sc.record_layer_temps(tstep)

    def pvt_calculations(self, direct_irrad, diffuse_irrad):
        self.pvtc.calc_solar_geometry(self.time_steps)
        self.pvtc.calc_solar_irrad(direct_irrad, diffuse_irrad)
        self.pvtc.calc_solar_enthalpy()
        self.pvtc.calc_solar_temp_lift()

    def house_calculations(self, weather_data):
        self.hc.calc_heating(weather_data)
        self.hc.calc_heating_temp_drop(self.time_step)


class WeatherData:
    """Data class that loads weather data from specific city"""

    def __init__(self, city):
        self.ambient_air_temps: None
        self.irradiance_direct_norm: None
        self.irradiance_diffuse_horiz: None

        with open("climate_data.pickle", "rb") as file:
            climate_data = pickle.load(file)

        self.ambient_air_temps = np.array(climate_data[city].TAir)
        self.irradiance_direct_norm = np.array(climate_data[city].IDirNorm)
        self.irradiance_diffuse_horiz = np.array(climate_data[city].IDiffHor)


def main():

    house_data = HouseData(read_variables("GUI.xlsx", "house_data"))
    house_calculations = HouseCalculator(house_data)
    hc = house_calculations

    storage_data = StorageData(read_variables("GUI.xlsx", "storage_data"))
    storage_calculations = StorageCalculator(storage_data)
    sc = storage_calculations

    pvt_data = PVTData(read_variables("GUI.xlsx", "pvt_data"))
    pvt_calculations = PVTCalculator(pvt_data)
    pvtc = pvt_calculations

    system_mediator = SysMediator(read_variables("GUI.xlsx", "mediator_data"), sc, pvtc, hc)
    sm = system_mediator

    weather_data = WeatherData(sm.city)
    wea = weather_data

    sm.house_calculations(wea.ambient_air_temps)
    sm.pvt_calculations(wea.irradiance_direct_norm, wea.irradiance_diffuse_horiz)

    for tstep in range(sm.no_of_time_steps):
        sm.storage_calculations(
            pvt_data.solar_temp_lift_hourly_avg, house_data.heating_temp_drop, wea.ambient_air_temps, tstep
        )

    plt.plot(storage_data.sl_temps_across_time)
    plt.title(f"Hourly layer temperatures: {sm.city}")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Layer temperature (degC)")
    # plt.plot(ambient_air_temps)
    plt.show()


if __name__ == "__main__":
    main()
