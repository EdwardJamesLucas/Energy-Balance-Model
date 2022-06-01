import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle

from physics.PhysicalConstants import WATER_DENSITY, WATER_HEAT_CAPACITY, WATER_CONDUCTION_COEFF


def read_variables(excel_file_name: str, excel_sheet_name) -> pd.DataFrame:
    df = pd.read_excel(excel_file_name, excel_sheet_name)
    return df


def celsius_to_kelvin(temp_in_celsius):
    return temp_in_celsius + 273.0


def kelvin_to_celsius(temp_in_kelvin):
    return temp_in_kelvin - 273.0


class StorageData:
    def __init__(self, df):
        self.tank_height = df.tank_height[0]
        self.tank_footprint = df.tank_footprint[0]
        self.no_of_layers = df.no_of_layers[0]
        self.time_step = df.time_step[0]
        self.no_of_time_steps = df.no_of_time_steps[0]
        self.tank_u_value = df.tank_u_value[0]
        self.ambient_temp = df.ambient_temp[0]
        self.charging_massflow = df.charging_massflow[0]
        self.discharging_massflow = df.discharging_massflow[0]
        self.charging_temp = df.charging_temp[0]
        self.discharging_temp_reduction = df.discharging_temp_reduction[0]
        self.discharging_min_return_water_temp = df.discharging_min_return_water_temp[0]
        self.initially_charged = df.initially_charged[0]
        self.init_charging_temp = df.init_charging_temp[0]


class StorageCalculator:
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
            np.array([(self.sd.sl_mass * WATER_HEAT_CAPACITY) for x in range(self.sd.no_of_layers)]) * self.sd.sl_temps
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
            self.sd.charging = True
        else:
            self.sd.charging = False

    def calc_charging_enthalpy(self):
        if self.sd.charging:
            self.sd.charging_enthalpy = np.array(
                [self.sd.charging_mass * WATER_HEAT_CAPACITY * celsius_to_kelvin(self.sd.charging_temp)]
            )
        else:
            self.sd.discharging_enthalpy = np.array(
                [
                    self.sd.discharging_mass
                    * WATER_HEAT_CAPACITY
                    * max(
                        celsius_to_kelvin(self.sd.discharging_min_return_water_temp),
                        (self.sd.sl_temps[0] - self.sd.discharging_temp_reduction),
                    )
                ]
            )

    def record_layer_temps(self, tstep):
        self.sd.sl_temps_across_time[tstep, :] = kelvin_to_celsius(self.sd.sl_temps)

    def calc_layer_enthalpy_after_charging(self):
        for i, (av, bv) in enumerate(zip(self.sd.sl_mass, self.sd.sl_temps)):
            if self.sd.charging:
                self.sd.sl_enthalpy_added[i] = self.sd.charging_mass * bv * WATER_HEAT_CAPACITY
                self.sd.sl_enthalpy_rest[i] = (av - self.sd.charging_mass) * bv * WATER_HEAT_CAPACITY
            else:
                self.sd.sl_enthalpy_added[i] = self.sd.discharging_mass * bv * WATER_HEAT_CAPACITY
                self.sd.sl_enthalpy_rest[i] = (av - self.sd.discharging_mass) * bv * WATER_HEAT_CAPACITY

        if self.sd.charging:
            self.sd.sl_enthalpy = (
                np.concatenate([self.sd.charging_enthalpy, self.sd.sl_enthalpy_added[:-1]])
            ) + self.sd.sl_enthalpy_rest
        else:
            self.sd.sl_enthalpy = (
                np.concatenate([self.sd.sl_enthalpy_added[1:], self.sd.discharging_enthalpy])
            ) + self.sd.sl_enthalpy_rest

    def calc_layer_losses(self):
        self.sd.sl_enthalpy_losses = (
            self.sd.sl_enthalpy_losses_area
            * self.sd.tank_u_value
            * (celsius_to_kelvin(self.sd.ambient_temp) - self.sd.sl_temps)
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

    def calc_layer_temps(self):
        self.sd.sl_temps = (self.sd.sl_enthalpy + self.sd.sl_enthalpy_losses + self.sd.sl_enthalpy_conduction) / (
            self.sd.sl_mass * WATER_HEAT_CAPACITY
        )

    def reorder_layer_temps(self):
        self.sd.sl_temps = np.sort(self.sd.sl_temps)[::-1]

    def grab_air_temp(self, ambient_air_temps, tstep):
        self.sd.ambient_temp = ambient_air_temps[tstep]


class PVTData:
    def __init__(self) -> None:
        self.max_efficiency = 0.75
        self.loss_coeff_1 = 3.37
        self.loss_coeff_2 = 0.0104
        self.massflow = 0.02
        self.time_step = 3600
        self.no_of_time_steps = 8760
        self.inclination_angle_south = 30
        self.inclination_angle_west = 90
        self.inclination_angle_east = 30
        self.surface_area_horiz = 0
        self.surface_area_south = 160
        self.surface_area_west = 0
        self.surface_area_east = 0
        self.max_return_temp = 95
        self.collector_temp = 14
        self.solar_latitute = 46.6


class PVTCalculator:
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
        self.pvt.solar_enthalpy = self.pvt.solar_irrad * self.pvt.time_step #* self.pvt.max_efficiency
        print(self.pvt.solar_enthalpy)

    def calc_solar_temp_lift(self):
        self.pvt.solar_temp_lift = self.pvt.solar_enthalpy / (self.pvt.massflow * self.pvt.time_step * WATER_HEAT_CAPACITY)


class HouseData:
    def __init__(self, df):
        self.is_included = df.is_included[0]
        self.heating_demand = df.heating_demand[0]
        self.footprint = df.footprint[0]
        self.target_temp = df.target_temp[0]


class HouseCalculator:
    def __init__(self, house_data):
        self.hd = house_data

    def calc_heating(self, ambient_air_temps):
        self.hd.temp_grad = np.array([max(0, x) for x in self.hd.target_temp - ambient_air_temps])
        self.hd.summed_temp_grad = self.hd.temp_grad.sum()
        self.hd.heating = (
            self.hd.heating_demand * 1000 * self.hd.footprint * (self.hd.temp_grad / self.hd.summed_temp_grad)
        )
        self.hd.heating_enthalpy = self.hd.heating * 3600


def main():
    # .\.venv\Scripts\activate
    with open("climate_data.pickle", "rb") as file:
        climate_data = pickle.load(file)

    df = read_variables("GUI.xlsx", "storage_data")
    storage_data = StorageData(df)
    storage_calculations = StorageCalculator(storage_data)
    sc = storage_calculations

    ambient_air_temps = np.array(climate_data[df.city[0]].TAir)
    irradiance_direct_norm = np.array(climate_data[df.city[0]].IDirNorm)
    irradiance_diffuse_horiz = np.array(climate_data[df.city[0]].IDiffHor)
    timesteps = np.array([x + 1 for x in range(df.no_of_time_steps[0])])

    pvt_data = PVTData()
    pvt_calculations = PVTCalculator(pvt_data)
    pvtc = pvt_calculations
    pvtc.calc_solar_declination(timesteps)
    pvtc.calc_solar_time(timesteps)
    pvtc.calc_solar_elevation()
    pvtc.calc_solar_azimuth()
    pvtc.calc_solar_irrad(irradiance_direct_norm, irradiance_diffuse_horiz)
    pvtc.calc_solar_enthalpy()
    pvtc.calc_solar_temp_lift()

    house_data = HouseData(read_variables("GUI.xlsx", "house_data"))
    house_calculations = HouseCalculator(house_data)
    hc = house_calculations
    hc.calc_heating(ambient_air_temps)
    print(house_data.heating)

    # dummy variable representing enthalpy passed to(+ve)/from(-ve) storage
    recieved_enthalpy = (
        [0 for x in range(int(storage_data.no_of_time_steps / 4))]
        + [1 for x in range(int(storage_data.no_of_time_steps / 2))]
        + [0 for x in range(int(storage_data.no_of_time_steps / 4))]
    )

    # enthalpies = {"charging": np.zeros(8760), "discharging": np.zeros(8760)}

    for tstep in range(storage_data.no_of_time_steps):
        sc.is_charging(recieved_enthalpy[tstep])
        sc.calc_charging_enthalpy()
        # if storage_data.charging:
        #     enthalpies["charging"][tstep] = storage_data.charging_enthalpy
        # else:
        #     enthalpies["discharging"][tstep] = storage_data.discharging_enthalpy
        # print(storage_data.charging_mass, storage_data.charging_enthalpy, storage_data.charging_temp)
        sc.calc_layer_enthalpy_after_charging()
        sc.grab_air_temp(ambient_air_temps, tstep)
        sc.calc_layer_losses()
        sc.calc_layer_conduction()
        sc.calc_layer_temps()
        sc.reorder_layer_temps()
        sc.record_layer_temps(tstep)

    # print(enthalpies["charging"].sum())
    # print(enthalpies["discharging"].sum())

    print(storage_data.sl_temps_across_time[-1, :].mean())
    plt.plot(storage_data.sl_temps_across_time)
    plt.plot(ambient_air_temps)
    plt.show()


if __name__ == "__main__":
    main()
