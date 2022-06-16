"""Component Class: Storage"""

import math
import numpy as np

from components.physical_constants import (
    WATER_DENSITY,
    WATER_HEAT_CAPACITY,
    WATER_CONDUCTION_COEFF,
)

from components.helper_funcs import celsius_to_kelvin, kelvin_to_celsius


class StorageData:
    """Storage class for storing storage data from GUI."""

    def __init__(self, df):
        self.tank_height = df.tank_height[0]
        self.tank_footprint = df.tank_footprint[0]
        self.no_of_layers = int(df.no_of_layers[0])
        self.time_step = None
        self.no_of_time_steps = None
        self.tank_u_value = df.tank_u_value[0]
        self.charging_massflow = df.charging_massflow[0]
        self.discharging_massflow = df.discharging_massflow[0]
        self.charging_temp = df.charging_temp[0]
        self.discharging_temp_reduction = df.discharging_temp_reduction[0]
        self.discharging_min_return_water_temp = df.discharging_min_return_water_temp[0]
        self.initially_charged = df.initially_charged[0]
        self.init_charging_temp = df.init_charging_temp[0]
        self.ambient_temp: None
        self.sl_temps: None
        self.sl_temps_across_time: None
        self.boiler_needed = False


class StorageCalculator:
    """Calculator class for performing calculations on storage data."""

    def __init__(self, storage_data: StorageData, simulator):
        self.sd = storage_data
        self.sd.time_step = simulator.time_step
        self.sd.no_of_time_steps = simulator.no_of_time_steps
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

    def set_charging_mass(self, massflow=0, time_step=0):
        self.sd.charging_mass = massflow * time_step

    def set_discharging_mass(self, massflow=0, time_step=0):
        self.sd.discharging_massflow = massflow
        self.sd.discharging_mass = self.sd.discharging_massflow * time_step

    def charge_storage(self):
        """SHOULD THIS METHOD BE GENERIC? AND ABLE TO SET CHARING ENTHALPY EQUAL TO 0"""
        self.sd.charging_enthalpy = np.array(
            [self.sd.charging_mass * WATER_HEAT_CAPACITY * celsius_to_kelvin(self.sd.charging_temp)]
        )

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

    def calc_discharging_enthalpy(self, discharging_temp_reduction: float):
        if celsius_to_kelvin(self.sd.discharging_min_return_water_temp) > (
            self.sd.sl_temps[0] - discharging_temp_reduction
        ):
            self.sd.boiler_needed = True
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

    def record_layer_temps_init(self, repetition_period):
        self.sd.sl_temps_across_time = np.ones((repetition_period, self.sd.no_of_time_steps, self.sd.no_of_layers))

    def record_layer_temps(self, tstep, repetition_period):
        self.sd.sl_temps_across_time[repetition_period, tstep, :] = kelvin_to_celsius(self.sd.sl_temps)

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
        self.sd.sl_temps = np.array((self.sd.sl_enthalpy) / (self.sd.sl_mass * WATER_HEAT_CAPACITY))

    def reorder_layer_temps(self):
        self.sd.sl_temps = np.sort(self.sd.sl_temps)[::-1]

    def check_capacity(self):
        """Returns fraction of utilised storage capacity"""
        return kelvin_to_celsius(np.mean(self.sd.sl_temps)) / self.sd.charging_temp

    def charge_layers(self, layer_weighting: np.ndarray, charging_enthalpy: float):
        """Add charging enthalpy to layers according to weighting array"""
        if charging_enthalpy:
            self.sd.sl_enthalpy_change_charging += layer_weighting * charging_enthalpy

    def discharge_layers(self, layer_weighting: np.ndarray, discharging_enthalpy: float):
        """Remove discharging enthalpy from layers according to weighting array"""
        self.sd.sl_enthalpy_change_discharging -= layer_weighting * discharging_enthalpy
