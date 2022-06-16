"""Component Module: Photovoltaic Thermal Panel"""

import math
import numpy as np

from components.physical_constants import WATER_HEAT_CAPACITY


class PVTData:
    """Data class for storing PVT data from GUI."""

    def __init__(self, df):
        self.max_efficiency = df.max_efficiency[0]
        self.loss_coeff_1 = df.loss_coeff_1[0]
        self.loss_coeff_2 = df.loss_coeff_2[0]
        self.massflow = df.massflow[0]
        self.inclination_angle_south = df.inclination_angle_south[0]
        self.inclination_angle_west = df.inclination_angle_west[0]
        self.inclination_angle_east = df.inclination_angle_east[0]
        self.surface_area_horiz = df.surface_area_horiz[0]
        self.surface_area_south = df.surface_area_south[0]
        self.surface_area_west = df.surface_area_west[0]
        self.surface_area_east = df.surface_area_east[0]
        self.max_return_temp = df.max_return_temp[0]
        self.collector_temp = df.collector_temp[0]
        self.solar_latitude: float
        self.solar_temp_lift_hourly_avg: np.ndarray
        self.htf_temp: None

        self.time_step: None
        self.no_of_time_steps: None


class PVTCalculator:
    """Calculator class for performing calculations on PVT data."""

    def __init__(self, pvt_data: PVTData, weather, simulator):
        self.pvt = pvt_data
        self.pvt.solar_latitude = weather.solar_latitude
        self.pvt.time_step = simulator.time_step
        self.pvt.no_of_time_steps = simulator.no_of_time_steps
        self.calc_total_area()

    def calc_total_area(self):
        """Calculator method for total area of pvt."""
        self.pvt.surface_area_total = (
            self.pvt.surface_area_south + self.pvt.surface_area_west + self.pvt.surface_area_east
        )

    def calc_solar_declination(self, tstep):
        """Calculator method for solar declination at current timestep."""
        self.pvt.solar_declination = 23.45 * np.sin(2 * math.pi * (284.5 + np.floor(tstep / 24)) / 365)

    def calc_solar_time(self, tstep):
        """Calculator method for solar time (in hours) at current timestep."""
        self.pvt.solar_time = (tstep - (np.floor(tstep / 24) * 24) - 12.5) * 15

    def calc_solar_elevation(self):
        """Calculator method for solar elevation at current timestep."""
        self.pvt.solar_elevation = np.degrees(
            np.array(
                [
                    max(0, x)
                    for x in (
                        np.arcsin(
                            (
                                np.sin(np.radians(self.pvt.solar_latitude))
                                * np.sin(np.radians(self.pvt.solar_declination))
                            )
                            + (
                                np.cos(np.radians(self.pvt.solar_latitude))
                                * np.cos(np.radians(self.pvt.solar_declination))
                                * np.cos(np.radians(self.pvt.solar_time))
                            )
                        )
                    )
                ]
            )
        )

    def calc_solar_azimuth(self):
        """Calculator method for solar azimuth angle at current timestep."""
        self.pvt.solar_azimuth = np.degrees(
            np.arccos(
                np.array(
                    [
                        min(1, x)
                        for x in (
                            (
                                np.sin(np.radians(self.pvt.solar_elevation))
                                * np.sin(np.radians(self.pvt.solar_latitude))
                                - np.sin(np.radians(self.pvt.solar_declination))
                            )
                            / (
                                np.cos(np.radians(self.pvt.solar_elevation))
                                * np.cos(np.radians(self.pvt.solar_latitude))
                            )
                        )
                    ]
                )
            )
        )

    def calc_solar_geometry(self, timesteps):
        """Combination method for solar calculations."""
        self.calc_solar_declination(timesteps)
        self.calc_solar_time(timesteps)
        self.calc_solar_elevation()
        self.calc_solar_azimuth()

    def calc_solar_irrad(self, IDirNorm, IDiffHorz):
        """Calculator method for solar irradiation at current timestep."""
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

    def calc_solar_energy(self):
        """Calculator method for solar energy at current timestep."""
        self.pvt.solar_energy = self.pvt.solar_irrad * self.pvt.time_step  # * self.pvt.max_efficiency

    def calc_solar_temp_lift(self):
        """Calculator method for solar temperature lift of HTF at current timestep."""
        self.pvt.solar_temp_lift = self.pvt.solar_energy / (
            self.pvt.massflow * self.pvt.time_step * WATER_HEAT_CAPACITY
        )
        self.pvt.solar_temp_lift_hourly_avg = np.repeat(
            (np.sum(self.pvt.solar_energy.reshape(365, 24), 1) / 24)
            / (self.pvt.massflow * self.pvt.time_step * WATER_HEAT_CAPACITY),
            24,
        )

    def check_htf_temp(self, top_layer_temp, tstep):
        """Method to check if pvt HTF is hot enough to be passed to the top layer of the storage."""
        return self.pvt.htf_temp + self.pvt.solar_temp_lift_hourly_avg[tstep] >= top_layer_temp

    def replenish_htf(self, bottom_layer_temp):
        """Set HTF temperature to storage bottom layer temperature as HTF is replenished from storage."""
        self.pvt.htf_temp = bottom_layer_temp

    def reciruclate_htf(self, tstep):
        """HTF is heated and recirculated within pvt for next time step."""
        self.pvt.htf_temp += self.pvt.solar_temp_lift_hourly_avg[tstep]
