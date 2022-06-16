"""Component Class: Photovoltaic Panel"""
import math
import numpy as np


class PV:
    """Component Class: Photovoltaic Panel"""

    def __init__(self, df, weather, simulator):

        # component attributes
        self.cell_efficiency = df.cell_efficiency[0]
        self.inclination_angle_south = df.inclination_angle_south[0]
        self.inclination_angle_west = df.inclination_angle_west[0]
        self.inclination_angle_east = df.inclination_angle_east[0]
        self.surface_area_horiz = df.surface_area_horiz[0]
        self.surface_area_south = df.surface_area_south[0]
        self.surface_area_west = df.surface_area_west[0]
        self.surface_area_east = df.surface_area_east[0]
        self.solar_latitude = weather.solar_latitude
        self.solar_declination: np.ndarray
        self.solar_time: np.ndarray
        self.solar_elevation: np.ndarray
        self.solar_azimuth: np.ndarray
        self.solar_irrad: np.ndarray
        self.solar_energy: np.ndarray

        # time resolution data from simulator class
        self.time_step = simulator.time_step
        self.no_of_time_steps = simulator.no_of_time_steps

        # calculate total area from attributes
        self.calc_total_area()

    def calc_total_area(self):
        """total area of pv."""
        self.surface_area_total = self.surface_area_south + self.surface_area_west + self.surface_area_east

    def calc_solar_declination(self, tstep):
        """solar declination at current timestep."""
        self.solar_declination = 23.45 * np.sin(2 * math.pi * (284.5 + np.floor(tstep / 24)) / 365)

    def calc_solar_time(self, tstep):
        """solar time (in hours) at current timestep."""
        self.solar_time = (tstep - (np.floor(tstep / 24) * 24) - 12.5) * 15

    def calc_solar_elevation(self):
        """solar elevation at current timestep."""
        self.solar_elevation = np.degrees(
            np.array(
                [
                    max(0, x)
                    for x in (
                        np.arcsin(
                            (np.sin(np.radians(self.solar_latitude)) * np.sin(np.radians(self.solar_declination)))
                            + (
                                np.cos(np.radians(self.solar_latitude))
                                * np.cos(np.radians(self.solar_declination))
                                * np.cos(np.radians(self.solar_time))
                            )
                        )
                    )
                ]
            )
        )

    def calc_solar_azimuth(self):
        """solar azimuth angle at current timestep."""
        self.solar_azimuth = np.degrees(
            np.arccos(
                np.array(
                    [
                        min(1, x)
                        for x in (
                            (
                                np.sin(np.radians(self.solar_elevation)) * np.sin(np.radians(self.solar_latitude))
                                - np.sin(np.radians(self.solar_declination))
                            )
                            / (np.cos(np.radians(self.solar_elevation)) * np.cos(np.radians(self.solar_latitude)))
                        )
                    ]
                )
            )
        )

    def calc_solar_geometry(self, timesteps):
        """Combining solar calculations."""
        self.calc_solar_declination(timesteps)
        self.calc_solar_time(timesteps)
        self.calc_solar_elevation()
        self.calc_solar_azimuth()

    def calc_solar_irrad(self, IDirNorm, IDiffHorz):
        """solar irradiation at current timestep."""
        self.solar_irrad = (
            self.surface_area_horiz * (IDiffHorz + IDirNorm * np.sin(np.radians(self.solar_elevation)))
            + self.surface_area_south
            * (
                IDiffHorz * 0.5 * (1 + np.cos(np.radians(self.inclination_angle_south)))
                + IDirNorm
                * np.array(
                    [
                        max(0, x)
                        for x in (
                            np.sin(np.radians(self.inclination_angle_south))
                            * np.cos(np.radians(self.solar_elevation))
                            * np.cos(np.radians(self.solar_azimuth))
                            + np.cos(np.radians(self.inclination_angle_south))
                            * np.sin(np.radians(self.solar_elevation))
                        )
                    ]
                )
            )
            + self.surface_area_west
            * (
                IDiffHorz * 0.5 * (1 + np.cos(np.radians(self.inclination_angle_west)))
                + IDirNorm
                * np.array(
                    [
                        max(0, x)
                        for x in (
                            np.sin(np.radians(self.inclination_angle_west))
                            * np.cos(np.radians(self.solar_elevation))
                            * np.cos(np.radians(self.solar_azimuth))
                            + np.cos(np.radians(self.inclination_angle_west))
                            * np.sin(np.radians(self.solar_elevation))
                        )
                    ]
                )
            )
            + self.surface_area_east
            * (
                IDiffHorz * 0.5 * (1 + np.cos(np.radians(self.inclination_angle_east)))
                + IDirNorm
                * np.array(
                    [
                        max(0, x)
                        for x in (
                            np.sin(np.radians(self.inclination_angle_east))
                            * np.cos(np.radians(self.solar_elevation))
                            * np.cos(np.radians(self.solar_azimuth))
                            + np.cos(np.radians(self.inclination_angle_east))
                            * np.sin(np.radians(self.solar_elevation))
                        )
                    ]
                )
            )
        )

    def calc_solar_energy(self):
        """solar enthalpy at current timestep."""
        self.solar_energy = self.solar_irrad * self.time_step * self.cell_efficiency
