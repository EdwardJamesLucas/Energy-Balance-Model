"""Component Module: System Mediator"""

import numpy as np

import components.pv as pv
import components.pvt as pvt
import components.hp as hp
import components.storage as st
import components.house as hs
import components.simulator as sim

from components.helper_funcs import celsius_to_kelvin


class SysMediator:
    """Overarching class that operates the individual calculator classes in accordance with the system layout."""

    def __init__(
        self,
        df,
        storage_calculator: st.StorageCalculator,
        pvt_calculator: pvt.PVTCalculator,
        house_calculator: hs.HouseCalculator,
        component_pv: pv.PV,
        component_hp: hp.HeatPump,
        simulator: sim.Simulator,
    ):
        # Nicknames
        self.sc = storage_calculator
        self.pv = component_pv
        self.hp = component_hp
        self.pvtc = pvt_calculator
        self.hc = house_calculator
        self.sim = simulator

        # System connections
        self.storage_is_included = df.storage_is_included[0]
        self.pvt_is_included = df.pvt_is_included[0]
        self.house_is_included = df.house_is_included[0]
        self.pvt_charges_storage = df.pvt_charges_storage[0]
        self.house_discharges_storage = df.house_discharges_storage[0]
        self.hp_charges_storage = df.hp_charges_storage[0]
        self.connection_pv_hp = df.connection_pv_hp[0]
        self.city = df.city[0]

        # Simulation properties
        self.no_of_time_steps = simulator.no_of_time_steps
        self.repetition_period = simulator.repetition_period
        self.time_step = simulator.time_step
        self.time_steps = np.array([x + 1 for x in range(simulator.no_of_time_steps)])

        # Energy balance validation starting values
        self.validation_storage_eb_losses = 0
        self.validation_storage_eb_discharged = 0
        self.validation_storage_eb_charged = 0

        # Logging
        self.record_sl_temps_init(self.repetition_period)
        self.print_system_layout()

    def print_system_layout(self):
        """Helper function to indicate component interaction. Will be removed in final version."""
        if self.pvt_charges_storage:
            print("PVT is charging storage when able.")
        if self.house_discharges_storage:
            print("House is discharging storage when able.")
        if self.hp_charges_storage:
            print("Heat Pump charges storage when able.")

    def charge_storage_with_pvt(self, tstep):
        """Charges storage if HTF is hotter than top layer temperature; otherwise recirculates HTF"""
        if self.pvtc.check_htf_temp(self.sc.sd.sl_temps[0], tstep):
            self.sc.set_charging_mass(self.pvtc.pvt.massflow)
            self.sc.charge_storage()
            self.pvtc.replenish_htf(self.sc.sd.sl_temps[-1])
            self.sc.calc_layer_enthalpy_change_charging()
        else:
            self.pvtc.reciruclate_htf(tstep)
            self.sc.set_charging_mass()
            self.sc.charge_storage()

    def discharge_storage_with_house(self, tstep):
        """Method to discharge storage using current massflow to, and temperature drop across, House"""
        self.sc.set_discharging_mass(self.hc.hd.massflow[tstep])
        self.sc.calc_discharging_enthalpy(self.hc.hd.temp_drop)
        self.sc.calc_layer_enthalpy_change_discharging()

    def charge_storage_with_hp(self, ambient_temps, tstep):
        """Adjust charging enthalpy with heat from heat pump"""

        if self.sc.check_capacity() < 0.95:
            self.hp.run_heatpump(ambient_temps[tstep], self.sc.sd.sl_temps[-1], self.pv.solar_power[tstep])
            self.sc.set_charging_mass(self.hp.massflow)
            self.sc.charge_storage(self.hp.condenser_temp)
            self.sc.calc_layer_enthalpy_change_charging()

    def energy_balance_storage(self, ambient_temps, tstep):
        """Conduct energy balance across storage to calculate new layer temperatures."""

        # Set default charging and discharging values to 0
        self.sc.sd.sl_enthalpy_change_charging = np.array([0.0 for x in range(self.sc.sd.no_of_layers)])
        self.sc.sd.sl_enthalpy_change_discharging = np.array([0.0 for x in range(self.sc.sd.no_of_layers)])

        # Charging of storage
        if self.pvt_charges_storage:
            self.charge_storage_with_pvt(tstep)
        if self.hp_charges_storage and self.connection_pv_hp:
            self.charge_storage_with_hp(ambient_temps, tstep)

        # Discharging of storage
        if self.house_discharges_storage:
            self.discharge_storage_with_house(tstep)

        # Calculate energy balance
        self.sc.calc_layer_losses(ambient_temps[tstep])
        self.sc.calc_layer_conduction()
        self.sc.calc_layer_enthalpy()
        self.sc.calc_layer_temps()
        self.sc.reorder_layer_temps()

    def energy_balance_pvt(self, direct_irrad, diffuse_irrad):
        """Energy balance across PVT to calculate temperature lift in HTF."""
        self.pvtc.calc_solar_geometry(self.time_steps)
        self.pvtc.calc_solar_irrad(direct_irrad, diffuse_irrad)
        self.pvtc.calc_solar_energy()
        self.pvtc.calc_solar_temp_lift()

    def energy_balance_pv(self, direct_irrad, diffuse_irrad):
        """Energy balance across PV to calculate electrical energy produced."""
        self.pv.calc_solar_geometry(self.time_steps)
        self.pv.calc_solar_irrad(direct_irrad, diffuse_irrad)
        self.pv.calc_solar_energy()
        self.pv.calc_solar_power()

    def energy_balance_house(self, weather_data):
        """Energy balance across house based on heating demand to calculate mass flow of circulated storage media."""
        self.hc.calc_heating(weather_data)
        self.hc.calc_heating_mass_flow(self.time_step)

    def record_sl_temps_init(self, repetition_period):
        """Helper function to initialise arrays for recording storage layer temperatures for plotting."""
        self.sc.record_layer_temps_init(repetition_period)

    def record_sl_temps(self, tstep, repetition_period):
        """Helper function to record storage layer temperatures of current time step for plotting."""
        self.sc.record_layer_temps(tstep, repetition_period)

    def validate_energy_balance(self):
        "Temporary helper method to keep track of energy balance"
        self.validation_storage_eb_losses += sum(self.sc.sd.sl_enthalpy_losses)
        self.validation_storage_eb_charged += sum(self.sc.sd.sl_enthalpy_change_charging)
        self.validation_storage_eb_discharged += sum(self.sc.sd.sl_enthalpy_change_discharging)
