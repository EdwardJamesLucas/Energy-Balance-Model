"""Energy Balance Model"""
# Edward James Lucas 2022

import numpy as np
import matplotlib.pyplot as plt

import components.pv as pv
import components.pvt as pvt
import components.hp as hp
import components.storage as st
import components.house as hs
import components.sys_mediator as sysm
import components.weather as wth
import components.simulator as sim

from components.helper_funcs import read_variables, celsius_to_kelvin, kelvin_to_celsius


def main():
    """Main Function"""

    # Pulling the data and setting up the data and calculator objects.
    simulator = sim.Simulator(read_variables("GUI.xlsx", "simulator_data"))
    weather_data = wth.WeatherData(read_variables("GUI.xlsx", "weather_data"))

    house_data = hs.HouseData(read_variables("GUI.xlsx", "house_data"))
    house_calculations = hs.HouseCalculator(house_data, simulator)

    storage_data = st.StorageData(read_variables("GUI.xlsx", "storage_data"))
    storage_calculations = st.StorageCalculator(storage_data, simulator)

    pvt_data = pvt.PVTData(read_variables("GUI.xlsx", "pvt_data"))
    pvt_calculations = pvt.PVTCalculator(pvt_data, weather_data, simulator)

    component_pv = pv.PV(read_variables("GUI.xlsx", "pv_data"), weather_data, simulator)
    component_hp = hp.HeatPump(read_variables("GUI.xlsx", "hp_data"), simulator)

    # note to self: will need to reformulate this to account for optional arguments
    system_mediator = sysm.SysMediator(
        read_variables("GUI.xlsx", "mediator_data"),
        storage_calculations,
        pvt_calculations,
        house_calculations,
        component_pv,
        component_hp,
        simulator,
    )

    # "Fill" PVT Panel from Storage bottom layer if required
    if system_mediator.pvt_charges_storage:
        system_mediator.pvtc.replenish_htf(storage_data.sl_temps[-1])

    # component energy balances which can be pre-computed for each timestep
    system_mediator.energy_balance_pvt(weather_data.irradiance_direct_norm, weather_data.irradiance_diffuse_horiz)
    system_mediator.energy_balance_pv(weather_data.irradiance_direct_norm, weather_data.irradiance_diffuse_horiz)
    system_mediator.energy_balance_house(weather_data.ambient_air_temps)
    print(f"Thermal Autarky set at: {house_data.autarky_thermal*100} %")
    print(f" Max massflow between house and storage: {round(max(house_data.massflow), 4)} kg/s")

    track_hp_elec = []
    track_cop = []

    # Time series simulation
    for repetition_period in range(system_mediator.repetition_period):
        for tstep in range(system_mediator.no_of_time_steps):

            if repetition_period == 0:
                track_hp_elec.append(system_mediator.hp.electricity_consumed())
                track_cop.append(system_mediator.hp.cop)

            # Storage
            system_mediator.energy_balance_storage(weather_data.ambient_air_temps, tstep)
            system_mediator.record_sl_temps(tstep, repetition_period)
            system_mediator.validate_energy_balance()
            # if tstep == 1721:
            #     break

    # Plots and print-outs
    show_plot = True
    if show_plot:
        values = np.zeros(
            (system_mediator.no_of_time_steps * system_mediator.repetition_period, storage_data.no_of_layers)
        )
        for repetition_period in range(system_mediator.repetition_period):
            values[
                (0 + system_mediator.no_of_time_steps * repetition_period) : (
                    system_mediator.no_of_time_steps + system_mediator.no_of_time_steps * (repetition_period)
                ),
                :,
            ] = storage_data.sl_temps_across_time[
                repetition_period,
                :,
                :,
            ]

        plt.plot(values)
        plt.title(f"Hourly layer temperatures: {system_mediator.city}")
        plt.xlabel("Time (hrs)")
        plt.ylabel("Layer temperature (degC)")
        plt.show()
        print(
            f"{round(system_mediator.validation_storage_eb_losses/3600000)} kWh Storage losses.",
            f"{round(system_mediator.validation_storage_eb_discharged/3600000)} kWh discharged.",
            f"{round(system_mediator.validation_storage_eb_charged/3600000)} kWh charged.",
            f"{round((system_mediator.validation_storage_eb_charged+system_mediator.validation_storage_eb_losses+system_mediator.validation_storage_eb_discharged)/3600000)} kWh Difference.",
            f"Average Tank Temp: {round(np.mean(storage_data.sl_temps_across_time),1)}",
            f"Boiler required: {storage_data.boiler_needed}",
        )

    plt.plot(system_mediator.pv.solar_energy)
    plt.plot(track_hp_elec)
    plt.show()

    plt.plot(track_cop)
    plt.show()


if __name__ == "__main__":
    main()
