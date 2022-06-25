"""Component Class: Heat Pump"""
from components.physical_constants import WATER_HEAT_CAPACITY
from components.physical_constants import AIR_HEAT_CAPACITY

from components.helper_funcs import celsius_to_kelvin


class HeatPump:
    """Component Class: Heat Pump"""

    def __init__(self, df, simulator):

        # Component attributes
        self.condenser_duty = df.condenser_duty[0]
        self.compressor_duty = df.compressor_duty[0]
        self.evaporator_duty = self.condenser_duty - self.compressor_duty
        self.load_fraction: float
        self.second_law_efficiency = df.second_law_efficiency[0]
        self.part_load_efficiency = df.part_load_efficiency[0]
        self.is_air_source = df.is_air_source[0]
        self.massflow_air = df.massflow_air[0]
        self.massflow_max = df.massflow_max[0]
        self.massflow_min = df.massflow_min[0]
        self.massflow: float
        self.max_cop = df.max_cop[0]
        self.base_cop = df.base_cop[0]
        self.cop = self.base_cop  # default value
        self.temp_lift = df.temp_lift[0]
        self.min_approach_temp = df.min_approach_temp[0]
        self.evaporator_temp_drop: float
        self.condenser_temp: float
        self.evaporator_temp: float
        self.condenser_temp_max = celsius_to_kelvin(df.condenser_temp_max[0])
        self.evaporator_temp_max = celsius_to_kelvin(df.evaporator_temp_max[0])
        self.evaporator_temp_min = celsius_to_kelvin(df.evaporator_temp_min[0])

        # Time resolution data from simulator class
        self.time_step = simulator.time_step
        self.no_of_time_steps = simulator.no_of_time_steps

    def run_heatpump(self, source_temp_in, sink_temp_in, power_supplied=None) -> float:
        """Calculate outlet temperature and massflow of storage media through heat pump within time step at full load. If eletricity is supplied for partial load, calculate these values for partial load"""
        self.calc_evaporator_duty(power_supplied)
        self.calc_evaporator_temp_drop()
        self.calc_evaporator_temp(source_temp_in)
        self.calc_condenser_temp()

        if self.is_condenser_hot_enough(sink_temp_in) and self.check_evaporator_temp():
            self.calc_massflowrate(sink_temp_in)
        else:
            self.massflow = 0

    def check_evaporator_temp(self):
        """Returns True if the evaporator temp is within its safe operating limits"""
        return self.evaporator_temp <= self.evaporator_temp_max and self.evaporator_temp >= self.evaporator_temp_min

    def calc_evaporator_duty(self, power_supplied=None) -> float:
        """Calculate the heatpump evaporator duty at partial duty based on the available electricity from the PV panel. If no electricity value is provided, this method returns the evaporator duty based on full load"""
        self.load_fraction = min(1, power_supplied / self.compressor_duty)
        if power_supplied and self.load_fraction < 1:
            self.evaporator_duty = (
                (self.condenser_duty / self.compressor_duty) * self.part_load_efficiency - 1
            ) * power_supplied
        elif power_supplied and self.load_fraction == 1:
            self.evaporator_duty = self.condenser_duty - self.compressor_duty
        else:
            self.load_fraction = 1
            self.evaporator_duty = self.condenser_duty - self.compressor_duty

    def calc_evaporator_temp_drop(self) -> float:
        """Calculate the source temperature drop across the evaporator"""
        if self.is_air_source:
            self.evaporator_temp_drop = self.evaporator_duty / (AIR_HEAT_CAPACITY * self.massflow_air)
        else:
            self.evaporator_temp_drop = 0

    def calc_evaporator_temp(self, source_temp_in) -> float:
        """Calculate the required evaporator temperature to extract all the evaporator duty from the source at its given temperature"""
        self.evaporator_temp = celsius_to_kelvin(source_temp_in) - self.evaporator_temp_drop - self.min_approach_temp

    def calc_condenser_temp(self) -> float:
        """Calculate the condenser temperature based on the current evaporator temperature and the temperature lift"""
        self.condenser_temp = min(self.evaporator_temp + self.temp_lift, self.condenser_temp_max)

    def is_condenser_hot_enough(self, sink_temp_in) -> bool:
        """Check if condenser can pass heat to sink"""
        return self.condenser_temp >= sink_temp_in

    def calc_massflowrate(self, sink_temp_in) -> float:
        """Calculate the mass flowrate through the heat pump on the storage side, based on the condenser load, temperature, and the sink inlet temperature. Value is subject to boundary conditions set by massflow_max and massflow_min"""

        self.massflow = (self.evaporator_duty + self.load_fraction * self.compressor_duty) / (
            WATER_HEAT_CAPACITY * (self.condenser_temp - sink_temp_in)
        )

        if self.massflow > self.massflow_max:
            self.massflow = self.massflow_max
        if self.massflow < self.massflow_min:
            self.massflow = 0

    def calc_cop(self) -> float:
        """Calculate the real COP of the Heat Pump. Return value is subject to boundary conditions set by max_cop."""  # should this return a COP value instead of altering a cop attribute of the HeatPump object?
        self.cop = self.second_law_efficiency * (self.condenser_temp / (self.condenser_temp - self.evaporator_temp))

        # Boundary conditions on COP
        if self.cop > self.max_cop:
            self.cop = self.max_cop
        elif self.cop < 0:
            self.cop = self.max_cop
        elif self.cop < 1:
            self.cop = 1

    def electricity_consumed(self) -> float:
        """Calculate the amount of electrical energy (in joules) consumed within the current time step"""
        return (self.condenser_duty / self.cop) * self.time_step
