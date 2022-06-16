"""Component Class: Heat Pump"""


class HeatPump:
    """Component Class: Heat Pump"""

    def __init__(self, df, simulator):

        # Component attributes
        self.capacity = df.capacity[0]
        self.second_law_efficiency = df.second_law_efficiency[0]
        self.is_air_source = df.is_air_source[0]
        self.cop = 3

        # Time resolution data from simulator class
        self.time_step = simulator.time_step
        self.no_of_time_steps = simulator.no_of_time_steps

    def calc_cop(self, source_temp_in, sink_temp_in):
        """Method for calculating COP"""  # should this return a COP value instead of altering a cop attribute of the HeatPump object?
        self.cop = self.second_law_efficiency * (sink_temp_in / (sink_temp_in - source_temp_in))

        # Boundary conditions on COP
        if self.cop > 6:
            self.cop = 6
        elif self.cop < 0:
            self.cop = 6
        elif self.cop < 1:
            self.cop = 1

    def electricity_consumed(self):
        """Calculate the amount of electricty (in joules) consumed within the current time step"""
        return (self.capacity / self.cop) * self.time_step

    def run_heatpump(self, electricity_supplied=None):
        """Return heat supplied by heat pump within time step at full load. If eletricity is supplied for partial load, calculate heat using the cop"""
        if electricity_supplied * self.cop > self.capacity:
            return self.capacity * self.time_step
        else:
            return electricity_supplied * self.cop * self.time_step
