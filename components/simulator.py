"""Component Module: Simulator"""


class Simulator:
    """Class for simulation time resolution properties"""

    def __init__(self, df):
        self.time_step = df.time_step[0]
        self.no_of_time_steps = df.no_of_time_steps[0]
        self.repetition_period = df.repetition_period[0]
