"""Helper Module: helper functions"""

import pandas as pd


def read_variables(excel_file_name: str, excel_sheet_name) -> pd.DataFrame:
    """Create pandas dataframe by reading GUI excel spreadsheet"""
    df = pd.read_excel(excel_file_name, excel_sheet_name)
    return df


def celsius_to_kelvin(temp_in_celsius):
    """Convert degrees Celsius to degrees Kelvin"""
    return temp_in_celsius + 273.0


def kelvin_to_celsius(temp_in_kelvin):
    """Convert degrees Kelvin to degrees Celsius"""
    return temp_in_kelvin - 273.0

def joules_to_kilowatthrs(energy_in_joules):
    """Convert joules to kilowatt hours"""
    return energy_in_joules / 3600000
