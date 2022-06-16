"""Helper Module: helper functions"""

import pandas as pd


def read_variables(excel_file_name: str, excel_sheet_name) -> pd.DataFrame:
    """Helper function to read in pandas dataframe from GUI excel spreadsheet"""
    df = pd.read_excel(excel_file_name, excel_sheet_name)
    return df


def celsius_to_kelvin(temp_in_celsius):
    """Helper function to convert C to K"""
    return temp_in_celsius + 273.0


def kelvin_to_celsius(temp_in_kelvin):
    """Helper function to convert K to C"""
    return temp_in_kelvin - 273.0
