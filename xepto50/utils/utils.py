import pandas as pd
import numpy as np


def unit_conversion(concentration_unit):
    """
    This function takes a preselected name of concentration unit and convert it to a numerical value
    Args:
        concentration_unit: Preselected concentration units.

    Returns:
        conversion_unit
    """

    # Check concentration_unit
    valid_units = ['Molar', 'Millimolar', 'Micromolar', 'Nanomolar', 'Picomolar']
    if concentration_unit not in valid_units:
        raise ValueError(f"Unknown concentration unit: {concentration_unit}")
    conversion_unit = {'Molar': 1, 'Millimolar': 1e-3, 'Micromolar': 1e-6,
                       'Nanomolar': 1e-9, 'Picomolar': 1e-12}[concentration_unit]
    return conversion_unit


def create_df(drug_concentrations, inhibition_values, inhibition_sem=None, conversion_unit=1e-6):
    """
    To create a df and validate data
    Args:
        drug_concentrations: drug concentrations as array or series
        inhibition_values: inhibition values
        inhibition_sem: Standard error of mean for inhibition values, can be None
        conversion_unit: factor to be used to convert drug concentrations

    Returns:
        df: pandas dataframe with three or four columns 'drug_concentration', 'inhibition', 'log10dose'
         and 'inhibition_sem', and sorted by 'drug_concentration'.
    """

    # Check the length of the series and create a dataframe
    if len(drug_concentrations) != len(inhibition_values):
        raise ValueError("Provided concentration data and corresponding response data does not match.")
    df = pd.DataFrame({
        'drug_concentration': drug_concentrations.astype(float),
        'inhibition': inhibition_values.astype(float)
    })
    if inhibition_sem is not None:
        df['inhibition_sem'] = inhibition_sem.astype(float)
    df['log10dose'] = np.log10(df['drug_concentration'] * conversion_unit)
    # sort df by drug_concentration
    df = df.sort_values(by='drug_concentration')
    # Handle duplicated inhibition values
    if df['inhibition'].duplicated().any():
        increment_values = np.arange(0, len(df) * 0.01, 0.01)
        df['inhibition'] += increment_values
    return df


# Function to identify and remove outliers in each row


def remove_row_outliers(row, threshold=1):
    '''

    Args:
        row: A set of data to be used for finding outliers
        threshold: To a value

    Returns:
        row: The data after removal of outliers.

    '''
    Q1 = row.quantile(0.25)
    Q3 = row.quantile(0.75)
    IQR = Q3 - Q1

    # Identify outliers
    outliers = (row < (Q1 - threshold * IQR)) | (row > (Q3 + threshold * IQR))

    # Replace outliers with NaN
    row[outliers] = np.nan
    return row
