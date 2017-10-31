import pandas as pd
from utilities import timestamps

def import_catalog(location, timestamp_column='decimal_year', **kwargs):
    """
    imports column names and returns dataframe with timestamp index
    
    accepts kwargs for pandas.read_csv

    location : ?
    timestamp_column : str
    kwargs : pandas.read_csv kwargs
    """
    # TODO : provide ability to parse header files
    timestamp_conversion = {'decimal_year':timestamps.convert_decimal_year_to_numpy_datetime64
                           ,'epoch_time':timestamps.convert_epoch_to_numpy_datetime64
                           ,'none':None}
    
    df = pd.read_csv(location, **kwargs)
    if timestamp_conversion[timestamp_column] is None:
        pass
    else:
        df['timestamp'] = df[timestamp_column].apply(timestamp_conversion[timestamp_column])
    df = df.set_index('timestamp')
    return df

def export_catalog(dataframe, **kwargs):
    """
    exports data as csv

    dataframe : pandas.DataFrame
    kwargs : pandas.DataFrame.to_csv kwargs
    """
    dataframe.to_csv(**kwargs)