import datetime
from dateutil import tz
import calendar
import numpy as np

def convert_decimal_year_to_numpy_datetime64(decimal_year):
    """
    Convert the weird decimal year format of ZMAP origin to datetime64
    
    Assumes UTC
    """
    # get year = integer fraction of decimal_year
    # NOTE: if decimal year is very close to the next higher integer value,
    # rounding takes place
    year_fraction, year = np.modf(decimal_year)
    startyear_dt = datetime.datetime(int(year), 1, 1)

    # get seconds that have passed in fraction of the current year
    if calendar.isleap(startyear_dt.year):
        year_seconds = year_fraction * 86400.0 * 366
    else:
        year_seconds = year_fraction * 86400.0 * 365

    dt = startyear_dt + datetime.timedelta(seconds=year_seconds)
    return np.datetime64(dt)

def convert_epoch_to_numpy_datetime64(epoch_time):
    """
    Convert epoch time to numpy.datetime64
    
    Assumes UTC
    """
    dt = datetime.datetime.fromtimestamp(ts).replace(tzinfo=tz.tzutc())
    return np.datetime64(dt)
    
def convert_multiple_column_timestamp_to_numpy_datetime64(ts, precision='[s]'):
    """
    Converts multiple column timestamp to single numpy.datetime64 column
    
    precision takes any numpy.datetime64 precision operator
    """
    ts = ts.astype(np.float64)
    ts[np.isnan(ts)] = 0
    ts = ts.astype(int)
    yr = ts.yr
    mo = ts.mo
    dy = ts.dy
    hr = ts.hr
    mi = ts.mi
    sc = ts.sc
    if hr > 23:
        print(hr)
    dt = datetime.datetime(year=yr, month=mo, day=dy, hour=hr, minute=mi, second=sc)
    return np.datetime64(dt, utc=True).astype('datetime64{p}'.format(p=precision))#.astype(str)