
def shoot(lon, lat, azimuth, maxdist=None):
    """Shooter Function
    Original javascript on http://williams.best.vwh.net/gccalc.htm
    Translated to python by Thomas Lecocq

    lon : float
    lat : float
    azimuth : float
    return : list
    """
    glat1 = lat * np.pi / 180.
    glon1 = lon * np.pi / 180.
    s = maxdist / 1.852
    faz = azimuth * np.pi / 180.
 
    EPS= 0.00000000005
    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
        alert("Only N-S courses are meaningful, starting at a pole!")
 
    a=6378.13/1.852
    f=1/298.257223563
    r = 1 - f
    tu = r * np.tan(glat1)
    sf = np.sin(faz)
    cf = np.cos(faz)
    if (cf==0):
        b=0.
    else:
        b=2. * np.arctan2 (tu, cf)
 
    cu = 1. / np.sqrt(1 + tu * tu)
    su = tu * cu
    sa = cu * sf
    c2a = 1 - sa * sa
    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    x = (x - 2.) / x
    c = 1. - x
    c = (x * x / 4. + 1.) / c
    d = (0.375 * x * x - 1.) * x
    tu = s / (r * a * c)
    y = tu
    c = y + 1
    while (np.abs (y - c) > EPS):
 
        sy = np.sin(y)
        cy = np.cos(y)
        cz = np.cos(b + y)
        e = 2. * cz * cz - 1.
        c = y
        x = e * cy
        y = e + e - 1.
        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
              d / 4. - cz) * sy * d + tu
 
    b = cu * cy * cf - su * sy
    c = r * np.sqrt(sa * sa + b * b)
    d = su * cy + cu * sy * cf
    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    c = cu * cy - su * sy * cf
    x = np.arctan2(sy * sf, c)
    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    d = ((e * cy * c + cz) * sy * c + y) * sa
    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
 
    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
 
    glon2 *= 180./np.pi
    glat2 *= 180./np.pi
    baz *= 180./np.pi
 
    return (glon2, glat2, baz)
 
def equi(ax, m, centerlon, centerlat, radius, *args, **kwargs):
    """
    plots circle on matplotlib basemap map

    m : mpl_toolkits.Basemap
    centerlon : float
    centerlat : float
    centerlat : float
    radius : float
    type args : args for matplotlib.Axes
    type kwargs : kwargs for matplotlib.Axes
    return: None
    """
    glon1 = centerlon
    glat1 = centerlat
    X = []
    Y = []
    for azimuth in range(0, 360):
        glon2, glat2, baz = shoot(glon1, glat1, azimuth, radius)
        X.append(glon2)
        Y.append(glat2)
    X.append(X[0])
    Y.append(Y[0])
 
    X,Y = m(X,Y)
    ax.plot(X,Y,**kwargs)

def plot_circle_on_map(m, centerlon, centerlat, radius, **kwargs):
    """
    Wrapper for equi
    
    m : mpl_toolkits.Basemap
    centerlon : float
    centerlat : float
    radius : float
    kwargs : axis kwargs
    return : equi
    """
    return equi(m, centerlon, centerlat, radius, **kwargs)

def plot_line_on_map(m, point_1, point_2, s=5, l=5, color='red'):
    """
    plots line on matplotlib basemap map

    m : mpl_toolkits.Basemap
    point_1 : list
    point_2 : list
    s : float
    l : float
    color : str
    return : None
    """
    x = (point_1[0], point_2[0])
    y = (point_1[1], point_2[1])
    x, y = m(x, y)
    m.plot(x, y, marker='D', color=color, markersize=s, linewidth=l)
    
def plot_text_on_map(m, ax, lat, lon, text, fontsize=15):
    """
    plots text on matplotlib basemap map

    m : mpl_toolkits.Basemap
    ax : mpl figure axes
    lat : float
    lon : float
    text : str
    fontsize : int
    """
    x, y = m(lon, lat)
    ax.text(s=text, x=x, y=y, fontsize=fontsize)

def draw_screen_poly(lats, lons, m, ax):
    """
    draws polygon on matplotlib basemap map

    lats : list
    lons : list
    m : mpl_toolkits.Basemap
    ax : mpl figure axes
    return : None
    """
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = pg( xy, facecolor='None', edgecolor='red', linestyle='--', linewidth=5 )
    ax.add_patch(poly)
    
def plot_seismicity_map(dataframe, lon_lat_min_max=None, **kwargs):
    """
    Plots seismicity map for provided catalog
    
    Plots depth as color, magnitude proportional
    to size. Plots down to 100km. >100km depths
    are reduced to 100km in color chart.
    
    Fails on catalogs with >1 million events
    (perhaps smaller).
    
    dataframe : pandas.DataFrame
    lon_lat_min_max : list
    kwargs : figure axes kwargs
    """
    df = dataframe.copy()
    if lon_lat_min_max is None:
        lat_min = np.floor(df.lat.min())
        lat_max = np.ceil(df.lat.max())
        lon_min = np.floor(df.lon.min())
        lon_max = np.ceil(df.lon.max())
    else:
        lon_min, lon_max, lat_min, lat_max = lon_lat_min_max

    fig, ax = plt.subplots(1, figsize=(8,8))
    m = Basemap(projection='merc'
           ,llcrnrlat=lat_min
           ,urcrnrlat=lat_max
           ,llcrnrlon=lon_min
           ,urcrnrlon=lon_max
           ,resolution='i'
           ,area_thresh=1000
           ,ax=ax)

    lat_labels = np.arange(lat_min, lat_max, int((lat_max - lat_min) / 4) + 1)
    lon_labels = np.arange(lon_min, lon_max, int((lat_max - lat_min) / 4) + 1)

    m.drawparallels(lat_labels, labels=lat_labels)
    m.drawmeridians(lon_labels, labels=lon_labels)
    
    m.drawcoastlines()
    m.fillcontinents(color='0.72', zorder=0)
    
    x, y = m(df.lon.values, df.lat.values)
    # TODO : make the color user changeable
    cbar = ax.scatter(x, y, c=df.depth.values, s=1*np.exp(df.mag.values/2.), edgecolor='None'
                  , cmap='rainbow', alpha=0.5, vmin=0, vmax=100)
    c1 = fig.colorbar(cbar, label='depth (km)',fraction=0.0346, pad=0.084)
    c1.ax.invert_yaxis()
    
    return m, fig, ax