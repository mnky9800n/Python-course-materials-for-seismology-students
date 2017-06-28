
import numpy as np

#: Earth radius in km.
EARTH_RADIUS = 6371.0

#: Maximum elevation on Earth in km.
EARTH_ELEVATION = -8.848


def geodetic_distance(lons1, lats1, lons2, lats2, diameter=2 * EARTH_RADIUS):
    """
    Calculate the geodetic distance between two points or two collections
    of points.
    Parameters are coordinates in decimal degrees. They could be scalar
    float numbers or numpy arrays, in which case they should "broadcast
    together".
    Implements http://williams.best.vwh.net/avform.htm#Dist
    :returns:
        Distance in km, floating point scalar or numpy array of such.
    """
    lons1, lats1, lons2, lats2 = _prepare_coords(lons1, lats1, lons2, lats2)
    distance = np.arcsin(np.sqrt(
        np.sin((lats1 - lats2) / 2.0) ** 2.0
        + np.cos(lats1) * np.cos(lats2)
        * np.sin((lons1 - lons2) / 2.0) ** 2.0
    ))
    return diameter * distance


def azimuth(lons1, lats1, lons2, lats2):
    """
    Calculate the azimuth between two points or two collections of points.
    Parameters are the same as for :func:`geodetic_distance`.
    Implements an "alternative formula" from
    http://williams.best.vwh.net/avform.htm#Crs
    :returns:
        Azimuth as an angle between direction to north from first point and
        direction to the second point measured clockwise in decimal degrees.
    """
    lons1, lats1, lons2, lats2 = _prepare_coords(lons1, lats1, lons2, lats2)
    cos_lat2 = np.cos(lats2)
    true_course = np.degrees(np.arctan2(
        np.sin(lons1 - lons2) * cos_lat2,
        np.cos(lats1) * np.sin(lats2)
        - np.sin(lats1) * cos_lat2 * np.cos(lons1 - lons2)
    ))
    return (360 - true_course) % 360


def distance(lons1, lats1, depths1, lons2, lats2, depths2):
    """
    Calculate a distance between two points (or collections of points)
    considering points' depth.
    Calls :func:`geodetic_distance`, finds the "vertical" distance between
    points by subtracting one depth from another and combine both using
    Pythagoras theorem.
    :returns:
        Distance in km, a square root of sum of squares of :func:`geodetic
        <geodetic_distance>` distance and vertical distance, which is just
        a difference between depths.
    """
    hdist = geodetic_distance(lons1, lats1, lons2, lats2)
    vdist = depths1 - depths2
    return np.sqrt(hdist ** 2 + vdist ** 2)


def _prepare_coords(lons1, lats1, lons2, lats2):
    """
    Convert two pairs of spherical coordinates in decimal degrees
    to numpy arrays of radians. Makes sure that respective coordinates
    in pairs have the same shape.
    """
    lons1 = np.radians(lons1)
    lats1 = np.radians(lats1)
    assert lons1.shape == lats1.shape
    lons2 = np.radians(lons2)
    lats2 = np.radians(lats2)
    assert lons2.shape == lats2.shape
    return lons1, lats1, lons2, lats2


def great_circle(alon, alat, a_xoffset, a_yoffset, adepth
                 , blon, blat, b_xoffset, b_yoffset, bdepth):
    """
    calculates great circle including offsets
    """
    alat_offset = a_yoffset / 111199. + alat
    alon_offset = a_xoffset / 111199. * np.cos(np.deg2rad(alat_offset)) + alon

    blat_offset = b_yoffset / 111199. + blat
    blon_offset = b_xoffset / 111199. * np.cos(np.deg2rad(blat_offset)) + blon

    return distance(alon_offset, alat_offset, adepth
                    , blon_offset, blat_offset, bdepth) * 1e3


def spherical_to_cartesian(lons, lats, depths):
    """
    Return the position vectors (in Cartesian coordinates) of list of spherical
    coordinates.

    For equations see: http://mathworld.wolfram.com/SphericalCoordinates.html.

    Parameters are components of spherical coordinates in a form of scalars,
    lists or numpy arrays. ``depths`` can be ``None`` in which case it's
    considered zero for all points.

    :returns:
        ``numpy.array`` of 3d vectors representing points' coordinates in
        Cartesian space. The array has the same shape as parameter arrays.
        In particular it means that if ``lons`` and ``lats`` are scalars,
        the result is a single 3d vector. Vector of length ``1`` represents
        distance of 1 km.

    See also :func:`cartesian_to_spherical`.
    """
    phi = np.radians(lons)
    theta = np.radians(lats)
    if depths is None:
        rr = EARTH_RADIUS
    else:
        rr = EARTH_RADIUS - np.array(depths)
    cos_theta_r = rr * np.cos(theta)
    xx = cos_theta_r * np.cos(phi)
    yy = cos_theta_r * np.sin(phi)
    zz = rr * np.sin(theta)
    vectors = np.array([xx.transpose(), yy.transpose(), zz.transpose()]) \
        .transpose()
    return vectors


def cartesian_distance(alon, alat, a_xoffset, a_yoffset, adepth
                       , blon, blat, b_xoffset, b_yoffset, bdepth):
    """
    Calculates the cartesian distance between two points including offsets.
    """
    a = spherical_to_cartesian(lons=alon, lats=alat, depths=adepth)
    b = spherical_to_cartesian(lons=blon, lats=blat, depths=bdepth)

    a[0] = a[0] + a_xoffset * 1e3
    a[1] = a[1] + a_yoffset * 1e3

    b[0] = b[0] + b_xoffset * 1e3
    b[1] = b[1] + b_yoffset * 1e3

    dist = np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) * 1e3
    return dist

def gc_dist(a, b):
    """
    Wrapper for great_circle.
    :param a: 
    :param b: 
    :return: 
    """
    alon, alat, a_xoffset, a_yoffset, adepth = a
    blon, blat, b_xoffset, b_yoffset, bdepth = b
    return great_circle(alon, alat, a_xoffset, a_yoffset, adepth
                 , blon, blat, b_xoffset, b_yoffset, bdepth)

def c_dist(a, b):
    """
    Wrapper for cartesian_distance.

    :param a: vector
    :param b: 
    :return: 
    """
    alon, alat, a_xoffset, a_yoffset, adepth = a
    blon, blat, b_xoffset, b_yoffset, bdepth = b
    return cartesian_distance(alon, alat, a_xoffset, a_yoffset, adepth
                 , blon, blat, b_xoffset, b_yoffset, bdepth)
