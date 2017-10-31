# shapely polygon selection

from shapely.geometry import Point
import numpy as np
from scipy import spatial

def cartesian_distance_between_two_three_vectors(vector_a, vector_b):
    """
    Calculates cartesian distance between two 3d points
    """
    dist = spatial.distance.cdist([vector_a], [vector_b])
    return dist

def distance_between_two_coordinates(lat_a, lon_a, lat_b, lon_b):
    """
    calculates the distance between two LAT/LON coordinates
    
    accurate up to 200km distances
    
    returns single LON/LAT vector difference
    """
    lat_diff = lat_a - lat_b
    lon_diff = lon_a - lon_b
    lat_km = 111.19 * lat_diff
    
    lat_radians = np.deg2rad(lat_a)
    lon_km = 111.19 * lon_diff * np.cos(lat_radians)
    distance = np.array([lon_km, lat_km])
    x, y = distance
    return np.sqrt(x**2 + y**2)

def get_node_data(node, radius, data, m):
    """
    returns data within a circle with given radius
    """
    node_lon = node[0]
    node_lat = node[1]

    distance_from_node = (radius * 1.2) / 111.19

    lon_bounds = [node_lon - distance_from_node, node_lon + distance_from_node]
    lat_bounds = [node_lat - distance_from_node, node_lat + distance_from_node]

    df = data.copy()
    df = df[df.lon.between(lon_bounds[0], lon_bounds[1]) &
            df.lat.between(lat_bounds[0], lat_bounds[1])]
    if df.shape[0] > 0:
        df['distance'] = df.apply(lambda row: distance_between_two_coordinates(row['lat'], row['lon'], node_lat, node_lon),
                              axis=1)
        df = df[df.distance <= radius]
        return df
    else:
        df['distance'] = radius + 1
        df = df[df.distance <= radius]
        return df
