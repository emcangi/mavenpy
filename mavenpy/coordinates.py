import numpy as np


def spherical_to_cartesian(radius, theta, phi):

    """ Calculates Cartesian coordinates for spherical
     coordinates.

     radius: float, 0 < r
     theta: float (rad), 0 <= theta <= pi
     phi: float (rad), 0 <= phi <= 2pi

     Returns a tuple (x, y, z)."""

    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    return x, y, z


def cartesian_to_spherical(x, y, z):

    """ Calculates spherical coordinates for Cartesian
     coordinates.

     x: float
     y: float
     z: float

     Returns a tuple (r, theta (rad), phi (rad))."""

    radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    # Returns theta between 0 and pi radians, where pi is in -z direction.
    theta = np.arccos(z / radius)
    # Returns phi between 0 and 2 pi radians, where pi is in -y direction.
    phi = (np.arctan2(y, x) + 2 * np.pi) % (2 * np.pi)
    return radius, theta, phi


def rlatlon_to_geographic(r, lat, lon, re=3389.0, rp=3389.0):

    theta = np.radians(90.0 - lat)
    e = 1.0 - (rp / re)
    r0 = re * (1 - e * (np.cos(theta)) ** 2)
    altitude = r - r0

    return altitude, lat, lon


def cartesian_to_geographic(x, y, z, re=3389.0, rp=3389.0):

    # Calculate altitudes assuming Mars is spherical.
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    # Spherical coordinates of collision locations
    theta = np.arccos(z / r)
    phi = (np.arctan2(y, x) + 2 * np.pi) % (2 * np.pi)

    # Which can be turned into lat / lon coordinates
    lat = np.degrees(np.pi / 2 - theta)
    lon = np.degrees(phi) % 360.0

    # Radii -> altitude
    e = 1.0 - (rp / re)
    r0 = re * (1 - e * (np.cos(theta)) ** 2)
    altitude = r - r0

    return altitude, lat, lon


def cartesian_to_rlatlon(x, y, z):

    # Calculate altitudes assuming Mars is spherical.
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    # Spherical coordinates of collision locations
    theta = np.arccos(z / r)
    phi = (np.arctan2(y, x) + 2 * np.pi) % (2 * np.pi)

    # Which can be turned into lat / lon coordinates
    lat = np.degrees(np.pi / 2 - theta)
    lon = np.degrees(phi) % 360.0
    return r, lat, lon


def rlatlon_to_cartesian(radius, latitude, longitude):

    """ Calculates Cartesian coordinates for lat lon
     coordinates.

     latitude: float (deg), -90 <= lat <= 90
     longitude: float (deg), 0 <= lon <= 360
     radius: float, 0 < r

     Returns a tuple (x, y, z)."""

    theta = np.radians(90.0 - latitude)
    phi = np.radians(longitude)

    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    return x, y, z


def cartesian_to_spherical_vector(vx, vy, vz, theta, phi):
    r_unit_vector =\
        vx * np.sin(theta) * np.cos(phi) +\
        vy * np.sin(theta) * np.sin(phi) +\
        vz * np.cos(theta)

    t_unit_vector =\
        vx * np.cos(theta) * np.cos(phi) +\
        vy * np.cos(theta) * np.sin(phi) -\
        vz * np.sin(theta)

    p_unit_vector = -vx * np.sin(phi) + vy * np.cos(phi)

    return r_unit_vector, t_unit_vector, p_unit_vector
