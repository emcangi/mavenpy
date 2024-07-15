from collections.abc import Iterable

import numpy as np
import spiceypy

from . import coordinates
from .spice import dt_to_et


# Adaptation of mvn_sep_anc_fov_mars_fraction.pro,
# which is used to calculate the SEP-??_FRAC_FOV_MARS
# variable saved in the ancillary SEP data files.


SEP_fov_NAIF_ID =\
    {"SEP1A_front": -202126, "SEP1B_front": -202127,
     "SEP1A_back": -202128, "SEP1B_back": -202129,
     "SEP2A_front": -202121, "SEP2B_front": -202122,
     "SEP2A_back": -202123, "SEP2B_back": -202124}

# No difference between 1A_front and 1B_front

SEP_fovs = {"1F": "SEP1A_front", "1R": "SEP1A_back",
            "2F": "SEP2A_front", "2R": "SEP2A_back"}


def get_fov_pixels(sep_sensor, sep_look_direction, pixel_width_deg=1.5):

    '''Returns a list of (x, y, z) coordinates that describe the boresight
    for a given SEP sensor-look direction.
    '''

    sep_s_ld = "{}{}".format(sep_sensor, sep_look_direction)

    naif_code = SEP_fov_NAIF_ID[SEP_fovs[sep_s_ld]]

    # Ask Spice getfov for the corner vectors of the FOV
    # Returns: shape (e.g. rectangle), frame, boresight, # edges
    #     bounds (array of 4 3-D vectors for rectangle)
    shape, frame_i, boresight, n_edges, corners = spiceypy.getfov(
        naif_code, 4)

    # Make a range of pixels around the boresight stretching to the corners
    # r is 1
    r, theta_bounds, phi_bounds = coordinates.cartesian_to_spherical(
        corners[:, 0], corners[:, 1], corners[:, 2])
    # print(np.degrees(theta_bounds))
    # print(np.degrees(phi_bounds))

    # Get pixel width in radians:
    pixel_width_rad = np.radians(pixel_width_deg)

    # Discretize the cone angle into pixels with width of ~1.5 deg apart:
    # # Method from SEP IDL routine, not needed but kept for documentation:
    # theta_range = max(theta_bounds) - min(theta_bounds)
    # ntheta = np.ceil(theta_range/pixel_width_rad)
    # dtheta = theta_range/ntheta
    # theta_edges_array = theta_bounds[0] + dtheta*np.arange(ntheta + 1)
    # theta_centers_array = theta_edges_array[:-1] +\
    #     np.ediff1d(theta_edges_array)/2
    theta_centers = np.arange(
        min(theta_bounds), max(theta_bounds), pixel_width_rad)
    # print(np.degrees(theta_centers))

    # Do the same for the phi angle:
    if phi_bounds[3] > phi_bounds[0]:
        phi_centers = np.arange(
            min(phi_bounds), max(phi_bounds), pixel_width_rad)
    else:
        phi_centers = np.arange(
            max(phi_bounds), min(phi_bounds) + 2*np.pi, pixel_width_rad) % 360
    # print(np.degrees(phi_centers))

    theta_centers_2d, phi_centers_2d = np.meshgrid(theta_centers, phi_centers)

    # fov_x, fov_y, fov_z = spherical_to_cartesian(
    #     1, theta_centers_array[np.newaxis, :],
    #     phi_centers_array[:, np.newaxis])
    fov_x, fov_y, fov_z = coordinates.spherical_to_cartesian(
        1, theta_centers_2d, phi_centers_2d)

    return fov_x, fov_y, fov_z, theta_centers_2d, phi_centers_2d


def fraction_Mars_in_FOV(sensor, look_direction, time_utc,
                         pixel_width_deg=1.5):
    '''Returns the fraction of the field of view
    of a given SEP sensor is taken up by (sunlit) Mars

    time_utc: str or datetime or list of datetimes
    '''

    # Convert to ephemeris time:
    ephemeris_time = dt_to_et(time_utc)

    # Spice frame that all SEP boresight vectors are in:
    frame = "MAVEN_SEP{}".format(sensor)

    # Get the pixel ranges that will be iterated across
    # to see if Mars intersects:
    fov_x, fov_y, fov_z, fov_theta_deg, fov_phi_deg = get_fov_pixels(
        sensor, look_direction, pixel_width_deg=pixel_width_deg)

    # Flatten the arrays for iterating:
    fov_x = fov_x.flatten()
    fov_y = fov_y.flatten()
    fov_z = fov_z.flatten()

    if not isinstance(ephemeris_time, Iterable):
        ephemeris_time = [ephemeris_time]

    n_time = len(ephemeris_time)
    n_fov_vec = len(fov_x)

    mars_in_fov = np.zeros(shape=(n_time, n_fov_vec))
    sunlit_mars_in_fov = np.zeros(shape=(n_time, n_fov_vec))

    # Disables the error when no intercept:
    with spiceypy.no_found_check():
        for j, et_j in enumerate(ephemeris_time):
            for i, fov_xyz_i in enumerate(zip(fov_x, fov_y, fov_z)):
                # print(i, fov_xyz_i)
                # Just pxform takes 6 seconds/100 time segments:
                # sincpt (Surface intercept):
                # Takes 6 seconds/100 time segments
                spoint, trgepc, srfvec, intercept = spiceypy.sincpt(
                    'Ellipsoid', 'Mars', et_j,
                    'IAU_MARS', 'NONE', 'MAVEN', frame,
                    fov_xyz_i)

                mars_in_fov[j, i] = intercept

                if intercept:
                    # Adds 0.5 sec/100 time segments

                    trgepc, srfvec, phase_angle,\
                        solar_zenith_angle, emission_angle = spiceypy.ilumin(
                            'Ellipsoid', 'MARS', et_j, 'IAU_MARS',
                            'NONE', 'MAVEN', spoint)

                    sunlit_mars_in_fov[j, i] =\
                        (solar_zenith_angle < np.pi/2) *\
                        np.cos(solar_zenith_angle)

    weight = np.sin(fov_theta_deg.flatten())[np.newaxis, :]

    frac_fov = np.sum(mars_in_fov*weight, axis=1)/np.sum(weight)
    frac_illum_fov = np.sum(sunlit_mars_in_fov*weight, axis=1)/np.sum(weight)

    return frac_fov, frac_illum_fov
