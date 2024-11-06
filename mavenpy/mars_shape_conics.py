from . import helper

import numpy as np
from scipy.interpolate import CubicSpline

# Routines for calculating commonly used conic section
# fits for the bow shock and the magnetic pileup boundary,
# and determining when a given (x, y, z) MSO coordinate
# is in the solar wind / MPB / optical shadow.

# Conic section fits based on maven_orbit_tplot.pro.

bow_shock_parameters =\
    {"Trotignonetal2006": {"x0": 0.6, "e": 1.026, "L": 2.081},
     "Vignesetal2000_Slavin": {"x0": 0.72, "e": 1.02, "L": 1.93},
     "Vignesetal2000_directfit": {"x0": 0.64, "e": 1.03, "L": 2.04}}

MPB_parameters =\
    {"Trotignonetal2006_twoconics":
        {"x0": (0.64, 1.60), "e": (0.77, 1.009), "L": (1.08, 0.528),
         "x": (1, -1)}}

# Colors used in IDL for each region:
region_colors = {"sw": "k", "pileup": "orange", "shadow": "b", "sheath": "g"}

# Mars radius used for identifying SW / sheath regions:
Mars_radius = 3389.5


def conic_section(x_offset, L, eccentricity, phi):
    '''Returns Cartesian coordinates (x, y, z) in Rm
    for a conic section given the x coordinate of the
    focus (x_offset), semilatus rectum (L), eccentricity,
    and angle along the x-axis (phi = 0 at the subsolar point).

    Follows the equation used in Trotignon et al. [2006]
    and Slavin and Holzer [1981]:

    y^2 + z^2 = (e^2 - 1) (x - x_f)^2 - 2eL(x - x_f) + L^2

    x_f: focus position along the x-axis
    L: semi-latus rectum
    e: eccentricity

    This can be parameterized in polar coordinates s.t.
    r cos(theta) = x - x_f and r sin(theta) = y = z,
    giving r (1 + e cos(theta)) = L.'''

    r = L / (1.0 + eccentricity * np.cos(phi))
    x = x_offset + r * np.cos(phi)
    y = r * np.sin(phi)
    z = r * np.sin(phi)

    return x, y, z


def bow_shock(reference='Trotignonetal2006'):

    '''Returns x, y, z coordinates of the bow shock
    in Mars radii.'''

    param = bow_shock_parameters[reference]
    phi = np.radians(np.arange(-150, 151, 1))

    return conic_section(param["x0"], param["L"], param["e"], phi)


def MPB(reference="Trotignonetal2006_twoconics"):

    param = MPB_parameters[reference]

    if "twoconics" in reference:

        # Get conic sections for both
        phi = np.radians(np.linspace(-160, 160, 321))
        x1, y1, z1 = conic_section(
            param["x0"][0], param["L"][0], param["e"][0], phi)
        x2, y2, z2 = conic_section(
            param["x0"][1], param["L"][1], param["e"][1], phi)

        p = np.where(x1 >= 0)[0]
        n = np.where(x2 < 0)[0]

        # Select only x > 0 from the first conic section
        # and x < 0 for the second and merge:
        x = np.concatenate((x2[n], x1[p]))
        y = np.concatenate((y2[n], y1[p]))
        z = np.concatenate((z2[n], z1[p]))

        # Sort all arrays s.t. z is increasing:
        sorted_index = np.argsort(z)
        x = x[sorted_index]
        y = y[sorted_index]
        z = z[sorted_index]

    return x, y, z


# Routines that return indices for a given
# (x, y, z) where in the optical shadow,
# solar wind, or pileup.


def km_to_Rm(x_km, y_km, z_km, Rm_km=Mars_radius):

    x_Rm = x_km / Rm_km
    y_Rm = y_km / Rm_km
    z_Rm = z_km / Rm_km

    return x_Rm, y_Rm, z_Rm


def shadow_indices(x, y, z, Rm=Mars_radius):
    x, y, z = km_to_Rm(x, y, z, Rm_km=Rm)
    s = np.sqrt(y ** 2 + z ** 2)
    shadow_index = np.where((s < 1) & (x < 0))[0]
    return shadow_index


def solar_wind_indices(x, y, z, reference='Trotignonetal2006',
                       Rm=Mars_radius):

    x, y, z = km_to_Rm(x, y, z, Rm_km=Rm)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    altitude = (r - 1.0) * Mars_radius
    s = np.sqrt(y ** 2 + z ** 2)

    param = bow_shock_parameters[reference]
    x0 = param["x0"]
    ecc = param["e"]
    L = param["L"]
    phm = np.radians(160)

    phi = np.arctan2(s, (x - x0))
    rho_s = np.sqrt((x - x0) ** 2 + s ** 2)
    # shock = L/(1. + ecc*np.cos(phi))
    shock = np.where(
        phi < phm, L / (1.0 + ecc * np.cos(phi)),
        L / (1.0 + ecc * np.cos(phm)))
    sw_index = np.where(rho_s >= shock)[0]
    not_sw_index = np.where(rho_s < shock)[0]

    return sw_index, altitude, not_sw_index


def pileup_indices(x, y, z, reference="Trotignonetal2006_twoconics",
                   Rm=Mars_radius):

    x, y, z = km_to_Rm(x, y, z, Rm_km=Rm)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    altitude = (r - 1.0) * Mars_radius
    s = np.sqrt(y ** 2 + z ** 2)

    rho_p = np.ones(x.size)
    MPB = np.ones(x.size)

    param = MPB_parameters[reference]
    x0_p, x0_n = param["x0"]
    ecc_p, ecc_n = param["e"]
    L_p, L_n = param["L"]

    indx = np.where(x >= 0)[0]
    phi = np.arctan2(s, (x - x0_p))
    rho_p[indx] = np.sqrt((x[indx] - x0_p) ** 2 + s[indx] ** 2)
    MPB[indx] = L_p / (1.0 + ecc_p * np.cos(phi[indx]))

    # plt.figure()
    # plt.plot(rho_p)
    # plt.plot(MPB)

    indx = np.where(x < 0)[0]
    phi = np.arctan2(s, (x - x0_n))
    phm = np.radians(160)

    rho_p[indx] = np.sqrt((x[indx] - x0_n) ** 2 + s[indx] ** 2)
    MPB[indx] = np.where(
        phi[indx] < phm,
        L_n / (1.0 + ecc_n * np.cos(phi[indx])),
        L_n / (1.0 + ecc_n * np.cos(phm))
    )

    # plt.figure()
    # plt.plot(rho_p)
    # plt.plot(MPB)
    # plt.show()

    sheath_index = np.where(rho_p >= MPB)[0]
    pileup_index = np.where(rho_p < MPB)[0]

    return sheath_index, pileup_index, altitude


def region_index(x, y, z,
                 region_names=('sw', 'sheath', 'pileup', 'shadow')):

    sw = solar_wind_indices(x, y, z)[0]
    shadow = shadow_indices(x, y, z)
    sheath, pileup = pileup_indices(x, y, z)[:2]

    sheath = np.array([i for i in sheath if i not in sw])
    pileup = np.array([i for i in pileup if i not in shadow])

    index_by_region = {}
    if "sw" in region_names:
        index_by_region["sw"] = sw
    if "sheath" in region_names:
        index_by_region["sheath"] = sheath
    if "pileup" in region_names:
        index_by_region["pileup"] = pileup
    if "shadow" in region_names:
        index_by_region["shadow"] = shadow

    return index_by_region


def region_separation(x, y, z, arr):
    '''Break array 'arr' (e.g. spacecraft altitude)
    into solar wind / pileup / shadow
    / sheath components based on spacecraft x, y, z.'''

    region_names = ("sw", "sheath", "pileup", "shadow")

    index_by_region = region_index(x, y, z, region_names=region_names)

    arr_shape = arr.shape
    arr_dim = len(arr_shape)

    region_arr = {}

    for name_i in region_names:
        arr_i = np.zeros(shape=arr_shape) + np.nan

        index_i = index_by_region[name_i]

        if arr_dim > 1:
            broadcast_index_i = helper.broadcast_index(
                arr, x, matching_axis_index=index_i,
                other_axis_index=Ellipsis)
        else:
            broadcast_index_i = index_i
        arr_i[broadcast_index_i] = arr[broadcast_index_i]
        region_arr[name_i] = arr_i

    return region_arr


def cartesian_spline(t, x, y, z):
    x_spline = CubicSpline(t, x)
    y_spline = CubicSpline(t, y)
    z_spline = CubicSpline(t, z)

    return x_spline, y_spline, z_spline
