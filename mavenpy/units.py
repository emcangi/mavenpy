from collections.abc import Iterable

import numpy as np

from .helper import broadcast_index
from .constants import hc_eVum, hc_Jm, eV_to_joules


# Unit conversion functions:

ev2joules = 1.60217646e-19  # Joules in 1 eV
Avagadro = 6.022e23  # #/mol
speed_of_light_kms = 299792.458  # km/s
boltzmann_constant_JK = 1.38e-23  # J/K
boltzmann_constant_eVK = 8.617e-5  # eV/K
mu_0 = 1.256e-6  # Vacuum permissivity

h_Js = 6.626e-34  # J s
h_eVs = 4.135e-15  # eV s
hbar_Js = 1.054e-34  # J * s
hbar_eVs = 6.582e-16  # eV * s

# Subatomic masses:
electron_mass = 9.10938215e-31  # kg
proton_mass = 1.67262178e-27  # kg

# Molecular mass in g/mol:
mass_amu_dict =\
    {"H": 1.00794, "H2": 2.01568, "He": 4.002602, "CO2": 44.01,
     "CO": 28.010, "N2": 28.02, "O1": 15.9994, "O2": 31.9988, "Ar": 39.948,
     "H2O": 18.01528, "CH4": 16.04}

# Get the molecular mass in kg
mass_kg_dict = {n: i/Avagadro*1e-3 for n, i in mass_amu_dict.items()}
mass_kg_dict["e-"] = electron_mass
mass_kg_dict["H+"] = proton_mass

particle_num_charge = {"H+": 1, "H0": 0, "e-": -1, "e+": 1, "H-": -1}


def mass(mass_kg=None, name=None):
    '''Function to return mass given a mass or name,
    checks if name already in list.'''

    if not mass_kg and not name:
        raise ValueError("Must provide mass (kg) or name (str, e.g. H+).")
    elif mass_kg:
        return mass_kg
    elif name not in mass_kg_dict:

        if 'e' in name:
            return mass(name='e-')
        elif 'p' in name:
            return mass(name='H+')

        raise ValueError(
            "Mass not known for name '{}', supply mass"
            " in kilograms instead.".format(name))
    else:
        return mass_kg_dict[name]


def lorentz_factor_from_velocity(v_kms):
    '''Returns relativstic Lorentz factor (gamma) given a
    scalar velocity in km/s'''

    # g = 1 / sqrt(1 - v^2/c^2)
    return 1 / np.sqrt(1 - (v_kms / speed_of_light_kms) ** 2)


def rest_mass(mass_kg, unit='J'):
    '''Returns rest mass in units of Joules
    '''

    # M0 = Mc^2 (kgm2s2 = J)
    M0 = mass_kg * (speed_of_light_kms * 1e3) ** 2

    if unit == 'eV':
        M0 = M0 / ev2joules

    return M0


def energy(v_kms, mass_kg=None, name=None):
    """ Calculate (kinetic) energy of a particle given mass and velocity.
    Accounts for relativistic effects.

    v_kms: float, (km/s)
    mass: float, (kg), optional
    name: string, particle name to lookup mass.

    Returns energy, (eV) """

    # E = kg * (km/s)**2 * (1e3)**2 / (J/eV)
    #   = kg * (km*1e3/s)**2 / (J/eV)
    #   = kg * (m/s)**2 / (J/eV)
    #   = eV
    # Get mass:
    mass_kg = mass(mass_kg=mass_kg, name=name)

    # Calculate lorentz factor:
    g = lorentz_factor_from_velocity(v_kms)

    # Rest mass:
    rest_mass_eV = rest_mass(mass_kg, unit='eV')

    # Non relativstic commented out:
    # return 0.5 * mass_kg * (v_kms * 1e3)**2 / ev2joules

    return rest_mass_eV * (g - 1)


def scalar_velocity(energy_eV, mass_kg=None, name=None):
    '''Returns a scalar velocity for an energy and species or mass'''

    # Get the mass:
    mass_kg = mass(mass_kg=mass_kg, name=name)

    # Non-relativistic commented out:
    # v_kms = 1e-3 * math.sqrt(2. * energy_eV * ev2joules / mass)

    # Rest mass:
    rest_mass_eV = rest_mass(mass_kg, unit='eV')

    # Get ratio of energy to rest mass:
    ratio = energy_eV / rest_mass_eV

    v_kms = speed_of_light_kms * np.sqrt(
        1 - (1 / (1 + ratio)) ** 2)

    return v_kms


# Pressures (for plasma pressure balance):

def dynamic_pressure(density_cm3, velocity_kms, mass_kg=None, name=None):
    '''Returns in nPa the dynamic plasma pressure
    (`r u^2 = m_i n_i u_i^2)'''

    # Get the mass:
    mass_kg = mass(mass_kg=mass_kg, name=name)

    # P_dyn = m_H * n * u^2
    #   = kg cm-3 * (km/s)^2 = kg m-3 * (10^2)^3 * (10^3 m/s)^2
    #   = # kg / s^2 / m *1e12 = Pa
    P_dyn = density_cm3 * velocity_kms**2 * mass_kg * 1e21

    return P_dyn


def thermal_pressure(density_cm3, temperature_eV):
    '''Returns thermal pressure in nPa
    (P_th = n_i k T_i)'''

    # Thermal pressure
    # P_th = N k T = eV/cm3 = J/m3 (eV2joules * (100)^3)
    P_therm = (density_cm3 * 1e6) * (temperature_eV * ev2joules) * 1e9

    return P_therm


def magnetic_pressure(B_nT):
    '''Returns magnetic pressure in nPa
    (P_mag = B^2 / 2 mu_0)'''

    # MAgnetic pressure
    # N/A^2 * (nT)^2 = 1e-18 N T^2 /A^2

    P_mag = (B_nT * 1e-9)**2/(2 * mu_0) * 1e9

    return P_mag


def B_angles(bx, by, bz):

    # Calculate clock angle and cone angle
    clock_angle = np.degrees(np.arctan2(by, bz)) % 360
    cone_angle = np.degrees(np.arctan2(by, bx)) % 360

    return clock_angle, cone_angle


def wavelength_to_energy(wavelength_nm, energy_unit='eV'):
    '''Calculates the energy (eV) of a photon
    of wavelength l (nm) based on E = hc/l
    where hc is the planck constant.'''

    if isinstance(wavelength_nm, Iterable) and not\
            isinstance(wavelength_nm, np.ndarray):
        wavelength_nm = np.array(wavelength_nm)
    wavelength_nm = wavelength_nm.flatten()

    if energy_unit == "eV":
        energy = hc_eVum * 1e3 / wavelength_nm
    elif energy_unit == "J":
        energy = hc_Jm * 1e9 / wavelength_nm

    return energy


def change_unit(spectra, wavelength_nm=None, energy_eV=None,
                initial_unit="W/m2/nm", target_unit="#/cm2/s/nm"):

    # Convert a spectra of arbitrary dimensionality from a given unit
    # (irradiance, W per m2 per nm; or energy flux, eV/s per cm2 per nm)
    # to photon flux (N photons / s / cm^2 / nm).

    # Get the numerator and divisor for each unit:
    initial_unit = initial_unit.replace("^", "")
    initial_unit = initial_unit.split("/")
    initial_numerator = initial_unit[0]
    initial_divisors = initial_unit[1:]

    target_unit = target_unit.replace("^", "")
    target_unit = target_unit.split("/")
    target_numerator = target_unit[0]
    target_divisors = target_unit[1:]

    # Check if energy or wavelength needed:
    conditions =\
        (("#" in initial_numerator and "eV" in target_numerator),
         ("eV" in initial_numerator and "#" in target_numerator),
         ("eV" in initial_divisors and "nm" in target_divisors),
         ("nm" in initial_divisors and "eV" in target_divisors))

    if sum(conditions) > 0:
        if energy_eV is None and wavelength_nm is None:
            raise IOError(
                "Need to provide energy or "
                "wavelength to convert to # flux!")
        elif wavelength_nm is not None:
            if isinstance(wavelength_nm, Iterable) and not\
                    isinstance(wavelength_nm, np.ndarray):
                wavelength_nm = np.array(wavelength_nm)
            wavelength_nm = wavelength_nm.flatten()
            energy_eV = wavelength_to_energy(
                wavelength_nm, energy_unit='eV')

        # print(energy_eV)

        dE = np.ediff1d(energy_eV, to_end=(energy_eV[-1] - energy_eV[-2]))
        dl = np.ediff1d(
            wavelength_nm, to_end=(wavelength_nm[-1] - wavelength_nm[-2]))
        dE_dl = -dE/dl
        # dE_dl = (hc_eVum * 1e3)/wavelength_nm**2

        # from matplotlib import pyplot as plt
        # plt.plot(wavelength_nm, -dE/dl, label='-(d(hc/λ)/dλ)')
        # plt.plot(wavelength_nm,
        #          (hc_eVum * 1e3)/wavelength_nm**2, label='hc/λ^2')
        # plt.yscale('log')
        # plt.xscale('log')
        # plt.legend()
        # plt.ylabel('dE/dλ, eV/nm')
        # plt.xlabel('λ, nm')
        # plt.show()

    # If need to convert the numerator:
    scaling_factor = 1
    if initial_numerator != target_numerator:
        # For mW, need to scale down or up accordingly:
        if initial_numerator == "mW":
            scaling_factor /= 1e3
            initial_numerator = "W"
        if target_numerator == "mW":
            scaling_factor *= 1e3
            target_numerator = "W"

        # For energy unit conversion (eV/s <-> J):
        if initial_numerator == "eV" and target_numerator == "W":
            scaling_factor *= eV_to_joules
        elif initial_numerator == "W" and target_numerator == "eV":
            scaling_factor /= eV_to_joules

        # For energy flux (eV or W) <-> flux (#):
        # I = flux / energy_per_photon = F / (hc/λ)
        #   = (W/cm^2/nm) / J = n_photons/cm^2/s/nm
        #   = (eV/cm^2/nm/s) / eV = n_photons/cm^2/s/nm
        if ("#" in target_numerator) ^ ("#" in initial_numerator):
            # Check if energy or wavelength provided to calculate
            # the photon energy:
            if initial_numerator == "W" or target_numerator == "W":
                # Get the photon energy in J
                photon_energy = eV_to_joules*energy_eV
            else:
                photon_energy = energy_eV

            if "#" in target_numerator:
                # W or eV -(divide by energy per photon)-> #/s
                scaling_factor /= photon_energy
            if "#" in initial_numerator:
                # #/s -(multiply by energy per photon)-> W or eV
                scaling_factor *= photon_energy

    # Now compare the divisors, and if different, scale:
    if initial_divisors != target_divisors:
        # Area scaling:
        if "m2" in initial_divisors and "cm2" in target_divisors:
            # converts 1/m^2 -> 1/cm^2
            scaling_factor = 1e-4 * scaling_factor
        elif "cm2" in initial_divisors and "m2" in target_divisors:
            # converts 1/cm^2 -> 1/m^2
            scaling_factor = 1e4 * scaling_factor
        # Energy normalization:
        # dE/dl = d(hc/l)/dl = -hc/l^2
        # dE = (-hc / l^2) dl
        # so F ~ 1/dE(eV) ~ l^2 /hc (1/-dl(nm))
        if "eV" in initial_divisors and "nm" in target_divisors:
            scaling_factor *= dE_dl
        elif "nm" in initial_divisors and "eV" in target_divisors:
            scaling_factor *= 1/dE_dl

    if isinstance(scaling_factor, Iterable):
        broadcast_index_i = broadcast_index(
            spectra, wavelength_nm)[1]
        scaling_factor = scaling_factor[broadcast_index_i]

    scaled_spectra = spectra * scaling_factor

    # Convert a [spectral] irradiance (mW/m^2[/nm]) into an energy flux
    # (eV/s/cm^2[/nm])
    # mW/m^2[/nm] * (1 m/ 10 cm)^2 / (1.602 x 10^-19 J /eV) =  eV/s/cm^2[/nm]

    return scaled_spectra
