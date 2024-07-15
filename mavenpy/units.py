from collections.abc import Iterable

import numpy as np

from .helper import broadcast_index
from .constants import hc_eVum, hc_Jm, eV_to_joules


# Unit conversion functions:


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
