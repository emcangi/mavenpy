
import os

from .read import read_cdf, read_sav
from .helper import invert_energy_axis, process_data_dict


swea_units =\
    {"time_unix": "s, POSIX seconds after 1970", "epoch": "utc",
     "num_accum": "# of 2-sec accumulations per observation",
     "energy": "eV",
     "de_over_e": "unitless, FWHM energy resolution",
     "counts": "#",
     "diff_en_fluxes": "eV/(cm2 s sr eV)",
     "binning": "unitless, energy array index stride",
     "elev": "deg, instrument elevation angle (energy-varying)",
     "azim": "deg, instrument azimuthal angle",
     "pa": "deg., 64 x 16 array of center pitch angles (time-varying)",
     "d_pa": "deg., pitch angle range spanned by each bin"
     }


#########################################
#           SWEA routines               #
#########################################


def read(data_file, dataset_type="None",
         fields=None, spectra_type="eflux", lib='cdflib',
         include_unit=True):

    """Returns SWEA data.

    This routine is based on the SWEA SIS document and the MAVEN SPEDAS
    routine 'mvn_swe_load_l2.pro'. We thank and credit the
    authors of both.

    dataset_type: name of dataset to be pulled. options:
        - pad: electron counts/differential energy flux binned
            into 64-energy x 16-azimuth angle x 6-elevation angle
            array with a binning factor indicating whether adjacent
            energy-varying elements were duplicated (binning = 2)
            or quadriplicated (binning = 4)
        - pad: electron pitch angle distribution with dimensions
            64 (energies) x 16 (pitch angles). pitch angles
            vary as a function of time and energy bin,
            so has same dim as the flux
        - spec: angle-averaged electron energy spectra
           for 64 energies, effective accumul time ~0.4 sec
           - Note: while "arc" versions of 3d and pad exist,
           an "arc" version of spec doesnt exist since telemetry
           bandwidth ended up higher than planned/
        - scpot: 

    data_file: name of file to be read
    """

    if dataset_type == "None":
        filename = os.path.split(data_file)[1]
        if "spec" in filename:
            dataset_type = "s"
        elif "3d" in filename:
            dataset_type = "3"
        elif "pad" in filename:
            dataset_type = "p"

    if not fields:
        if dataset_type == "s":
            # spectra
            fields = ("time_unix", "epoch",
                      "diff_en_fluxes", "counts",
                      "energy", "de_over_e",
                      "num_accum")

        if dataset_type == "p":
            # PAD
            fields = ("time_unix", "epoch",
                      "diff_en_fluxes", "counts",
                      "energy", "de_over_e",
                      "pa", "d_pa",
                      "binning")

        if dataset_type == "3":
            # 3D
            fields = ("time_unix", "epoch",
                      "diff_en_fluxes", "counts",
                      "energy", "elev", "azim",
                      "binning")

    # Pull data from the SWIA Level 2 CDF
    data_cdf = read_cdf(data_file, fields, lib=lib)

    data_cdf = invert_energy_axis(
        data_cdf, energy_names=('energy',),
        energy_dependent_var_axis=-1)

    if include_unit:
        data_cdf = process_data_dict(data_cdf, units=swea_units)

    return data_cdf


#########################################
#      Spacecraft potential             #
#########################################


scpot_method =\
    {-1: "Invalid", 0: "Manual", 1: "pot_swelpw", 2: "pot_swepos",
     3: "pot_sweneg", 4: "pot_sta", 5: "pot_sweshdw"}


def read_scpot(filename, source=('comp',)):
    '''Retrieve SWEA Level-3 spacecraft potential files in the original IDL sav
    format. Note this is very slow.

    source: defaults to 'comp' (for the merged SC potential).'''
    if source == "all":
        source = ("comp", "swelpw", "swepos", "sweneg", "sta", "sweshdw")

    scpot = read_sav(filename)

    scpot_struct = {}
    for name in scpot:
        name_split = name.split("_")
        if name_split[1] not in source:
            continue

        # Skip the unit name (V) and methods category
        if name == "pot_comp_units_name" or name == "pot_comp_pot_name":
            continue

        if "time" in name:
            unit = "unix"
        elif "potential" in name:
            unit = "V"
        elif name == "pot_comp_method":
            unit = [i.decode() for i in scpot["pot_comp_pot_name"][0]]
        else:
            raise IOError()

        scpot_struct[name] = (scpot[name], unit)

    return scpot_struct
