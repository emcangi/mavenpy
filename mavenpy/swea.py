
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
     "d_pa": "deg., pitch angle range spanned by each bin",
     "quality": "unitless, quality flag that marks if affected by "
                "sporadic low energy (below 28 eV) anomaly, "
                "(0 - affected, 1 - unknown [common in sheath],"
                " 2 - not affected)",
     "secondary": "eV/(cm2 s sr eV), secondary electron contamination"
     }


edependent_vars =\
    ("diff_en_fluxes", "counts", "secondary", "variance")
e_vars = ("g_engy", "de_over_e", "energy", "en_label")

scpot_method =\
    {-1: "Invalid", 0: "Manual", 1: "pot_swelpw", 2: "pot_swepos",
     3: "pot_sweneg", 4: "pot_sta", 5: "pot_sweshdw"}


#########################################
#           SWEA routines               #
#########################################


def read(data_file, dataset_type="None",
         fields=None, spectra_type="eflux", lib='cdflib',
         include_unit=True, source=('comp',)):

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
        - scpot: spacecraft potential as a function of time,
          combined from multiple datasources.

    data_file: name of file to be read
    """

    if dataset_type == "None":
        # Cut the path down to just the filename:
        filename = os.path.split(data_file)[1]
        if "spec" in filename:
            dataset_type = "spec"
        elif "3d" in filename:
            dataset_type = "3d"
        elif "pad" in filename:
            dataset_type = "pad"
        elif "scpot" in filename:
            dataset_type = "scpot"

    if dataset_type == "scpot":
        return read_scpot(data_file, include_unit=include_unit, source=source)

    if not fields:
        if dataset_type == "spec":
            # spectra
            fields = ("time_unix", "epoch",
                      "diff_en_fluxes", "counts",
                      "energy", "de_over_e",
                      "num_accum", "quality", "secondary")

        if dataset_type == "pad":
            # PAD
            fields = ("time_unix", "epoch",
                      "diff_en_fluxes", "counts",
                      "energy", "de_over_e",
                      "pa", "d_pa",
                      "binning", "quality")

        if dataset_type == "3d":
            # 3D
            fields = ("time_unix", "epoch",
                      "diff_en_fluxes", "counts",
                      "energy", "elev", "azim",
                      "binning", "quality")

    # Pull data from the SWIA Level 2 CDF
    data_cdf = read_cdf(data_file, fields, lib=lib)
    # print(data_cdf.keys())
    # for n in data_cdf:
    #     print(n, data_cdf[n].shape)
    # input()

    # All energy-dependent axes are ordered from highest to lowest
    # since the instrument works by cycling through the energies
    # in that order. However, for our uses, it's best to flip those
    # axes so they are increasing.
    data_cdf = invert_energy_axis(
        data_cdf, energy_names=(i for i in fields if i in e_vars),
        energy_dependent_var_names=(i for i in fields if i in edependent_vars),
        energy_dependent_var_axis=-1)

    # Option to include unit:
    if include_unit:
        data_cdf = process_data_dict(data_cdf, units=swea_units)

    return data_cdf


#########################################
#      Spacecraft potential             #
#########################################


def read_scpot(filename, source=('comp',), include_unit=True):
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

        if include_unit:
            scpot_struct[name] = (scpot[name], unit)
        else:
            scpot_struct[name] = scpot[name]

    return scpot_struct
