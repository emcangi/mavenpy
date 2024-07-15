from .read import read_cdf
from .helper import process_data_dict


# EUV read routines

# Version / revision information for EUV:
# - v## refers to the file version, which changes if CDF fields are modified.
#    As of 8/23/23, it is locked at v04.
# - r## refers to the version of data available for that version.
# want the newest version & revision.

units = {"epoch": "utc", "x": "unx", "flag": "usability flag",
         "maven_sun_distance": "km", "v": "nm",
         "y": "W/m2/nm", "uncertainty": "%",
         "irradiance": "W/m2/nm"}

column_to_field_name = {"time": "x", "wavelength": "v", "spectra": "y"}
field_name_to_column = {v: k for k, v in column_to_field_name.items()}
modified_labels = column_to_field_name.keys()


def read(dataset_filename, lib='cdflib',
         column_names=("time", "epoch", "wavelength", "spectra", "flag"),
         include_unit=True, relabel_columns=True):

    """Loads daily files of EUV Level 3 data.

    lib: string, 'cdflib' or 'spacepy'

    This routine is based on the EUV SIS document and the MAVEN SPEDAS
    routine 'mvn_euv_l3_load.pro'. We thank and credit the authors of both.

    L3 data from EUV consists of outputs of the FISM-M model based on
    EUV calibrated band irradiance. This includes an averaged spectra
    (day or minute) as a function of time, a quality flag,
    and the Mars-sun distance at time of observation.

    wavelength: wavelengths that the spectra is determined on in
        units of nanometer, (190)
    spectra: averaged spectra in units of Watts per square meter per nanometer,
             2D dataset of dimensions (190 x n_time)
    epoch: datetime of observation, (n_time)
    unx_time: Number of seconds elapsed since 1970-01-01/00:00:00, (n_time)
    flag: integer representing data useability
         (0 - good, 1 - occultation, 2 - no pointing, 3 - sun partial in FOV,
          4 - sun not in FOV, 5 - windowed,  6 - eclipse, 7 - spare), (n_time)
    maven_sun_distance: distance between MAVEN and sun in kilometers, (n_time)

    """

    field_names = [column_to_field_name[i] if i in modified_labels
                   else i for i in column_names]
    data = read_cdf(dataset_filename, field_names, lib=lib)

    # Data relabeled according to new column names
    # Include unit if requested:

    data = process_data_dict(
        data, units=(units if include_unit else None),
        alias=(field_name_to_column if relabel_columns else None))

    return data
