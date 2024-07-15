import os

import numpy as np


l2_csn_dtype = [('epoch_utc', 'U19'),
                ('time_unix', '<f8'),
                ('t_sclk', '<f8'),
                ("t_tid", "f8"),
                ("tid", "u4"),
                ("orbit", "u4"),
                ("focusmode", "U3"),
                ("alt", "<f8"),
                ("lst", "<f8"),
                ("long", "<f8"),
                ("lat", "<f8"),
                ("sza", "<f8"),
                ("mass", "u4"),
                ("species", "U3"),
                ("cps_dt_bkd", "f8"),
                ("abundance", "f8"),
                ("precision", "f8"),
                ("quality", "U2")]


l2_ion_dtype = [('epoch_utc', 'U19'),
                ('time_unix', '<f8'),
                ('t_sclk', '<f8'),
                ("t_tid", "f8"),
                ("tid", "u4"),
                ("orbit", "u4"),
                ("focusmode", "U3"),
                ("alt", "<f8"),
                ("lst", "<f8"),
                ("long", "<f8"),
                ("lat", "<f8"),
                ("sza", "<f8"),
                ("ion_mass", "u4"),
                ("cps_dt", "f8"),
                ("abundance", "f8"),
                ("sensitivity", "f8"),
                ("SC_potential", "f8"),
                ("precision", "f8"),
                ("quality", "U2")]


l3_rsn_dtype = [('epoch_utc', 'U19'),
                ('time_unix', '<f8'),
                ('t_sclk', '<f8'),
                ("t_tid", "f8"),
                ("tid", "u4"),
                ("orbit", "u4"),
                ("focusmode", "U3"),
                ("alt", "<f8"),
                ("mass", "u4"),
                ("species", "U3"),
                ("density_bins", "f8"),
                ("quality", "U2")]

l3_rsh_dtype = [('epoch_utc', 'U19'),
                ('time_unix', '<f8'),
                ('t_sclk', '<f8'),
                ("t_tid", "f8"),
                ("tid", "u4"),
                ("orbit", "u4"),
                ("exo_alt", "f8"),
                ("mass", "u4"),
                ("species", "U3"),
                ("scale_height", "f8"),
                ("scale_height_error", "f8"),
                ("temperature", "f8"),
                ("temperature_error", "f8"),
                ("fit_residual", "f8"),
                ("quality", "U2")]


def filename_to_datatype(name):
    split_str = name.split("_")
    level = split_str[2]
    datatype = "-".join(split_str[3].split("-")[:2])

    return level, datatype


def read(filename, dataset_type=None, fields=None, include_unit=None, remove_errant_data=True):

    '''Routine to read Level 2 and 3 CSV files
    produced by the NGIMS team

    Note: Level 3 has significant errors and doesn't make any sense.
        Disrecommend use.

    dataset_type: string, can be 'csn-abund',
        'ion-abund', 'res-den', 'res-sht'.
        NOTE: res-den is not currently in use and DISRECOMMENDED
        for use.
    remove_errant_data: a keyword that when True, removes all rows containing
        -999 as abundance
      '''

    if not dataset_type:
        level, dataset_type = filename_to_datatype(os.path.split(filename)[1])

    # print(filename, level, dataset_type)

    if dataset_type == 'csn-abund':
        dtype_i = l2_csn_dtype
    elif dataset_type == 'ion-abund':
        dtype_i = l2_ion_dtype
    elif dataset_type == 'res-den':
        dtype_i = l3_rsn_dtype
    elif dataset_type == 'res-sht':
        dtype_i = l3_rsh_dtype

    # Retrieve all fields
    if not fields:
        fields = [i[0] for i in dtype_i]

    # genfromtxt is faster (0.02 s for 1000x iterations,
    # v 3.7 s for 1000x iterations for loadtxt)
    data = np.genfromtxt(filename, delimiter=',', skip_header=1, dtype=dtype_i)
    # data = np.loadtxt(filename, delimiter=',', skiprows=1, dtype=dtype_i)

    if remove_errant_data:
        not_errant = np.where(data['abundance'] != -999.0)[0]
        # print(data.shape, not_errant.shape)
        data = data[not_errant]
        # input()

    return data
