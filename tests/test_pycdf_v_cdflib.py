import datetime as dt
import time
import numpy as np
from matplotlib import pyplot as plt
import mavenpy as mav


def shape(dataset):

    if "densities" in dataset:
        data_shape = dataset["densities"][0].shape
    elif "diff_en_flux" in dataset:
        data_shape = dataset["diff_en_flux"][0].shape
    elif "1_f_elec_flux" in dataset:
        data_shape = dataset["1_f_elec_flux"][0].shape
    elif "spectra" in dataset:
        data_shape = dataset["spectra"][0].shape

    return data_shape


def time_load(function, **kwargs):

    # Read w/ cdflib
    start = time.time()
    cdflib_data = mav.load_maven_data(
        data_dir, target_date, n_days, function,
        **kwargs, lib="cdflib")
    cdflib_dur = np.around(time.time() - start, 1)
    print("Cdflib: ", cdflib_dur, "sec; N shape: ", shape(cdflib_data))

    # Read w/ spacepy
    start = time.time()
    spacepy_data = mav.load_maven_data(
        data_dir, target_date, n_days, function,
        **kwargs, lib="spacepy")
    spacepy_dur = np.around(time.time() - start, 1)
    print("Spacepy: ", spacepy_dur, "sec; N shape: ", shape(spacepy_data))


# # Results from 8/10/2023:
# Reading SWIA moment data starting 2014/12/01 for 76 days...
# Cdflib:  0.5 sec; N shape:  (1105038,)
# Spacepy:  3.0 sec; N shape:  (1105038,)
# Reading SWEA spec data starting 2014/12/01 for 76 days...
# Cdflib:  6.8 sec; N shape:  (2062492, 64)
# Spacepy:  9.3 sec; N shape:  (2062492, 64)
# Reading SEP1 data starting 2014/12/01 for 76 days...
# Cdflib:  3.0 sec; N shape:  (792492, 15)
# Spacepy:  3.9 sec; N shape:  (792492, 15)
# Reading EUV data starting 2014/12/01 for 76 days...
# Cdflib:  0.5 sec; N shape:  (109440, 190)
# Spacepy:  0.6 sec; N shape:  (109440, 190)
# # cdflib wins on each count.

data_dir = "/Volumes/isilzha/data"
target_date = dt.datetime(2014, 12, 1)
end_date = dt.datetime(2015, 2, 15)
n_days = int((end_date - target_date).total_seconds() / 60 / 60 / 24)

start_date_str = target_date.strftime("%Y/%m/%d")

print("Reading SWIA moment data starting {} for {} days...".format(
    start_date_str, n_days))
time_load(mav.read_swia_moments, swia_qlevel=0.5)

print("Reading SWEA spec data from {} for {} days...".format(start_date_str, n_days))
time_load(mav.read_swea_l2_spectra)

print("Reading SEP1 data from {} for {} days...".format(start_date_str, n_days))
time_load(mav.read_sep_l2_spectra, particles=("elec", "ion"),
          look_directions=("r", "f"), spectra_type="flux", detector="1")

print("Reading EUV data from {} for {} days...".format(start_date_str, n_days))
time_load(mav.read_euv_l3_spectra)
