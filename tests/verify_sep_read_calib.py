import os
import datetime as dt

import numpy as np

import cdflib

import matplotlib.dates as mdates
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm

from mavenpy import sep
from mavenpy import helper


def datetime64_to_datetime(d):

    return dt.datetime.strptime(np.datetime_as_string(d,unit='s'), '%Y-%m-%dT%H:%M:%S')


def cdf_todatetime_error(filename):

    data_cdf = cdflib.CDF(filename)
    file = os.path.split(raw_filename)[-1]
    data_i = data_cdf.varget("epoch")
    time_unx = data_cdf.varget("time_unix")

    unfiltered = cdflib.cdfepoch.to_datetime(data_i)
    filtered_time_unx = time_unx[np.where(data_i > 0)[0]]
    filtered = cdflib.cdfepoch.to_datetime(data_i[np.where(data_i > 0)[0]])

    plt.figure(figsize=(10.5, 5))
    plt.title(file)
    plt.plot(time_unx, unfiltered,
             label='CDF_EPOCH = cdf.varget("epoch")')
    plt.plot(filtered_time_unx, filtered,
             label='CDF_EPOCH = cdf.varget("epoch")[epoch > 0]')
    plt.xlabel("Unix time (elapsed s since 1970-01-01/00:00 w/o leap seconds")
    plt.ylabel("cdflib.cdfepoch.to_datetime(CDF_EPOCH[TT2000])")
    plt.legend()
    # plt.yticks()
    plt.gca().yaxis.set_major_formatter(mdates.DateFormatter('%b-%d %H:%M'))
    # plt.show()


def compare_cdf_lib(raw_filename):

    libs = ('cdflib', 'spacepy')

    detectors = ['f']
    telescopes = ['a', 'b']
    plot_t_dp = 'A-F'

    fig, ax = plt.subplots(nrows=2, sharex=True, sharey=True)
    fig_0, ax_0 = plt.subplots(nrows=1)

    epoch_array = []
    time_unx_array = []

    for i, lib_i in enumerate(libs):
        raw = sep.read_raw(raw_filename, lib=lib_i,
                           dp_ids=detectors, t_ids=telescopes)

        # print(raw.keys())
        # print(raw["epoch"][:100])
        # print(raw["time_unix"][:10])
        epoch_i = raw["epoch"]
        time_unx_i = raw["time_unix"]
        en_i = raw["{}_energy".format(plot_t_dp)]
        rates_i = raw["{}_rates".format(plot_t_dp)]
        flux_i = raw["{}_flux".format(plot_t_dp)]
        eflux_i = en_i[np.newaxis, :]*flux_i

        p = ax[i].pcolormesh(epoch_i, en_i, rates_i.T,
                             norm=LogNorm())
        plt.colorbar(p)
        ax[i].set_yscale('log')

        colors = cm.rainbow(np.linspace(0, 1, len(en_i)))

        print(lib_i, 'epoch[0] =', epoch_i[0])
        print(lib_i, 'unx[0] =', time_unx_i[0])
        print(lib_i, 'epoch[-1] =', epoch_i[-1])
        print(lib_i, 'unx[-1] =', time_unx_i[-1])
        print("Energy: ", en_i)
        print("Eflux @ t = 0: ", eflux_i[0, :])

        # for j in range(len(en_i)):
        j = -1

        if lib_i == 'cdflib':
            ls_i = '-'
        else:
            ls_i = ':'

        ax_0.plot(epoch_i, rates_i[:, j],
                  color=colors[j], linestyle=ls_i,
                  label=helper.format_energy_as_string(en_i[j]))
        epoch_array.append(epoch_i)
        time_unx_array.append(time_unx_i)

    epoch_cdflib = epoch_array[libs.index("cdflib")]
    epoch_spacepy = epoch_array[libs.index("spacepy")]

    epoch_cdflib_dt = [datetime64_to_datetime(i) for i in epoch_cdflib]
    epoch_spacepy_dt64 = epoch_spacepy.astype("datetime64[ns]")

    delta_dt = [(i - j).total_seconds() for (i, j) in zip(epoch_spacepy, epoch_cdflib_dt)]
    delta_dt64 = (epoch_spacepy_dt64 - epoch_cdflib).astype('float')*1.e-9
    print(delta_dt[:10])
    print(delta_dt64[:10])

    plt.figure()
    plt.plot(time_unx_i, delta_dt)
    plt.plot(time_unx_i, delta_dt64)
    # plt.plot(time_unx_i, epoch_i)
    # plt.plot(time_unx_i + 100, raw["epoch"].astype("datetime64[ns]"))
    plt.ylabel("Epoch_Spacepy - Epoch_cdflib, s")
    plt.xlabel("Unix time (elapsed s since 1970-01-01/00:00 w/o leap seconds")
    # plt.show()


def compare_calib_file_v_calc(raw_file, cal_file, detector, data_name):

    calib_data_name = "{}_{}".format(detector, data_name)

    cal = sep.read_cal(cal_file, lib='cdflib', energy_unit='keV')
    print("Cal from file keys:", cal.keys())

    # Load from raw
    raw = sep.read_raw(raw_file, lib='cdflib')

    # Convert into calibrated data
    calib_from_raw = sep.raw_to_calibrated(
        raw, look_directions=("f", "r"), particle=("elec", "ion"), unit="flux")
    print("Cal from raw keys:", calib_from_raw.keys())

    # Get data for comparison
    epoch_local = calib_from_raw["epoch"]
    local_calib = calib_from_raw[data_name]
    raw_atten = calib_from_raw["attenuator_state"]
    prefix = "_".join((data_name.split("_")[:2]))
    en_i = calib_from_raw["{}_energy".format(prefix)]
    print(en_i)

    epoch_file = cal["1_epoch"][0]
    calib_atten = cal["1_atten"][0]
    file_calib = cal[calib_data_name][0]
    calib_unit = cal[calib_data_name][1]

    # Compare flux dims
    print('Shape flux (calib from raw):', local_calib.shape)
    print('Shape flux (calib from file):', file_calib.shape)

    # Make plot
    fig, ax = plt.subplots(nrows=2, sharex=True)
    ax[0].plot(epoch_local, local_calib[:, 2])
    ax[0].plot(epoch_file, file_calib[:, 2])
    # ax[0].set_yscale('log')
    ax[0].set_ylabel("Flux, {}".format(calib_unit))

    ax[1].plot(epoch_local, raw_atten, label='calibrated from raw')
    ax[1].plot(epoch_file, calib_atten, label='read from calibrated')
    ax[1].set_ylabel("Attenuator")
    ax[1].legend()

    colors = cm.rainbow(np.linspace(0, 1, len(en_i)))
    fig, ax = plt.subplots(nrows=2, sharex=True)
    for j in range(len(en_i)):

        if j == 0:
            ax[0].plot(epoch_local, local_calib[:, j],
                       color=colors[j], label='calibrated from raw',
                       marker='.')
            ax[0].plot(epoch_file, file_calib[:, j],
                       color=colors[j], label='read from calibrated')
        else:
            ax[0].plot(epoch_local, local_calib[:, j],
                       color=colors[j],
                       marker='.')
            ax[0].plot(epoch_file, file_calib[:, j],
                       color=colors[j])

        # delta = flux_cal_from_raw[:, j] - flux_from_cal[:, j]
        # indx = np.where(delta != 0)[0]

        # ax_0[0].plot(epoch_i[indx], delta[indx]/flux_cal_from_raw[indx, j], color=colors[j], marker='x')

        # ax_0[0].plot(epoch_i[indx], flux_cal_from_raw[indx, j], color=colors[j], marker='x')
        # ax_0[0].plot(epoch_i[indx], flux_from_cal[indx, j], color=colors[j], marker='.')

    # ax_0[0].set_yscale('log')
    ax[1].plot(calib_from_raw["epoch"], calib_from_raw["attenuator_state"], color='gray', marker='.')
    ax[1].plot(cal["1_epoch"][0], cal["1_atten"][0], color='gray')

    ax[1].plot(epoch_local, raw_atten, label='calibrated from raw',
               color='gray', marker='.')
    ax[1].plot(epoch_file, calib_atten,
               label='read from calibrated')
    ax[1].set_ylabel("Attenuator")
    ax[0].legend()

    # plt.show()

if __name__ == "__main__":

    raw_filename = "/Users/rjolitz/Downloads/test_maven_data/mvn_sep_l2_s1-raw-svy-full_20220801_v04_r03.cdf"
    cal_filename = "/Users/rjolitz/Downloads/test_maven_data/mvn_sep_l2_s1-cal-svy-full_20220801_v04_r03.cdf"

    # Routine to show if a NaN-riddled time axis
    # is converted to UTC time properly
    cdf_todatetime_error(raw_filename)

    # Compare spacepy and cdflib read functions
    compare_cdf_lib(raw_filename)

    # compare calibrated and raw -> calibrated data
    compare_calib_file_v_calc(raw_filename, cal_filename, "1", "r_ion_flux")
    plt.show()
