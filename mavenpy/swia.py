import os

import numpy as np


from .read import read_cdf
from . import helper


swia_units =\
    {"time_unix": "s, POSIX seconds after 1970", "epoch": "utc",
     "density": "#/cm3",
     "velocity": "km/s, S/c frame",
     "velocity_mso": "km/s, MSO frame",
     "quality_flag": "unitless",
     "decom_flag": "unitless",
     "atten_state": "unitless, close (2) or open (1)",
     "num_accum": "# of 4-sec accumulations per observation",
     "energy_spectra": "eV",
     "energy_coarse": "eV",
     "energy_fine": "eV",
     "de_over_e_spectra": "unitless",
     "de_over_e_coarse": "unitless",
     "de_over_e_fine": "unitless",
     "counts": "#",
     "diff_en_fluxes": "eV/(cm2 s sr eV)",
     "spectra_diff_en_fluxes": "eV/(cm2 s sr eV)",
     "theta_spectra": "deg., instrument deflection angle when attenuator open",
     "theta_coarse": "deg., instrument deflection angle when attenuator open",
     "theta_fine": "deg., instrument deflection angle when attenuator open",
     "theta_atten_spectra": "deg., instrument deflection"
        " angle when attenuator closed",
     "theta_atten_coarse": "deg., instrument deflection"
        " angle when attenuator closed",
     "theta_atten_fine": "deg., instrument deflection"
        " angle when attenuator closed",
     "phi_spectra": "deg., instrument",
     "phi_coarse": "deg., instrument",
     "phi_fine": "deg., instrument",
     "estep_first": "index",
     "dstep_first": "index"
     }


#########################################
#           SWIA routines               #
#########################################


def uncompress(uncompressed_dim, index_arrs, compressed_arr,
               iterate_over_time=None, iterate_over_index=None):

    '''Expand a compressed array that has a time-varying
    start index in a larger array.
    uncompressed_dim: tuple of lengths of
      each uncompressed array
    index_arrs: tuple of index arrays
    compressed_arr: n_time x n_narrow
    '''
    N = compressed_arr.shape[0]
    compressed_dim = compressed_arr.shape[1:]
    # print(uncompressed_arr.shape, compressed_arr.shape)

    if iterate_over_time:
        # Basic for loop:
        # Iterating over time
        uncompressed_arr = np.zeros(shape=(N, *uncompressed_dim))
        for i in range(N):
            index_i = [i]
            for index_arr, len_i, in zip(index_arrs, compressed_dim):
                if index_arr is "None":
                    index_i.append(slice(None))
                else:
                    i_i = index_arr[i]
                    i_f = i_i + len_i
                    index_i.append(slice(i_i, i_f))
            uncompressed_arr[tuple(index_i)] = compressed_arr[i, ...]
    elif iterate_over_index:
        # Inverted for loop:
        uncompressed_arr = np.zeros(shape=(N, *uncompressed_dim))

        non_none_index = [i for i in index_arrs if i is not "None"]
        non_none_dim =\
            [j for (i, j) in zip(index_arrs, uncompressed_dim)
             if i is not "None"]

        aggregate_index_arr = non_none_index[0]*1
        for n_i, index_arr in zip(non_none_dim[:-1], non_none_index[1:]):
            # print(n_i, max(index_arr), min(index_arr))
            # input()
            aggregate_index_arr = n_i*aggregate_index_arr + index_arr
            # print(aggregate_index_arr.shape, min(aggregate_index_arr),
            #       max(aggregate_index_arr))

        unique_index = np.unique(aggregate_index_arr)

        for index_i in unique_index:
            time_index_i = np.where(aggregate_index_arr == index_i)[0]
            N_time_index = len(time_index_i)

            if N_time_index != 0:

                subset_i = compressed_arr[time_index_i, ...]

                index = [time_index_i]

                for index_arr, n_c in zip(index_arrs, compressed_dim):
                    if index_arr is "None":
                        slice_i = slice(None)
                    else:
                        n_index_arr = index_arr[time_index_i]
                        n_index = n_index_arr[0]
                        slice_i = slice(n_index, n_index + n_c)
                    index.append(slice_i)

                index = tuple(index)
                uncompressed_arr[index] = subset_i

    return uncompressed_arr


def read(file_path, dataset_type="None", fields=None, lib='cdflib',
         preprocess=True, include_unit=True, **kwargs):

    """Returns SWIA data.

    This routine is based on the SWIA SIS document and the MAVEN SPEDAS
    routine 'mvn_swia_load_l2_data.pro'. We thank and credit the
    authors of both.

    dataset_type: name of dataset to be pulled. options:
        - m: moment data, containing plasma moments (ion density, velocity,
            pressure, temperature) calculated onboard
        - f: fine high resolution spectra, 10 small anode (phi) counts/eflux
            for 48 energy steps and 12 deflection (theta)
            steps around the peak flux
        - c: coarse low resolution spectra, summed over anodes (phi),
            deflection steps (theta), and energy bins to produce counts/eflux
            for 16 phi bins x 48 energy bins x 4 deflection bins
        - s: onboard 1D spectra, not weighted for reduced geometric factor at
            high deflection angles. should not be used for science.
    sum_over_fov: True if summing over the deflection/anode angles

    data_file: name of file to be read
    swia_qlevel: optional, floating point, quality flag for SWIA cutoff.

    Returns dictionary of time, Np, Vsc, Vmso corresponding
    to requested dataset.
    """

    if dataset_type == "None":
        filename = os.path.split(file_path)[1]
        if "onboardsvyspec" in filename:
            dataset_type = "s"
        elif "onboardsvymom" in filename:
            dataset_type = "m"
        elif "coarse" in filename:
            dataset_type = "c"
        elif "fine" in filename:
            dataset_type = "f"

    if not fields:
        if dataset_type == "m":
            # Moment
            fields = ("time_unix", "density", "velocity", "velocity_mso",
                      "epoch", "quality_flag", "decom_flag", "atten_state")

        if dataset_type == "s":
            # spectra
            fields = ("time_unix", "epoch", "num_accum", "energy_spectra",
                      "decom_flag", "spectra_diff_en_fluxes", "spectra_counts",
                      "de_over_e_spectra")

        if dataset_type == "c":
            # coarse
            fields = ("time_unix", "epoch", "num_accum", "energy_coarse",
                      "diff_en_fluxes", "counts", "de_over_e_coarse")
            # fields = ("time_unix", "epoch", "num_accum", "energy_coarse",
            #           "atten_state"
            #           "theta_coarse", "theta_atten_coarse", "phi_coarse",
            #           "diff_en_fluxes", "counts", "de_over_e_coarse")

        if dataset_type == "f":
            # fine
            fields = ("time_unix", "epoch", "energy_fine",
                      "diff_en_fluxes", "counts",  "estep_first")
            # fields = ("time_unix", "epoch", "energy_fine",
            #           "atten_state"
            #           "theta_fine", "theta_atten_fine", "phi_fine",
            #           "diff_en_fluxes", "counts", "de_over_e_fine")

    # Pull data from the SWIA Level 2 CDF
    data = read_cdf(file_path, fields, lib=lib)

    if preprocess:
        data = process(data, dataset_type, fields=fields, **kwargs)

    # Assign unit
    if include_unit:
        data = helper.process_data_dict(data, units=swia_units)
    return data


def process(data_cdf, dataset_type, fields=None, swia_qlevel=0.5,
            sum_over_fov=False, avg_over_fov=True,
            uncompress_fine=False):

    # Times reported in L0 and L2 files correspond
    # to the start time of accumulation. Thus,
    # the more appropriate time should be plotted
    # according to the center time of the accumulation
    # (and only start_time + 2 seconds for products
    # that aren't accumulated over multiple sweeps like
    # fine and moment).

    if not fields:
        fields = list(data_cdf)

    if "m" in dataset_type or "f" in dataset_type:
        ctime_unx = data_cdf["time_unix"] + 2.0
    else:
        ctime_unx = data_cdf["time_unix"] + 4.0*data_cdf["num_accum"]/2

    # Data should be filtered to remove NaNs in the time
    # axis (induced by changing attenuator state),
    # in addition to low quality flag
    # and low decom_flag (combination of attenuator
    # and telemetry state). This is done in analogy
    # to the SPEDAS routine.
    non_nan_index = (~np.isnan(ctime_unx))

    if dataset_type == "m":
        acceptable_index = (
            (data_cdf["decom_flag"] >= swia_qlevel) &
            ((data_cdf["quality_flag"] >= swia_qlevel) &
             non_nan_index))
    elif dataset_type == "s":
        # If recovering onboard S/C spectra from SWIA (SWIS),
        # then we only need to remove low quality datapoints.
        acceptable_index = (
            (data_cdf["quality_flag"] >= swia_qlevel) &
            non_nan_index)
    else:
        acceptable_index = non_nan_index

    data_cdf["time_unix"] = ctime_unx

    data_cdf = helper.process_data_dict(data_cdf, conditional_array=acceptable_index)

    # SWIA efluxes and counts are a 4D structure [time, theta, phi, energy]
    # If requested, can sum over all theta & phi bins to find
    # total counts / average eflux.

    if uncompress_fine:
        energy = data_cdf["energy_fine"]
        estep_first = data_cdf["estep_first"].astype('int')
        dstep_first = data_cdf['dstep_first'].astype('int')
        phi_first = np.zeros(shape=data_cdf["time_unix"].shape).astype('int')
        theta = data_cdf['theta_fine']
        phi = data_cdf['phi_fine']

    if dataset_type == "c":
        norm = 64
    elif dataset_type == "f":
        norm = 120

    # Retrieve the fields corresponding to the
    # 4D or 2D arrays
    flux_names = [name for name in fields if
                  "diff_en_fluxes" in name or "count" in name]

    for name_arr_4d in flux_names:
        data_i = data_cdf[name_arr_4d]
        # Sum over the theta & fine bins, divide by the number
        # of theta and phi bins if requested:
        if sum_over_fov or avg_over_fov:
            data_i = np.sum(data_i, axis=(1, 2))
        if avg_over_fov:
            data_i = data_i / norm

        # SWIA fine spectra index is timevarying,
        # so if want to uncompress on this time step
        # can do so.
        if dataset_type == "f" and uncompress_fine:
            dim_all = (len(energy),)
            start_index = (estep_first,)
            if not (sum_over_fov and avg_over_fov):
                dim_all = (len(phi), len(theta), len(energy))
                start_index = (phi_first, dstep_first, estep_first)
            data_i = uncompress(
                dim_all, start_index, data_i, iterate_over_index=True)

        # Reassign the data_cdf variable
        data_cdf[name_arr_4d] = data_i

    # if uncompress_fine:
    #     del data_cdf["estep_first"]
    #     del data_cdf["dstep_first"]

    # By default, energy-varying data in SWIA
    # are arranged in order of decreasing energy.
    # **REVERSES** the direction of the energy axis so
    # energies is an increasing instead of decreasing array.
    theta_containing_axis = [i for i in fields if "theta" in i]
    evarying_axes = flux_names + theta_containing_axis

    data_cdf = helper.invert_energy_axis(
        data_cdf, energy_dependent_var_axis=-1,
        energy_dependent_var_names=evarying_axes)

    # print(final_data)

    return data_cdf
