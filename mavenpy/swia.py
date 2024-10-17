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
         include_unit=True, preprocess=True, over_fov="sum", swia_qlevel=0.5,
         uncompress_fine=True):

    """Returns SWIA data.

    This routine is based on the SWIA SIS document and the MAVEN SPEDAS
    routine 'mvn_swia_load_l2_data.pro'. We thank and credit the
    authors of both.

    file_path: SWIA file to be read.
    dataset_type: name of dataset to be pulled. options:
        - onboardsvymom: moment data, containing plasma moments
            (ion density, velocity, ressure, temperature) calculated onboard
        - onboardsvyspec: onboard 1D spectra, not weighted for
            reduced geometric factor at high deflection angles.
            should not be used for science.
        - fine(svy/arc)3d: fine high resolution spectra, 10 small
            node (phi) counts/eflux for 48 energy steps and 12
            deflection (theta) steps around the peak flux
        - coarse(svy/arc)3d: coarse low resolution spectra,
            summed over anodes (phi),
            deflection steps (theta), and energy bins to produce counts/eflux
            for 16 phi bins x 48 energy bins x 4 deflection bins

    over_fov: "sum" if summing over the deflection/anode angles,
        "avg" if averaging over, "None" if not doing anything.
    preprocess: Boolean, if on will remove NaNs, recenter time axis,
        uncompress fine data product, etc.
    swia_qlevel: optional, floating point, quality flag for SWIA cutoff.
    """

    if dataset_type == "None":
        filename = os.path.split(file_path)[1]
        dataset_type = filename.split("_")[3]

    if not fields:
        if dataset_type == "onboardsvymom":
            # Moment
            fields = ("time_unix", "density", "velocity", "velocity_mso",
                      "epoch", "quality_flag", "decom_flag", "atten_state")

        elif dataset_type == "onboardsvyspec":
            # spectra
            fields = ("time_unix", "epoch", "num_accum", "energy_spectra",
                      "decom_flag", "spectra_diff_en_fluxes", "spectra_counts",
                      "de_over_e_spectra")

        elif "coarse" in dataset_type:
            # coarse
            if over_fov == "avg" or over_fov == "sum":
                fields = ("time_unix", "epoch", "num_accum", "energy_coarse",
                          "diff_en_fluxes", "counts", "de_over_e_coarse")
            else:
                fields = ("time_unix", "epoch", "num_accum", "energy_coarse",
                          "atten_state", "theta_coarse", "theta_atten_coarse",
                          "phi_coarse", "diff_en_fluxes", "counts",
                          "de_over_e_coarse")

        elif "fine" in dataset_type:
            # fine
            if over_fov == "avg" or over_fov == "sum":
                fields = ("time_unix", "epoch", "energy_fine",
                          "diff_en_fluxes", "counts",  "estep_first")
            else:
                fields = ("time_unix", "epoch", "energy_fine", "atten_state",
                          "theta_fine", "theta_atten_fine", "phi_fine",
                          "diff_en_fluxes", "counts", "de_over_e_fine",
                          "estep_first", "dstep_first")
        else:
            raise ValueError(
                "No dataset of name '{}' recognized, use instead: "
                "onboardsvyspec, onboardsvymom, coarse(svy/arc)3d, "
                "or fine(svy/arc3d).")

    # Pull data from the SWIA Level 2 CDF
    data = read_cdf(file_path, fields, lib=lib)

    if preprocess:
        data = process(data, dataset_type, fields=fields,
                       over_fov=over_fov,
                       uncompress_fine=uncompress_fine,
                       swia_qlevel=swia_qlevel)

    # Assign unit
    if include_unit:
        data = helper.process_data_dict(data, units=swia_units)
    return data


def process(data_cdf, dataset_type, fields=None, swia_qlevel=0.5,
            over_fov="sum", uncompress_fine=False):

    # Times reported in L0 and L2 files correspond
    # to the start time of accumulation. Thus,
    # the more appropriate time should be plotted
    # according to the center time of the accumulation
    # (and only start_time + 2 seconds for products
    # that aren't accumulated over multiple sweeps like
    # fine and moment).

    if not fields:
        fields = list(data_cdf)

    if dataset_type == "onboardsvymom" or "fine" in dataset_type:
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

    if dataset_type == "onboardsvymom":
        acceptable_index = (
            (data_cdf["decom_flag"] >= swia_qlevel) &
            ((data_cdf["quality_flag"] >= swia_qlevel) &
             non_nan_index))
    elif dataset_type == "onboardsvyspec":
        # If recovering onboard S/C spectra from SWIA (SWIS),
        # then we only need to remove low quality datapoints.
        acceptable_index = (
            (data_cdf["quality_flag"] >= swia_qlevel) &
            non_nan_index)
    else:
        acceptable_index = non_nan_index

    data_cdf["time_unix"] = ctime_unx

    data_cdf = helper.process_data_dict(
        data_cdf, conditional_array=acceptable_index)

    # SWIA efluxes and counts are a 4D structure [time, theta, phi, energy]
    # If requested, can sum over all theta & phi bins to find
    # total counts / average eflux.

    if "fine" in dataset_type and uncompress_fine:
        energy = data_cdf["energy_fine"]
        estep_first = data_cdf["estep_first"].astype('int')
        if not (over_fov == "sum" or over_fov == "avg"):
            dstep_first = data_cdf['dstep_first'].astype('int')
            phi_first = np.zeros(shape=data_cdf["time_unix"].shape).astype('int')
            theta = data_cdf['theta_fine']
            phi = data_cdf['phi_fine']

    # No need to define norm if not either coarse or fine,
    # since there are no 4d arrays in the other two:
    if "coarse" in dataset_type:
        norm = 64
    elif "fine" in dataset_type:
        norm = 120

    # Retrieve the fields corresponding to the
    # 4D or 2D arrays
    flux_names = [name for name in fields if
                  "diff_en_fluxes" in name or "count" in name]

    for name_arr_4d in flux_names:
        data_i = data_cdf[name_arr_4d]
        # Sum over the theta & fine bins, divide by the number
        # of theta and phi bins if requested:
        if over_fov == "sum" or over_fov == "avg":
            data_i = np.sum(data_i, axis=(1, 2))
            if over_fov == "avg":
                data_i = data_i / norm

        # SWIA fine spectra index is timevarying,
        # so if want to uncompress on this time step
        # can do so.
        if "fine" in dataset_type and uncompress_fine:
            dim_all = (len(energy),)
            start_index = (estep_first,)
            if not (over_fov == "sum" or over_fov == "avg"):
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
