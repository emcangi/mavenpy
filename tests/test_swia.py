import time
import argparse

import datetime as dt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from mavenpy import file_path, helper, retrieve, load, swia


def dimension_check(instrument_name, axis_data, data, axis_index=0):

    assert axis_data.shape[0] == data.shape[axis_index],\
        "{} dimensions ({}) don't match epoch dimensions ({})".format(
            instrument_name, data.shape, axis_data.shape)


# Load two days of instrument data for a given target_date and confirm that the
# dimensions are correct.
# Remember that the output of each data load
# function is a dictionary, where each key is a variable name from the original
# CDF/IDL sav/NPZ, and the corresponding value is a two-element tuple:
# the first element is the variable data (e.g. UTC times for the
# variable "epoch"), and the second element is a string defining the unit.

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data_directory", help="Directory containing MAVEN data.",
        required=True)
    parser.add_argument(
        "--start_date", help="Start date of retrieved data (YYYY-MM-DD).",
        type=lambda s: dt.datetime.strptime(s, "%Y-%m-%d"), required=True)
    # Download args.
    parser.add_argument(
        "--download", help="Download files to data_directory if not found.",
        action="store_true", default=False)
    parser.add_argument(
        "--verbose", help="Prints verbose logging of download.",
        action="store_true", default=False)
    parser.add_argument(
        "--remote", help="Remote to download data from (default: ssl_sprg).",
        default="ssl_sprg")
    parser.add_argument(
        "--local_source_tree",
        help="Directory tree of remote that is mirrored on local.",
        default="ssl_sprg")
    parser.add_argument("--username", help="Username for data download.")
    parser.add_argument("--password", help="Password for data download.")
    args = parser.parse_args()

    start = args.start_date
    data_dir = args.data_directory
    n_days = 2

    swia_datanames =\
        ("onboardsvymom", "coarsesvy3d", "coarsearc3d",
         "finesvy3d", "finearc3d")
    # swia_datanames =\
    #     ("onboardsvymom",)
    eflux_zlim = [1e4, 1e8]

    fig, ax = plt.subplots(
        nrows=len(swia_datanames), sharex=True, figsize=(8.5, 7.5))

    for i, name_i in enumerate(swia_datanames):
        print(name_i)

        if args.download:
            print("Checking for file updates / download if not found...")
            init_time = time.time()
            retrieve.sdc_retrieve(
                'swi', destination_dir=data_dir,
                start_date=start, n_days=n_days, verbose=args.verbose,
                mirrored_local_source_tree=args.local_source_tree,
                source=args.remote,
                username=args.username, password=args.password,
                dataset_name=name_i, ext='cdf')
            print("Done, took {} s.".format(time.time() - init_time))

        # Get the local files:
        swi_files = file_path.local_file_names(
            data_dir, 'swi', start_date=start, n_days=n_days,
            source=args.local_source_tree, dataset_name=name_i, ext='cdf')
        print("File names: ", swi_files)

        # Retrieve the data:
        print("Loading the data...")
        init_time = time.time()
        swi = load.load_data(swi_files, include_unit=True, swia_qlevel=0.5)
        print("Done, took {} s.".format(time.time() - init_time))

        # print(swi.keys())
        swi_epoch = swi["epoch"][0]

        if name_i == "onboardsvymom":
            n, n_unit = swi["density"]
            vxyz, v_unit = swi["velocity_mso"]
            dimension_check("Densities", swi_epoch, n, axis_index=0)
            dimension_check("Velocity_MSO", swi_epoch, vxyz, axis_index=0)
            vmag = np.sqrt(np.sum(vxyz**2, axis=1))

            ax[i].plot(swi_epoch, n, color='k')
            ax[i].set_yscale('log')
            ax[i].set_ylabel('n, {}'.format(n_unit))
            ax2 = ax[i].twinx()
            ax2.plot(swi_epoch, vmag, color='b')
            ax2.set_ylabel('|v|, {}'.format(v_unit))

        else:
            energy_name = [n for n in swi if "energy" in n][0]
            energy = swi[energy_name][0]
            diff_en_flux, diff_en_unit = swi["diff_en_fluxes"]

            print('diff en flux dims:', diff_en_flux.shape)
            print('epoch dims:', swi_epoch.shape)
            print('energy dims:', energy.shape)

            dimension_check("SWIA diff Eflux",
                            swi_epoch, diff_en_flux,
                            axis_index=0)

            if "fine" in name_i:
                estep_first = swi["estep_first"][0]
                diff_en_flux = swia.uncompress(
                    (len(energy),),  (estep_first,),
                    diff_en_flux, iterate_over_index=True)

            dimension_check("SWIA diff Eflux (energy)",
                            energy, diff_en_flux,
                            axis_index=-1)
            p = ax[i].pcolormesh(
                swi_epoch, energy, diff_en_flux.T,
                norm=LogNorm(vmin=eflux_zlim[0], vmax=eflux_zlim[1]))
            ax[i].set_yscale('log')
            ax[i].set_ylabel("{}\nenergy".format(name_i))
            helper.add_colorbar_outside(
                p, fig, ax[i], label='Eflux\n{}'.format(diff_en_unit))

    plt.show()
