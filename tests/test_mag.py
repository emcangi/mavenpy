import os
import sys
import argparse
import time
import itertools
import datetime as dt

import numpy as np
from matplotlib import pyplot as plt
from dateutil.relativedelta import relativedelta

from mavenpy import file_path, retrieve, load, spice, anc, specification, helper


mag_file_name = "MAG {level} {res} {coord} {ext} file"
alternate_coord_name = {"pc": "GEO", "ss": "MSO", "pl": "Payload"}
ylabel = "Magnetic field, nT\n({alt_coord}, {res}, {level})"


# Test to compare loading time / mem size of different Level 1/2
# MAG data (STS, SAV, CDF).

# Findings 8/26/23:
# For loading one file (8/24/2023 or 5/31/2023)
# Load times, from longest to shortest, for L1:
# - L1 sav full (in /l1/sav/full): ~60 sec
# - L1 sts full (w/ genfromtxt): 10 sec
# - L1 sav 1 sec (in /l1/1sec/full): 2.5 sec
# - L1 sts full (w/ loadtxt): 1 sec
# - L1 sav 30 sec (in /l1/sav/30sec): 0.08 sec
# - L1 CDF full (w/ cdflib): 0.05 sec
# - L1 sav full (in /l1_sav/): 0.01 sec

# Load times, from longest to shortest, for L2 (2022/09/30):
# - L2 sav full in /l2/sav/full: 130 sec
# - L2 sav 1 sec in /l2/sav/1sec: 8 sec
# - L2 sts full in /l2 (w/ loadtxt): 2 sec
# - L2 sav 30 sec in /l2/sav/30sec: 0.3 sec
# - L2 sts 1 sec in /l2 w/ loadtxt: 0.07 sec
# Overall, L2 takes longer due to the inclusion
# of fields for s/c position.

# genfromtxt >> loadtxt is expected because
# additional forloop in function.

# Why the huge difference between l1_sav and l1/sav?
# readsav loops over the records in an array.
# l1/sav has n_sample records, while
# l1_sav has 11 records corresponding
# to the array fields.


def dimension_check(instrument_name, axis_data, data, axis_index=0):

    assert axis_data.shape[0] == data.shape[axis_index],\
        "{} dimensions ({}) don't match epoch dimensions ({})".format(
            instrument_name, data.shape, axis_data.shape)


def check_mag(levels, res, coords, exts,
              data_dir, start, n_days,
              download_info_dict={},
              local_source_tree='',
              download_if_not_found=True, make_plot=None,
              testpayloadtoMSOrotation=False):

    for level_i, res_i, coord_i, ext_i in itertools.product(
            levels, res, coords, exts):

        params = {'ext': ext_i, 'res': res_i,
                  'level': level_i, 'coord': coord_i}
        name_i = mag_file_name.format(**params)

        try:
            specification.check_if_dataset_exist(
                'mag', level=level_i, ext=ext_i, coord=coord_i, res=res_i)
        except IOError:
            print(name_i, ": Doesn't exist, moving on to next.")
            continue

        print(name_i, ": Exists.")

        if download_if_not_found:
            # See if available on the remote:
            remote = download_info_dict['source']
            try:
                specification.check_if_dataset_on_remote(
                    'mag', remote, level=level_i, ext=ext_i,
                    auth_present=(download_info_dict['username'] != ''))
            except IOError:
                print(name_i,
                      ": Not available on {}".format(remote))
                continue

            print("Checking for file updates / download if not found...")
            init_time = time.time()
            retrieve.sdc_retrieve(
                'mag', destination_dir=data_dir,
                start_date=start, n_days=n_days, verbose=True,
                mirrored_local_source_tree=local_source_tree,
                **params, **download_info_dict)
            print("Done, took {} s.".format(time.time() - init_time))

        # Get the local files:
        mag_files = file_path.local_file_names(
            data_dir, 'mag', start_date=start, n_days=n_days,
            source=local_source_tree, **params)
        print("File names: ", mag_files)
        print("File size (MB): ",
              [np.around(os.path.getsize(i)*1e-6, 1) for i in mag_files])

        # Retrieve the data:
        print("Loading the data...")
        init_time = time.time()
        mag = load.load_data(
            mag_files, include_unit=True, **params, verbose=True)
        print("Read/Parse time: {} s.".format(time.time() - init_time))
        print("Mem size: ", sys.getsizeof(mag))

        # Get the Bx, By, Bz, and |B| components
        bx, by, bz = mag["Bx"][0], mag["By"][0], mag["Bz"][0]
        bmag = np.sqrt(bx**2 + by**2 + bz**2)
        mag_epoch = mag["epoch"][0]

        # Check that Bx and the time axis have same length
        dimension_check(name_i, mag_epoch, bx, axis_index=0)

        # Test a rotation from payload to MSO coordinates if desired:
        if coord_i == "pl" and testpayloadtoMSOrotation:
            print("Beginning PL -> MSO coordinate conversion...")
            init_time = time.time()
            bmso = spice.bpl_to_bmso(mag_epoch, bx, by, bz)
            print("Done, took {} s.".format(time.time() - init_time))
            bmag_mso = np.sqrt(np.sum(bmso**2, axis=1))

        # Make the plot if wanted:
        if make_plot:
            nrows = 1
            if coord_i == "pl" and testpayloadtoMSOrotation:
                nrows = 2
            fig, ax = plt.subplots(nrows=nrows, sharex=True)
            fig.suptitle(name_i)

            if not (testpayloadtoMSOrotation and coord_i == 'pl'):
                ax = (ax,)

            # Make the appropriate label
            # (relabeled according to the common coordinate system
            # name).
            label_coord_i = alternate_coord_name[coord_i]
            label_i = ylabel.format(
                level=level_i, res=res_i, alt_coord=label_coord_i)

            ax[0].plot(mag_epoch, bx, color='b', label='Bx')
            ax[0].plot(mag_epoch, by, color='g', label='By')
            ax[0].plot(mag_epoch, bz, color='r', label='Bz')
            ax[0].plot(mag_epoch, bmag, color='k', label='|B|')
            ax[0].set_ylabel(label_i)

            if coord_i == "pl" and testpayloadtoMSOrotation:
                ax[1].plot(mag_epoch, bmso[:, 0], color='b', label='Bx')
                ax[1].plot(mag_epoch, bmso[:, 1], color='g', label='By')
                ax[1].plot(mag_epoch, bmso[:, 2], color='r', label='Bz')
                ax[1].plot(mag_epoch, bmag_mso, color='k', label='|B|')
                label_bmso = ylabel.format(
                    level=level_i, res=res_i, alt_coord='MSO')
                ax[1].set_ylabel(label_bmso)

            # Add gray dashed line for zero, legend for components,
            # and set a common y range of -20 nT to +20 nT
            for ax_i in ax:
                ax_i.axhline(0, color='gray', linestyle='--')
                ax_i.legend()
                ax_i.set_ylim([-20, 20])
    if make_plot:
        plt.show()
    print("MAG Download / Data read / Plot test complete.")


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

    parser.add_argument(
        "--make_plot",
        help="Make a plot of Bx/By/Bz (and rotations if requested).",
        action="store_true", default=False)
    parser.add_argument(
        "--test_rotation",
        help="Rotate MAG data into different coordinate systems.",
        action="store_true", default=False)

    # Download args.
    parser.add_argument(
        "--download", help="Download files to data_directory if not found.",
        action="store_true", default=False)
    parser.add_argument(
        "--remote", help="Remote to download data from (default: ssl_sprg).",
        action="store_true", default="ssl_sprg")
    parser.add_argument("--username", help="Username for data download.")
    parser.add_argument("--password", help="Password for data download.")
    args = parser.parse_args()

    # Clean the date inputs:
    start, n_days, end_date = helper.sanitize_date_inputs(
        start_date=args.start_date, n_days=1)

    # If download, set up the username / password:
    if args.download:
        if not args.username:
            six_months_ago = dt.datetime.now() - relativedelta(months=6)
            if start > six_months_ago and args.remote == "ssl_sprg":
                raise NameError(
                    "Need username to download newer than 6 months.")
            username = None
            password = None
        else:

            if args.remote == "ssl_sprg":
                username = args.username
                password = "{}_pfp".format(username)
            elif args.remote == "lasp_sdc_team":
                raise Exception(
                    "LASP SDC requires Oauth, not currently implemented.")


    # My local directory structure is based on
    # the SSL SPRG dataset, but yours might
    # be based on the LASP SDC, or might be
    # an unstructured directory (in which case
    # set local_source_tree to '')
    local_source_tree = 'ssl_sprg'

    # Contains the auth tokens if needed:
    data_download_info =\
        {'source': args.remote,
         'username': args.username,
         'password': args.password}

    # Get spice kernels
    k = spice.load_kernels(
        args.data_directory,
        start_date=start, n_days=n_days,
        download_if_not_available=args.download,
        verbose=False)
    print("Loaded kernels.")

    # Retrieve orbit ephemeris for this time period
    eph = anc.read_orbit_ephemeris(
        args.data_directory,
        start_date=start, n_days=n_days,
        download_if_not_available=args.download)

    # MAG L2
    # L2 files are available for multiple coordinate
    # systems and don't necessarily require rotation.
    b_levels = ('l2',)
    b_res = ('1sec', '30sec')
    b_coords = ('pl', 'ss', 'pc')
    b_exts = ('sav', 'sts')

    # MAG L1
    # Generally L1 files are only available in
    # payload coordinates. Any cadence other than 'full'
    # ('30sec' and '1sec') is only available for IDL sav files.
    # Usually the sav files are the quickest to read
    # for the full resolution.
    b_levels = ('l1',)
    b_res = ('30sec', '1sec', 'full')
    # b_res = ('full',)
    b_coords = ('pl',)
    b_exts = ('sav',)
    # b_exts = ('sav', 'sts', 'cdf')

    check_mag(b_levels, b_res, b_coords, b_exts,
              args.data_directory, start, n_days,
              download_info_dict=data_download_info,
              local_source_tree=local_source_tree,
              download_if_not_found=args.download,
              make_plot=args.make_plot,
              testpayloadtoMSOrotation=args.test_rotation)
