import argparse
import time
import datetime as dt

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from dateutil.relativedelta import relativedelta

from mavenpy import ngims, retrieve, file_path, helper


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

    # Clean the date inputs:
    data_dir = args.data_directory
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

        print("Checking for file updates / download if not found...")
        init_time = time.time()
        retrieve.sdc_retrieve(
            'ngims', destination_dir=data_dir,
            start_date=start, n_days=n_days, verbose=True,
            mirrored_local_source_tree='ssl_sprg',
            source=args.remote, username=username, password=password,
            dataset_name="csn-abund-(.*)", ext='csv', level='l2')
        print("Done, took {} s.".format(time.time() - init_time))

    # Get the local files:
    ngi_files = file_path.local_file_names(
        data_dir, 'ngims', start_date=start, n_days=n_days,
        source=args.local_source_tree,
        dataset_name='csn-abund-(.*)', ext='csv',
        level='l2')
    print("File names: ", ngi_files)

    # Retrieve the data:
    print("Loading the data...")
    init_time = time.time()
    ngims = [ngims.read(i) for i in ngi_files]
    print("Done, took {} s.".format(time.time() - init_time))
    # input()

    orbnum = [dat["orbit"][0] for dat in ngims]
    print(orbnum)
    # input()

    # Make a discrete colorbar
    N_orbits = orbnum[-1] - orbnum[0]
    cmaplist = []
    colors_plot = {}

    for i in range(N_orbits + 1):
        orbnum_i = orbnum[0] + i
        print(orbnum_i)

        if orbnum_i in orbnum:
            color_i = plt.cm.jet(i/(N_orbits + 1))
        else:
            # gray out
            color_i = (.5, .5, .5, 1.0)

        cmaplist.append(color_i)
        colors_plot[orbnum_i] = color_i

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'jet', cmaplist, N=N_orbits)
    print(cmaplist)



    # input()

    species = "Ar"

    fig, ax = plt.subplots()

    for dat in ngims:

        epoch_i = dat["epoch_utc"]
        alt_i = dat["alt"]
        sp_i = dat["species"]
        qual = dat["quality"]
        abundance_i = dat["abundance"]
        orb = dat["orbit"][0]
        print(orb)
        # input()
        ar_index = np.where(sp_i == species)[0]

        # Verified / unverified
        ver_index = np.where((qual == "IV") | (qual == "OV"))[0]
        unver_index = np.where((qual != "IV") & (qual != "OV"))[0]

        # Only inbound
        ver_index = np.where(qual == "IV")[0]
        unver_index = np.where(qual != "IV")[0]

        # plt.scatter(epoch_i[ar_index].astype("datetime64[ns]"), abundance_i[ar_index], marker='.')
        # plt.plot(epoch_i[ar_index].astype("datetime64[ns]"), abundance_i[ar_index])
        # plt.scatter(abundance_i[ar_index], alt_i[ar_index], marker='.')
        color_i = colors_plot[orb]

        ax.scatter(abundance_i[np.intersect1d(unver_index, ar_index)],
                   alt_i[np.intersect1d(unver_index, ar_index)], marker='x',
                   color=color_i)
        ax.scatter(abundance_i[np.intersect1d(ver_index, ar_index)],
                   alt_i[np.intersect1d(ver_index, ar_index)], marker='.',
                   color=color_i, label="{}; IV or OV".format(orb))

    ax.set_xscale('log')
    # plt.yscale('log')
    ax.legend()
    ax.set_ylabel("Altitude, km")
    ax.set_xlabel("{} density, cm-3".format(species))

    plt.subplots_adjust(right=0.85)

    ax2 = fig.add_axes([0.87, 0.1, 0.03, 0.8])
    bounds = np.linspace(orbnum[0], orbnum[-1], N_orbits + 1)
    print(bounds)
    norm = mpl.colors.BoundaryNorm(bounds - 0.5, cmap.N)

    cb = mpl.colorbar.ColorbarBase(
        ax2, cmap=cmap, ticks=bounds, boundaries=bounds - 0.5, norm=norm,
        spacing='proportional', format='%1i')

    plt.show()
