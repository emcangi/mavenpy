import argparse
import itertools
import datetime as dt

from dateutil.relativedelta import relativedelta

from mavenpy import file_path, retrieve, helper, specification, spice


# Default datasets to retrieve
default_datasets = {
    ("euv", "l3"): "minute",
    ("iuv", "l2"): "limb",
    ("lpw", "l2"): "lpnt",
    ("ngi", "l2"): "csn-abund-(.*)",
    ("pfp", "l0"): "all",
    ("sep", "l2"): 's(.*)-cal-svy-full',
    ("sep", "l1"): "01hr",
    ("sta", "l2"): ("c0-64e2m", "c6-32e64m"),
    ("swe", "l2"): ("svyspec", "svypad"),
    ("swi", "l2"): ("onboardsvymom", "coarsesvy3d")}


default_mag = {"ext": "sav", "coord": "pl", "res": "30sec"}

if __name__ == "__main__":

    # Arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data_directory", help="Directory containing MAVEN data.")
    parser.add_argument(
        "--start_date",
        help="Start date of download range (YYYY-MM-DD).",
        required=True,
        type=lambda s: dt.datetime.strptime(s, "%Y-%m-%d"))
    parser.add_argument(
        "--end_date",
        help="End date of download range (YYYY-MM-DD).",
        type=lambda s: dt.datetime.strptime(s, "%Y-%m-%d"))
    parser.add_argument(
        "--n_days", help="Number of days since start date to recover date.",
        type=int)

    # Instrument keywords:
    parser.add_argument(
        "--instrument", help="Instrument name (e.g. SEP).", nargs='*')
    parser.add_argument(
        "--level",
        help="Instrument data level (e.g. 'l2' for level 2), "
             "will assume level 2 if not supplied.", nargs='*', default=['l2'])

    # Keyword to get Spice info in range:
    parser.add_argument(
        "--spice", help="Download spice files.", action='store_true')

    # Optional keywords for instrument download:
    parser.add_argument(
        "--dataset_name",
        help="Instrument data set name (e.g. 'csn-abund-*' for NGIMS L2, "
             "or 'c6-32e64m' for STATIC L2).",  nargs='*')
    # Required keywords for MAG:
    parser.add_argument(
        "--mag_res",
        help="Resolution of MAG file downloaded (30sec, 1sec, full).")
    parser.add_argument(
        "--mag_ext", help="Extension of MAG file downloaded (sts or sav).")
    parser.add_argument(
        "--mag_coord",
        help="Coordinate system of MAG file downloaded "
             "(GEO = pc, MSO = ss, PL = pl).")
    # Required keywords for IUVS
    parser.add_argument(
        "--orbit_segment", default='periapse',
        help="Orbit segment of to-be-downloaded IUVS observation.")
    parser.add_argument(
        "--imaging_mode", default='',
        help="Imaging mode of to-be-downloaded IUVS observation.")
    # Keywords for Spice access
    parser.add_argument(
        "--sc_position",
        help="Get Spice kernels for determining spacecraft position,"
             " defaults to True.",
        action='store_true',
        default=True)
    parser.add_argument(
        "--sc_pointing",
        help="Get Spice kernels for determining spacecraft pointing "
             "(can be skipped if not intending to do vector transforms)."
             "WARNING: there are a LOT of s/c pointing files, "
             "only set this if you need 'em because they will need "
             "a lot of download time.",
        action='store_true')
    parser.add_argument(
        "--app_pointing",
        help="Get Spice kernels for determining APP pointing "
             "(can be skipped if not intending to do vector transforms)."
             "WARNING: there are a LOT of APP pointing files, "
             "only set this if you need 'em because they will need "
             "a lot of download time.",
        action='store_true')
    parser.add_argument(
        "--prompt_for_spice_download",
        help="Raises a prompt for downloading Spice kernels, if desired.",
        action='store_true')


    # Keyword to include print statements
    parser.add_argument(
        "--verbose", help="Enable print statements to help debugging.",
        action='store_true')

    # List all available files in range after download:
    parser.add_argument(
        "--list", help="Print files for instruments.",
        action='store_true')

    # Download keywords:
    parser.add_argument(
        "--remote", help="Remote to download data from (default: ssl_sprg).",
        default="ssl_sprg")

    parser.add_argument(
        "--no_mirror",
        help="Boolean to disable mirroring of a remote directory tree."
             " If this is included in the commandline, will save all "
             "files to folder regardless of structure.",
        action="store_true")

    parser.add_argument(
        "--local_source_tree",
        help="Name of remote to mirror or the way the"
             " local tree is already mirrored. Assumes ssl_sprg by default,"
             " can be lasp_sdc_public, lasp_sdc_team, or empty.",
        default="ssl_sprg")
    parser.add_argument("--username", help="Username for data download.")
    parser.add_argument("--password", help="Password for data download.")

    args = parser.parse_args()

    # Check if instrument provided or spice selected:
    if not args.instrument and not args.spice:
        raise IOError("Provide either a MAVEN dataset "
                      "(--instrument <<INSTRUMENT_NAME>>)"
                      " or set spice argument (--spice).")

    # Get time range (checks end_date or n_days)
    start_utc, n_days, end_utc = helper.sanitize_date_inputs(
        start_date=args.start_date, end_date=args.end_date,
        n_days=args.n_days)

    # If download, set up the username / password:
    username = ""
    password = ""
    if args.username and args.remote == "ssl_sprg":
        username = args.username
        password = "{}_pfp".format(username)

    remote_info =\
        {"source": args.remote, "username": username,
         "password": password,
         "mirror_remote_tree": not args.no_mirror,
         "mirrored_local_source_tree": args.local_source_tree}

    # If data directory not provided, pull from the IDL cshrc if not defined
    if not args.data_directory:
        data_directory = file_path.get_IDL_data_dir()
    else:
        data_directory = args.data_directory

    if args.spice:
        mirror_spedas_dir_tree =\
            (args.local_source_tree == "ssl_sprg" and not args.no_mirror)

        k = spice.load_kernels(
            data_directory,
            start_date=start_utc, end_date=end_utc,
            n_days=n_days,
            download_if_not_available=True,
            load_spacecraft=args.sc_position,
            load_spacecraft_pointing=args.sc_pointing,
            load_APP=args.app_pointing,
            mirror_spedas_dir_tree=mirror_spedas_dir_tree,
            verbose=args.verbose,
            prompt_for_download=args.prompt_for_spice_download)

    # Check if instrument available:
    if args.instrument:
        # Args will supply a list
        instrument = args.instrument

        for instrument_i, level_i in itertools.product(instrument, args.level):

            # Get the lowercase tla
            tla_i = instrument_i.lower()[:3]

            # Get available formats:
            format_i = specification.formats[tla_i]
            levels_i = format_i["level"]

            # Is the level available for this?
            if level_i not in levels_i:
                raise IOError(
                    "Level '{}' not available for {}, "
                    "available options: {}".format(
                        level_i, instrument_i, levels_i))

            # Get the datasets available:
            dataset_params = {"level": level_i}
            if tla_i == "mag":
                # Update the default MAG parameters if provided
                res, ext, coord = args.mag_res, args.mag_ext, args.mag_coord
                if (res and ext) and coord:
                    default_mag = {"ext": ext, "coord": coord, "res": res}
                dataset_params.update(default_mag)
            else:
                if tla_i == "iuv":
                    default_iuv = {"orbit_segment": args.orbit_segment,
                                   "imaging_mode": args.imaging_mode}
                    dataset_params.update(default_iuv)
                available_datasets = format_i["datasets"][level_i]
                available_ext = format_i["ext"]

                if isinstance(available_ext, str):
                    ext_i = available_ext
                else:
                    ext_i = available_ext[level_i]
                dataset_params["ext"] = ext_i

                if args.dataset_name:
                    for d_name in args.dataset_name:
                        if d_name not in available_datasets:
                            raise IOError(
                                "Requested dataset not in available"
                                " datasets: {}".format(available_datasets))
                    dataset_i = args.dataset_name
                else:
                    # If only one available dataset, download that one:
                    if isinstance(available_datasets, str):
                        dataset_i = available_datasets
                    elif (tla_i, level_i) in default_datasets:
                        dataset_i = default_datasets[(tla_i, level_i)]

                dataset_params["dataset_name"] = dataset_i

            # print(ext_i, dataset_i)
            # input()

            retrieve.sdc_retrieve(
                tla_i, **remote_info, **dataset_params,
                destination_dir=data_directory,
                start_date=start_utc, end_date=end_utc, verbose=args.verbose)
            # input()
            print("{} files updated.".format(tla_i.upper()))

            if args.list:
                fnames = file_path.local_file_names(
                    data_directory, tla_i,
                    start_date=start_utc, end_date=end_utc,
                    **dataset_params)
