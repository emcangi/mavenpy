import re
import datetime as dt

from dateutil.relativedelta import relativedelta

from . import helper


Mars_orbit_insertion = '2014 9 22'

# Safe mode events:
# - 2015/04/04 - 04/13
# - 2015/06/29 - 07/01
# - 2019/09/05 - 09/13
# - 2022/02/22 - 04/21
safe_mode_periods = ['2015/04/04 - 2015/04/13',
                     '2015/06/29 - 2015/07/01',
                     '2019/09/05 - 2019/09/13',
                     '2022/02/22 - 2022/04/22',
                     '2023/02/16 - 2023/02/20']


# The options for remotes that are tested to work:
remote_options = ("lasp_sdc_team", "lasp_sdc_public", "ssl_sprg")

# Identify which instruments produce a file daily
# and which produce a file each orbit,
file_per_orbit = ("ngi", "iuv", "acc")
file_per_day = ("swe", "swi", "sep", "euv", "sta", "mag", "lpw", "pfp")


# Available instruments + levels + formats:
formats = {}

# Ancillary file info
formats["anc"] = {}
formats["anc"]["datasets"] = "monthly_ephemeris"
formats["anc"]["source"] =\
    {"monthly_ephemeris": "ssl_sprg"}
formats["anc"]["ext"] = "sav"

# ACC
formats["acc"] = {}
formats["acc"]["level"] = "l3"
formats["acc"]["ext"] = "tab"
formats["acc"]["datasets"] = {"l3": "pro-acc-p(.*)"}
formats["acc"]["source"] = {"l3": "ssl_sprg"}

# EUV
formats["euv"] = {}
formats["euv"]["level"] = ("l2", "l3")
formats["euv"]["ext"] = "cdf"
formats["euv"]["datasets"] =\
    {"l2": "bands", "l3": ("daily", "minute")}
formats["euv"]["source"] =\
    {"l2": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public"),
     "l3": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public")}

# IUVS
formats["iuv"] = {}
formats["iuv"]["level"] = ("l1c", "l2")
formats["iuv"]["ext"] = "fits.gz"
formats["iuv"]["datasets"] =\
    {"l1c": ("corona", "disk", "echelle", "limb", "occultation"),
     "l2": ("corona", "disk", "limb", "occultation")}
formats["iuv"]["source"] =\
    {"l1c": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public"),
     "l2": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public")}


# LPW
formats["lpw"] = {}
formats["lpw"]["level"] = "l2", "l0b"
formats["lpw"]["ext"] = "cdf"
formats["lpw"]["datasets"] =\
    {"l2": ("lpiv", "lpnt", "mrgscpot", "wn",
            "we12", "we12bursthf", "we12burstlf",
            "wspecact", "wspecpas"),
     "l0b": ("act", "adr", "atr", "euv",
             "hsk", "hsbmlf", "msbmmf",
             "spechfact", "speclfact", "speclfpas",
             "specmfact", "specmfpas"
             "swp1", "swp2")}
formats["lpw"]["source"] =\
    {"l2": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public"),
     "l0b": ("lasp_sdc_team", "lasp_sdc_public")}


# NGIMS
formats["ngi"] = {}
formats["ngi"]["level"] = ("l1b", "l2", "l3")
formats["ngi"]["ext"] = "csv"
formats["ngi"]["datasets"] =\
    {"l1b": ("cal-hk-(.*)", "cal-mkr-(.*)", "cal-msg-(.*)",
             "osnb-(.*)", "osion-(.*)"),
     "l2": ("csn-abund-(.*)", "ion-abund-(.*)", "(csn|ion)-(.*)"),
     "l3": ("res-den-(.*)", "res-sht-(.*)", "wind-(.*)")}
formats["ngi"]["source"] = ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public")


# PFP
formats["pfp"] = {}
formats["pfp"]["level"] = "l0"
formats["pfp"]["ext"] = "dat"
formats["pfp"]["source"] = {"l0": ("ssl_sprg", "lasp_sdc_team")}
formats["pfp"]["datasets"] = {"l0": ("all", "arc", "svy")}

# ROSE AKA RSE
formats["rse"] = {}
formats["rse"]["level"] = ("l2", "l3")
formats["rse"]["ext"] = "tab"
formats["rse"]["datasets"] = {"l2": ("dlf", "fup", "sky"), "l3": "edp"}
formats["rse"]["source"] =\
    {"l2": ("lasp_sdc_team", "lasp_sdc_public"),
     "l3": ("lasp_sdc_team", "lasp_sdc_public")}

# SEP
formats["sep"] = {}
formats["sep"]["level"] = ("l1", "l2", "l3", "anc")
formats["sep"]["ext"] = {"l1": "sav", "l2": "cdf", "l3": "sav",
                         "anc": ("sav", "cdf")}
formats["sep"]["datasets"] = {"l1": ("full", "01hr", "32sec", "5min"),
                              "l2": ("s(.*)-raw-svy-full",
                                     "s(.*)-cal-svy-full",
                                     "s(.*)-svy-full"),
                              "l3": "pad",
                              "anc": ""}
formats["sep"]["source"] =\
    {"l1": "ssl_sprg",
     "l3": ("ssl_sprg", "lasp_sdc_team"),
     "l2": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public"),
     "anc": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public")}


# STATIC
formats["sta"] = {}
formats["sta"]["level"] = ("l2", "l3", "iv1", "iv2", "iv3", "iv4")
formats["sta"]["ext"] =\
    {"l2": "cdf", "l3": ("sav", "tplot"),
     "iv1": "cdf", "iv2": "cdf", "iv3": "cdf", "iv4": "cdf"}
sta_appids =\
    ("2a-hkp", "c0-64e2m", "c6-32e64m", "c8-32e16d",
     "ca-16e4d16a", "d0-32e4d16a8m", "d1-32e4d16a8m",
     "d4-4d16a2m", "d6-events", "d7-fsthkp",
     "d8-12r1e", "d9-12r64e", "da-1r64e", "db-1024tof")
sta_calib_ids = ("c0-64e2m", "c6-32e64m", "c8-32e16d",
                 "ca-16e4d16a", "d0-32e4d16a8m", "d1-32e4d16a8m",)

formats["sta"]["datasets"] =\
    {"l2": sta_appids,
     "iv1": sta_calib_ids, "iv2": sta_calib_ids,
     "iv3": sta_calib_ids, "iv4": sta_calib_ids,
     "l3": ("density", "temperature")}
formats["sta"]["source"] =\
    {"l3": "ssl_sprg",
     "l2": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public"),
     "iv1": "ssl_sprg", "iv2": "ssl_sprg",
     "iv3": "ssl_sprg", "iv4": "ssl_sprg"}

# SWEA
formats["swe"] = {}
formats["swe"]["level"] = ("l2", "l3", "anc")
formats["swe"]["ext"] = {"l2": "cdf", "l3": "sav", "anc": "sav"}
formats["swe"]["datasets"] =\
    {"l2": ("svyspec", "svypad", "arcpad", "svy3d", "arc3d"),
     "l3": ("shape", "padscore", "scpot"),
     "anc": "quality"}
formats["swe"]["source"] =\
    {"l2": ("ssl_sprg", "lasp_sdc_public", "lasp_sdc_team"),
     "l3": ("ssl_sprg", "lasp_sdc_team"),
     "anc": "ssl_sprg"}

# SWIA
formats["swi"] = {}
formats["swi"]["level"] = "l2"
formats["swi"]["ext"] = "cdf"
formats["swi"]["datasets"] =\
    {"l2": ("onboardsvymom", "onboardsvyspec", "coarsearc3d", "coarsesvy3d",
            "finearc3d", "finesvy3d")}
formats["swi"]["source"] =\
    {"l2": ("ssl_sprg", "lasp_sdc_public", "lasp_sdc_team")}

# MAG
mag_avail_res = {
    ("l1", "sts"): "full",
    ("l1", "cdf"): "full",
    ("l1", "sav"): ("full", "1sec", "30sec"),
    ("l2", "sts"): ("full", "1sec"),
    ("l2", "sav"): ("full", "1sec", "30sec")}

formats["mag"] = {}
formats["mag"]["level"] = ("l1", "l2")
formats["mag"]["ext"] = {"l1": ("sts", "cdf", "sav"),
                         "l2": ("sav", "sts")}
formats["mag"]["source"] =\
    {"l1": {"sts": ("ssl_sprg", "lasp_sdc_team"),
            "sav": "ssl_sprg", "cdf": "ssl_sprg"},
     "l2": {"sts": ("ssl_sprg", "lasp_sdc_team", "lasp_sdc_public"),
            "sav": "ssl_sprg"}}

mag_coordinate_mapping = {"mso": "ss", "geo": "pc", "pl": "pl"}
mag_avail_coords = {
    ("l2", "30sec", "sav"): ("pl"),
    ("l2", "1sec", "sts"): ("pl", "pc", "ss"),
    ("l2", "1sec", "sav"): ("pl", "pc"),
    ("l2", "full", "sts"): ("pl", "pc", "ss"),
    ("l2", "full", "sav"): ("pl", "pc"),
    ("l1", "30sec", 'sav'): ("pl",),
    ("l1", "1sec", 'sav'): ("pl",),
    ("l1", "full", "sav"): ("pl",),
    ("l1", "full", "sts"): ("pl",),
    ("l1", "full", "cdf"): ("pl",)}


# Raw pfp data name
raw_pfp_name = "mvn_pfp_{data}_l0_{{yyyy}}{{mm}}{{dd}}_v[0-9][0-9][0-9].dat"
# For SWIA, SWEA, SEP, EUV, STATIC, LPW
daily_name = "mvn_{tla}_{level}_{dataset_name}_"\
    "{{yyyy}}{{mm}}{{dd}}_v[0-9][0-9]_r[0-9][0-9].{ext}"
# For IUVS, NGIMS, ACC
hourly_name = "mvn_{tla}_{level}_{dataset_name}_"\
    "{{yyyy}}{{mm}}{{dd}}T{{hhMMSS}}_v[0-9][0-9]_r[0-9][0-9].{ext}"

# Specialty data names
mag_nonsav_name = "mvn_mag_{alt_level}_{{yyyy}}{{doy}}{coord}{alt_res}_"\
    "{{yyyy}}{{mm}}{{dd}}_v[0-9][0-9]_r[0-9][0-9].{ext}"
mag_sav_name = "mvn_mag_{level}_{coord}_{res}_{{yyyy}}{{mm}}{{dd}}.sav"
swea_swi_regid_name = "mvn_swia_regid_{{yyyy}}{{mm}}{{dd}}_"\
    "v[0-9][0-9]_r[0-9][0-9].sav"

# i.e. mvn_sep_l1_20141231_1day.sav
# versioned after 2023/01 (mvn_sep_l1_20230630_v006.sav)
# sep_l1_name = "mvn_sep_l1_{yyyy}{mm}{dd}_*.sav"
sep_l1_name = "mvn_sep_l1_{yyyy}{mm}{dd}_(.*).sav"
sta_iv_name = 'mvn_sta_l2_{dataset_name}_{{yyyy}}{{mm}}{{dd}}_{iv_num}.cdf'


def during_safemode(time):
    '''Function to check if a time takes place during
    safe mode.'''

    datetime_dt = helper.date_to_dt(time)

    safemode = False
    safe_mode_dt =\
        [[parse(j) for j in i.split(" - ")] for i in
         safe_mode_periods]

    for missing_dt in safe_mode_dt:
        # print(missing_dt)
        if isinstance(missing_dt, list):
            if missing_dt[0] <= datetime_dt <= missing_dt[1]:
                safemode = True
                break
        elif missing_dt == datetime_dt:
            safemode = True
            break

    return safemode


def remote_base_directory(source):
    '''Return the base directory as a tuple
    where MAVEN data is saved for a given source'''

    if source not in remote_options:
        raise IOError(
            "Source '{}' not available, "
            "try instead: {}".format(source, remote_options))

    if source == "lasp_sdc_team":
        root = ("maven", "sdc", "team")
    elif source == "lasp_sdc_public":
        root = ("maven", "sdc", "public")
    elif source == "ssl_sprg":
        root = ("maven",)
    # elif source == "pds":
    #     # Not handled at this time.

    return root


def path(instrument_tla, level, ext="", dataset_name="", res=""):

    '''Returns the relative path for a given instrument dataset with a given
    level as a tuple, either to generate a local path or a URL to
    access the data remotely.

    instrument_tla: string, three-letter-acronym for a given instrument
    level: string, e.g. l1 or l2, level of data processing.
    dataset_name: optional, string, name of dataset e.g. 'pad' for SWEA
    ext: optional, string, file extension of dataset
        (some datasets like SEP and MAG contain data in directories
         named by the file extension)
    res: optional, string, resolution of the dataset, only for MAG
        (e.g. '1sec' or '30sec' or 'full')

    '''

    # Common default path:
    p = ("data", "sci", instrument_tla, level, "{yyyy}", "{mm}")

    # Special cases:
    # Ancillary monthly ephemeris:
    if instrument_tla == "anc":
        p = ("anc", "spice", "sav")

    if instrument_tla == "pfp":
        if dataset_name == "all":
            p = ("data", "sci", "pfp", "l0_all", "{yyyy}", "{mm}")
        else:
            p = ("data", "sci", "pfp", "l0")

    if instrument_tla == "iuv":
        p = ("data", "sci", "iuv", level, dataset_name, "{yyyy}", "{mm}")

    # Non-STS MAG data (variable file tree, prefer l1_sav over sav/l1 for
    # full resolution bc faster load.
    if instrument_tla == "mag" and ext != "sts":
        if (ext == "sav" and level == "l1") and res == "full":
            p = ("data", "sci", "mag", "l1_sav", "{yyyy}", "{mm}")
        elif ext == "sav":
            p = ("data", "sci", "mag", level, "sav", res, "{yyyy}", "{mm}")
        elif (ext == "cdf" and level == "l1"):
            p = ("data", "sci", "mag", "l1_cdf", "{yyyy}", "{mm}")

    # SWEA Level 3 has separate directories for each dataset_name
    if instrument_tla == "swe" and level == "l3":
        p = ("data", "sci", "swe", level, dataset_name, "{yyyy}", "{mm}")

    # Similar for SEP Level 3, with also a subdir for extension
    if instrument_tla == "sep":
        if level == "l3":
            p = ("data", "sci", "sep", level, dataset_name,
                 ext, "{yyyy}", "{mm}")
        elif level == "anc":
            p = ("data", "sci", "sep", level, ext, "{yyyy}", "{mm}")
        elif level == "l1":
            if dataset_name == "full":
                p = ("data", "sci", "sep", level, ext, "{yyyy}", "{mm}")
            else:
                p = ("data", "sci", "sep", level, "{}_{}".format(
                    ext, dataset_name), "{yyyy}", "{mm}")

    return p


def filename(instrument_tla, level="2", dataset_name=None, ext=None,
             coord=None, res=None):
    """Return the name of a MAVEN instrument datafile.

    instrument_tla: string, name of requested instrument data,
       e.g. 'sep' or 'swi'
    level: string, data product level, e.g. 'l2'
    dataset_name: string, name of requested dataset
    coord: string, coord system of requested filename, reqd
        for MAG and monthly ephemeris
    res: string, interval of requested filename, reqd for MAG"""

    if instrument_tla == "mag":
        # Need ext, res, and coord to make filename
        if not ((ext and res) and coord):
            raise IOError(
                "Need to define a res, coordinate system, and ext"
                " to construct filename.\n"
                "(run specification.check_if_dataset_exist to confirm"
                "exist)")
        # The MAG sav files have a different name format.
        # If getting the L1 sav file, there are two options,
        # mag/l1/sav/full (which have the mag_sav_name format)
        # or mag/l1_sav (has the mag_nonsave_name format).
        # The l1_sav files load faster with the scipy.io
        # readsav reader, so volunteer the nonsav format
        if ext == 'sav' and not (res == "full" and level == 'l1'):
            data_name = mag_sav_name.format(
                level=level, coord=coord, res=res)
        else:
            alt_level = ('ql' if level == 'l1' else level)
            alt_res = ("1s" if res == '1sec' else '')
            data_name = mag_nonsav_name.format(
                alt_level=alt_level, alt_res=alt_res, ext=ext, coord=coord)

    elif instrument_tla == "swe" and dataset_name == "swia_regid":
        data_name = swea_swi_regid_name
    elif instrument_tla == 'pfp':
        data_name = raw_pfp_name.format(data=dataset_name)
    elif instrument_tla == "sep" and level == 'l1':
        data_name = sep_l1_name
    elif instrument_tla == "sta" and "iv" in level:
        data_name = sta_iv_name.format(
            dataset_name=dataset_name, iv_num=level)
        print(data_name)
        input()

    elif instrument_tla in file_per_orbit:
        data_name = hourly_name.format(
            tla=instrument_tla, level=level,
            dataset_name=dataset_name, ext=ext)
    elif instrument_tla in file_per_day:
        if not ext and level == "l2":
            ext = "cdf"
        data_name = daily_name.format(
            tla=instrument_tla, level=level,
            dataset_name=dataset_name, ext=ext)

    return data_name


def check_if_dataset_on_remote(instrument, source, auth_present=False,
                               level='', ext='', dataset='',
                               start_utc=''):

    '''Inspects if a instrument dataset is available on a given
    remote, and if not, provides an error with the remotes it
    is available on.

    instrument: string, name of instrument
    source: string, name of remote to be downloaded from, e.g. 'ssl_sprg'
    auth_present: Boolean, True if username/password provided
    start_utc: Date time of the beginning of the data to download,
        required if looking up data on SSL_SPRG
    '''

    # First, make sure have credentials if accessing ssl sprg or
    # lasp sdc team:
    if not auth_present:

        if source == "lasp_sdc_team":
            raise TypeError(
                "Unable to access LASP SDC team website"
                " at this time since requires OAuth.")

        if source == "ssl_sprg":
            six_months_ago = dt.datetime.now() - relativedelta(months=6)
            if start_utc > six_months_ago:
                raise TypeError(
                    "Require authentication information (username and password)"
                    " for data on ssl_sprg within 6 months.")

    # Get the TLA for the instrument:
    tla = instrument[:3].lower()
    hosts = formats[tla]["source"]
    # print(source_avail)

    # If supplying a dataset name, that is the dict key for
    # the source
    if tla == "anc":

        if not dataset:
            raise IOError("Need a dataset name if requesting"
                          " ancillary info! Available datasets: {}".format(
                            formats["anc"]["datasets"]))

        level = dataset

    # Need an extension and level to get MAG data, since some datasets
    # are only available from specific hosts
    if tla == "mag" and (not ext or not level):
        raise IOError("Provide extension and level to verify"
                      " if MAG data is retrievable from {}.".format(source))

    # Some of the data is available
    if isinstance(hosts, dict):
        if level not in hosts:
            raise IOError(
                "No {} {} files available, select another level: ".format(
                    instrument, level, hosts.keys()))

        source_avail = hosts[level]
        if tla == "mag":
            source_avail = source_avail[ext]

        # print(source_avail, dataset, ext)
        # input()

        check_dataset_remote_exist_str_tuple(
            source, source_avail, instrument, level)

    else:
        check_dataset_remote_exist_str_tuple(
            source, hosts, instrument, level)


def check_dataset_remote_exist_str_tuple(source, source_available,
                                         instrument, level):

    if isinstance(source_available, str):
        if source not in source_available:
            raise IOError(
                "No {} {} files at {}, available instead at {}".format(
                    instrument, level, source, source_available))

    if isinstance(source_available, tuple):
        check_if_in = [1 if source in i else 0 for i in source_available]
        if sum(check_if_in) == 0:
            raise IOError(
                "No {} {} files at {}, available instead at {}".format(
                    instrument, level, source, source_available))


def missing_dataset_error_message(name, level, ext, datasets, dataset=''):

    if not dataset:
        return "Must provide a dataset name! Available {name} {level} {ext}"\
            " files for dataset(s): {ds}".format(
                name=name, level=level, ext=ext, ds=datasets)
    else:

        return "No available {name} {level} {ext} files for"\
            " dataset '{ds_i}', use '{ds} instead'".format(
                name=name, level=level, ext=ext, ds_i=dataset,
                ds=datasets)


def check_if_dataset_exist(instrument, ext="", level="", dataset="",
                           coord="", res=""):

    tla = instrument[:3].lower()
    NAME = instrument.upper()

    # Grab info for given instrument
    info = formats[tla]

    # For instrument types with levels, see if requested level
    # exists.
    if "level" in info:

        # Grab the available levels:
        levels = formats[tla]["level"]

        # If level is an empty string, raise error since level must be defined
        if not level:
            raise IOError(
                "Must provide a level! Available {NAME} files for level(s): "
                "{levels}".format(NAME=NAME, levels=levels))

        if level not in levels:
            raise IOError(
                "No available {NAME} {level} files, use level"
                " {levels} instead".format(
                    NAME=NAME, level=level, levels=levels))

    # Check if dataset column exists
    # Note that PFP and MAG lack data subsets bc...
    # they're the whole shebang from that instrument.
    # no subsets exist.
    if "datasets" in info:

        # Get available datasets
        ds = info["datasets"]
        # print(ds)

        # If dataset is an empty string, raise error since dataset must
        # be defined
        if not dataset:
            if isinstance(ds, dict):
                ds_i = ds[level]
            else:
                ds_i = ds
            ds_msg = missing_dataset_error_message(
                NAME, level, ext, ds_i, dataset="")
            raise IOError(ds_msg)

        # If dataset is defined and the known dataseet is
        # a string, do a search:
        if isinstance(ds, str):
            ds_msg = missing_dataset_error_message(
                NAME, level, ext, ds, dataset=dataset)
            # if ds != dataset:

            if not re.match(re.compile(ds), dataset):
                raise IOError(ds_msg)

        if isinstance(ds, dict):
            available_ds_for_level = ds[level]
            ds_msg = missing_dataset_error_message(
                NAME, level, ext, available_ds_for_level, dataset=dataset)
            if isinstance(available_ds_for_level, str):
                if not re.match(re.compile(available_ds_for_level), dataset) and dataset != available_ds_for_level:
                    raise IOError(ds_msg)
            elif isinstance(available_ds_for_level, tuple):
                check_if_in =\
                    [1 if re.match(re.compile(avail_ds), dataset) else 0
                     for avail_ds in available_ds_for_level]
                # print(check_if_in)
                # print(available_ds_for_level)
                # print(dataset)
                # print(dataset in available_ds_for_level)
                # input()

                if sum(check_if_in) == 0 and dataset not in available_ds_for_level:
                    raise IOError(ds_msg)

    # All instruments have an extension
    exts = formats[tla]["ext"]

    # If extension is an empty string, raise error since ext must be defined
    if not ext:

        if isinstance(exts, dict):
            ext_i = exts[level]
        else:
            ext_i = exts

        raise IOError(
            "Must provide a file extension! Available {NAME} files for ext(s):"
            " {exts}".format(NAME=NAME, exts=ext_i))

    # If the extension isn't the requested one:
    if exts != ext:
        if isinstance(exts, str):
            raise IOError(
                "No available {NAME} {level} files with format"
                " '{format}', use  '{formats} instead'".format(
                    NAME=NAME, level=level,
                    format=ext, formats=exts))

        if isinstance(exts, tuple):
            if ext not in exts:
                raise IOError(
                    "No available {NAME} {level} files with format"
                    " '{format}', use '{formats} instead'".format(
                        NAME=NAME, level=level,
                        format=ext, formats=exts))
        if isinstance(exts, dict):
            available_exts_per_level = exts[level]
            if ext not in available_exts_per_level:
                raise IOError(
                    "No available {NAME} {level} files with format"
                    " '{format}', use '{formats} instead'".format(
                        NAME=NAME, level=level,
                        format=ext, formats=exts))

    if tla == "mag":
        # print("MAG check")
        # Check if cadence exists for MAG level and extension
        available_res = mag_avail_res[(level, ext)]
        # print(available_res, res)

        if not res:
            raise IOError(
                "Require resolution to load MAG {level} {ext} files."
                " Options: {avail_res}".format(
                    level=level, ext=ext, avail_res=available_res))

        if res not in available_res:
            raise IOError(
                "No available MAG {level} {ext} files for "
                "resolution '{res}',"
                " use resolution {avail_res} instead".format(
                    level=level, ext=ext, res=res,
                    avail_res=available_res))

        # See if coordinate system needs to be mapped from common names
        # (e.g. 'mso' or 'geo')
        available_coords =\
            mag_avail_coords[(level, res, ext)]
        coord = coord.lower()
        if coord in mag_coordinate_mapping:
            coord = mag_coordinate_mapping[coord]

        # print(available_coords, coord)

        if not coord:
            raise IOError(
                "Choose coordinate system to load MAG {level} {ext} files."
                " Options: {coords}".format(
                    level=level, ext=ext, res=res,
                    coords=available_coords))

        if coord not in available_coords:
            raise IOError(
                "No MAG {level} {ext} file w/ res {res} "
                "for '{coord}' coordinate system, use instead: "
                "{available_coords}".format(
                    level=level, ext=ext, res=res,
                    coord=coord, available_coords=available_coords))
