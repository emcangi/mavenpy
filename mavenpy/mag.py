import os
import datetime as dt

import numpy as np

from .helper import UNX_to_UTC, UTC_to_UNX, process_data_dict
from .read import read_cdf, read_sav

###########################
#   MAG routines          #
###########################

b_units = {"epoch": "utc", "time": "unx", "Bx": "nT", "By": "nT", "Bz": "nT",
           "Brange": "nT",
           "sc_position_x": "km", "sc_position_y": "km", "sc_position_z": "km"}

# MAG data files are available in multiple formats
# depending on processing level:
# - ASCII STS (Standard Time Series)
#   - Level 1 & 2
#   Saved in the standard l1/l2 directory tree format:
#   e.g.  maven/data/sci/mag/l2/2015/01/,
#    mvn_mag_l2_2014265pl_20140922_v01_r02.sts
#   Note: L1 until 2016/01/22 has a "d"
#   between the four-digit year and three-digit DOY.
# - IDL SAV
#   - Level 1 & 2
#   - Located in two(!) places for Level 1,
# - CDF
#   - Level 1 only
# BE AWARE! Level 1 is uncorrected for various spacecraft shenanigans
# and should not be used for final publication.

# if reading L2 files:
# - 30sec cadence: sav from /l2/sav/30sec (PL coordinates only)
# - 1sec cadence: sts from /l2 (PL, PC, SS coordinates)
#       > l2/sav/1sec (PL, PC coordinates)
# - full cadence: sts from /l2 (PL, PC, SS coordinates)
#       > l2/sav/full (PL, PC coordinates)
# No L2 files? Cautiously look for L1 files next.
# if reading L1 files (all only in payload coords):
# - 30sec cadence: sav from /l1/sav/30sec
# - 1sec cadence: sav from /l1/sav/1sec
# - full cadence: L1_sav (PL only) > L1_CDF  > L1 (sts) >>> L1/sav/full


def read(file_path, level="", ext="", coord="", res="",
         include_sc_position=True, include_unit=True,
         verbose=False):

    '''Loads daily files of MAVEN MAG data.

    This routine is based on the MAG SIS document and the MAVEN SPEDAS
    routine 'mvn_mag_load.pro'. We thank and credit the authors of both.

    Unlike other MAVEN instruments, MAGs preferred data format is an ASCII
    text file postpended with '.sts'. Level 1 data files are also
    available as CDF and IDL sav files, while Level 2 data files are also
    available as IDL sav files. Level 1 data files are also only in the
    original cadence of the instrument.

    level: string, '1' or '2'

    '''

    # Retrieve ext from the filename if not provided:
    if not ext:
        filename = os.path.split(file_path)[1]
        ext = filename[-3:]
        level = filename.split("_")[2]

    if verbose:
        print("Reading MAG level {level} data with "
              "file extension {ext} that has coordinate"
              " system {coord} and resolution {res}...".format(
                level=level, ext=ext, coord=coord, res=res))

    # read the data with the appropriate read/parse
    # based on extension:
    if ext == "cdf":
        b_struc = parse_b_cdf(read_b_cdf(file_path))
    elif ext == "sts":
        b_struc = parse_sts(read_sts(file_path))
    elif ext == "sav":
        if level == "l2":
            if coord == 'pl':
                names_i = ['time', 'ob_bpl_x', 'ob_bpl_y', 'ob_bpl_z']
            else:
                names_i = ['time', 'ob_b_x', 'ob_b_y', 'ob_b_z']
            if include_sc_position:
                names_i.extend(['posn_x', 'posn_y', 'posn_z'])
        elif level == "ql" or (level == "l1" and res == "full"):
            names_i = ['time', 'ob_bpl']
        elif level == "l1":
            names_i = ['time', 'vec']

        b_sav_struct = read_sav(
            file_path, field_names=names_i, struct_name="data")
        if verbose:
            print("sav keys: ", b_sav_struct.keys())

        b_struc = parse_sav(b_sav_struct)
    else:
        raise IOError("No MAG file reader for file with extension '{ext}',"
                      " use 'cdf', 'sts', or 'sav' instead!")

    # Calculate posix time if not included:
    if "time" not in b_struc:
        if verbose:
            print("Posix time not included, calculating...")
        posix_time = UTC_to_UNX(b_struc["epoch"])
        if verbose:
            print("Done.")

        b_struc["time"] = posix_time

    # Include unit in assignment if requested:
    if include_unit:
        b_struc = process_data_dict(b_struc, units=b_units)

    return b_struc


def read_b_cdf(bfile):
    '''Reads the MAG data from a cdf file'''
    # The only fields in the CDF file are:
    # 'YEAR', 'DOY', 'HOUR', 'MIN', 'SEC', 'MSEC',
    # 'DDAY' (fractional day), 'X' (Bx), 'Y' (By), 'Z' (Bz),
    # 'RANGE'.
    field_names = ['year', 'doy', 'hour', 'min', 'sec',
                   'msec', 'x', 'y', 'z']

    return read_cdf(bfile, field_names, lib='cdflib')


def parse_b_cdf(b):
    '''Converts the MAG data from read_b_cdf into
    a dictionary that contains the UTC time'''

    time_utc = doy_to_utc(
        b['year'], b['doy'], b['hour'], b['min'], b['sec'], b['msec'])

    bx = b['x']
    by = b['y']
    bz = b['z']

    b_dict = {"epoch": time_utc, "Bx": bx, "By": by, "Bz": bz}

    return b_dict


def read_sts(filename, method='loadtxt'):

    '''Reads the MAG data in the sts plaintext files
    saved for level 1 ('ql') and level 2 ('l2'). Outputs
    a numpy array containing the same number of columns.

    method: string, optional, 'loadtxt' or 'genfromtxt'
       - defaults to loadtxt, which only loops once over
       the data, while genfromtxt does an additional
       formatting pass (useful for string data, which
       is not present in MAG)
    '''

    # Initially tried to skip fixed # of lines, but
    # the number of header lines differs based on
    # number of Spice kernels loaded to generate the data,
    # which is contained in the header.

    # This routine will open the file and seek to the beginning
    # of the data, as defined by the end of the "object" and
    # "end_object" tags in the header.
    with open(filename, "r") as fh:

        # Loop over lines in the header until reach the
        # data. This is done by counting the number of
        # "object" definitions and waiting until all
        # object definitions are closed with "end_object"
        line_i = fh.readline().lower()
        n_obj = line_i.startswith('object')
        # print(line_i, n_obj)

        while n_obj != 0:
            line_i = fh.readline().lower()
            if line_i.startswith('object'):
                n_obj += 1
            elif line_i.startswith('end_object'):
                n_obj -= 1

        # print(line_i, n_obj)
        # print(fh.readline().lower())

        # Now we're reading data. Pull them accordingly:
        data_lines = fh.readlines()
        # print(len(remaining_lines))
        # print(remaining_lines[0])
        # print(remaining_lines[-1])
        # input()

        if method == "loadtxt":
            b = np.loadtxt(data_lines, dtype='f8')
        elif method == "genfromtxt":
            b = np.genfromtxt(
                data_lines,
                delimiter=(6, 4, 3, 3, 3, 4, 14, 11, 10, 10, 4),
                dtype='f8')

    return b


def parse_sts(b, include_position=None):
    '''Converts the numpy array produced from read_sts
    into a dictionary structure, which also contains
    the UTC time as convrted by doy_to_utc.'''

    # L1 has 11 columns:
    # Year | DOY | Hour | Min | Sec | Msec | Decimal day |
    # (Outboard magnetic field columns)
    # OB_B_X (nT)| OB_B_Y (nT)| OB_B_Z (nT)| OB_B_RANGE (nT)|

    # L2 has 18 columns:
    # Year | DOY | Hour | Min | Sec | Msec | Decimal day |
    # (Outboard magnetic field columns)
    # OB_B_X (nT)| OB_B_Y (nT)| OB_B_Z (nT)| OB_B_RANGE (nT)|
    # (S/c position columns)
    # POSN_X (km) | POSN_Y (km) | POSN_Z (km)
    # (Outboard dynamic corrections)
    # OB (nT)| OB_BD_Y (nT)| OB_BD_Z (nT)| OB_BD_RANGE (nT)|

    time_utc = doy_to_utc(b[:, 0], b[:, 1], b[:, 2], b[:, 3], b[:, 4], b[:, 5])

    bx = b[:, 7]
    by = b[:, 8]
    bz = b[:, 9]

    b_dict = {"epoch": time_utc, "Bx": bx, "By": by, "Bz": bz}

    if include_position:
        b_dict["sc_position_x"] = b[:, 10]
        b_dict["sc_position_y"] = b[:, 11]
        b_dict["sc_position_z"] = b[:, 12]

    return b_dict


def parse_sav(b):
    '''Converts the dictionary read by read_sav
    into the dictionary form output by the above read function:'''
    # print(b.dtype)

    # Time conversions
    if "time" in b:
        time_unx = b['time']
        b["epoch"] = UNX_to_UTC(time_unx)

    if "doy" in b:
        time_utc = doy_to_utc(
            b.pop('year'), b.pop('doy'), b.pop('hour'),
            b.pop('min'), b.pop('sec'), b.pop('msec'))
        b["epoch"] = time_utc

    # Vector component split
    if "vec" in b:
        b_vec = b['vec']

        b["Bx"] = b_vec[:, 0]
        b["By"] = b_vec[:, 1]
        b["Bz"] = b_vec[:, 2]

        del b['vec']

    elif "ob_bpl_x" in b:
        b["Bx"] = b.pop('ob_bpl_x')
        b["By"] = b.pop('ob_bpl_y')
        b["Bz"] = b.pop('ob_bpl_z')
        if "ob_bpl_range" in b:
            b["Brange"] = b.pop("ob_bpl_range")

    elif "ob_b_x" in b:
        b["Bx"] = b.pop('ob_b_x')
        b["By"] = b.pop('ob_b_y')
        b["Bz"] = b.pop('ob_b_z')

    if "posn_x" in b:
        b["sc_position_x"] = b.pop('posn_x')
        b["sc_position_y"] = b.pop('posn_y')
        b["sc_position_z"] = b.pop('posn_z')

    return b


def doy_to_utc(year, doy, hour, minute, sec, msec):

    # Sometimes these fields are
    # improperly read as integers,
    # which will cause strange additive behaviors.
    # Convert them to floats to avoid this.
    year = year.astype('float')
    doy = doy.astype('float')
    hour = hour.astype('float')
    minute = minute.astype('float')
    sec = sec.astype('float')
    msec = msec.astype('float')

    # Get the total number of seconds that have transpired
    # since the year started.
    total_sec = (doy - 1)*24*60*60 + hour*60*60 +\
        minute*60 + sec + (1e-3*msec)
    time_utc =\
        [dt.datetime(int(y), 1, 1) + dt.timedelta(seconds=s) for (y, s) in
         zip(year, total_sec)]

    time_utc = np.array(time_utc)

    # print(yyyy[0], doy[0], hh[0], MM[0], ss[0], msec[0])
    # print(yyyy[-1], doy[-1], hh[-1], MM[-1], ss[-1], msec[-1])

    # plt.figure()
    # plt.plot(total_sec)
    # plt.show()

    return time_utc
