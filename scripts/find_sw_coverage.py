import argparse
import datetime as dt
import numpy as np
import os

from dateutil import rrule
from dateutil.parser import parse

from mavenpy import helper, specification
from mavenpy.file_path import local_file_names
from mavenpy.sep import read_anci
from mavenpy import spice, anc
from mavenpy.mars_shape_conics import solar_wind_indices

n_carrington_days = 27.2753


def retrieve_sw_indices(data_directory, date, data_source):

    if data_source == "sep_anci":
        # Retrieve a SEP ancillary file for a given day
        sep_anci_filename = local_file_names(
            data_directory, start_date=date, n_days=1, source='ssl_sprg')
        sep_anci = read_anci(sep_anci_filename)

        sc_time_unx = sep_anci["time"][0]
        sc_mso = sep_anci["maven_mso"][0]
        x, y, z = sc_mso[:, 0], sc_mso[:, 1], sc_mso[:, 2]
    elif data_source == "spice":
        # Generate a time list for each minute of a given day
        t_utc, sc_time_unx, x, y, z = spice.load_MAVEN_position(
            date, n_days=1, n_sample_points=(24*60))

    # Given the spacecraft position, retrieve the indices
    sw_indices, sw_alt, not_sw_indices = solar_wind_indices(x, y, z)
    n_sw = len(sw_indices)

    return n_sw, sw_indices, sc_time_unx[sw_indices]


def search_sw_coverage(data_directory, start_date_dt, end_date_dt,
                       data_source):

    '''Apply binary search on either daily ephemeris files or Spice
    information to find when solar wind coverage ends.'''

    # Get the # of observations on the first and last day of the
    # month:
    Nsw_start_day = retrieve_sw_indices(
        data_directory, start_date_dt, data_source)[0]
    Nsw_end_day = retrieve_sw_indices(
        data_directory, end_date_dt, data_source)[0]

    n_obs_in_sw = [Nsw_start_day, Nsw_end_day]

    # If there's no SW coverage or SW coverage at start and end of month,
    # return [Nsw_day1, Nsw_dayN] (observations of SW),
    # None (Time enter SW), None (Time exit SW)
    if sum(n_obs_in_sw) == 0 or (Nsw_start_day > 0 and Nsw_end_day > 0):
        return n_obs_in_sw, None, None

    u = end_date_dt
    l = start_date_dt

    time_enter_sw_utc = None
    time_leave_sw_utc = None

    while not time_enter_sw_utc and not time_leave_sw_utc:
        mid_day = (u - l) / 2 + l
        n_sw_i, sw_indices, sw_time_unx = retrieve_sw_indices(
            data_directory, mid_day, data_source
        )
        # print(mid_day)

        # Crossing into SW
        if Nsw_start_day == 0:
            if n_sw_i == 0:
                # Has crossed into sw, but not yet - set upper bound
                l = mid_day
            else:
                n_sw_i_m_1, sw_indices_i_m_1, sw_t_unx_i_m_1 = retrieve_sw_indices(
                    data_directory, mid_day - dt.timedelta(days=1), data_source
                )
                if n_sw_i_m_1 != 0:
                    # Already crossed in SW
                    # Resume search:
                    u = mid_day
                    continue
                else:
                    enter_sw_time = sw_time_unx[0]
                    time_enter_sw_utc = helper.UNX_to_UTC(enter_sw_time)
                    # break

        # Leaving SW
        if Nsw_end_day == 0:
            if n_sw_i == 0:
                # Has already left sw
                u = mid_day
            else:
                n_sw_i_m_1, sw_indices_i_m_1, sw_t_unx_i_m_1 = retrieve_sw_indices(
                    data_directory, mid_day + dt.timedelta(days=1), data_source
                )
                if n_sw_i_m_1 != 0:
                    # Resume search:
                    l = mid_day
                else:
                    exit_sw_time = sw_time_unx[-1]
                    time_leave_sw_utc = helper.UNX_to_UTC(exit_sw_time)

                    # break
    return n_obs_in_sw, time_enter_sw_utc, time_leave_sw_utc


def monthly_sw_coverage(data_directory, start_date_dt, end_date_dt):

    '''Using monthly ephemeris files, determine the number of
    solar wind observations between start_date and end_date, and if continuous
    daily coverage of the solar wind starts/ends during
    the month, return that time.'''

    # Retrieve month-long anci data file
    month_ephemeris = anc.read_spacecraft_ephemeris(
        data_directory, 'mso',
        start_date=start_date_dt, end_date=end_date_dt,
        restrict_to_timerange=False, verbose=False,
        download_if_not_available=True)

    # Get columns:
    sc_time_unx = month_ephemeris['t']
    sc_time_utc = month_ephemeris["time_utc"]
    x = month_ephemeris['x']
    y = month_ephemeris['y']
    z = month_ephemeris['z']

    # Given the spacecraft position, retrieve the indices
    sw_indices, sw_alt, not_sw_indices = solar_wind_indices(x, y, z)
    n_sw = len(sw_indices)

    # If there's no SW coverage, return [0, 0] (observations of SW),
    # None (Time enter SW), None (Time exit SW)
    if n_sw == 0:
        return [0, 0], None, None

    # Get the POSIX time index at the end of the first day
    # and beginning of the last day:
    first_index_unx = helper.find_closest_index_dt(
        start_date_dt + dt.timedelta(days=1), sc_time_unx)
    last_index_unx = helper.find_closest_index_dt(
        end_date_dt, sc_time_unx)

    # The Number of SW observations in the first day will
    # be the number of SW observations that precede
    # this end (of the first day-of-month) index:
    n_sw_in_first_day = np.where(
        sw_indices < first_index_unx)[0].size
    # Similarly, the Number of SW observations in the last day will
    # be the number of SW observations that follow
    # the beginning (of the last day-of-month) index:
    n_sw_in_last_day = np.where(
        sw_indices > last_index_unx)[0].size

    # Use these as the bracketing information:
    n_obs_in_sw = [n_sw_in_first_day, n_sw_in_last_day]

    # If there are no observations in the first day
    # but are SW obs on the last day, can get the time
    # of the entry as the first entry in sw_indices
    if n_sw_in_first_day == 0 and n_sw_in_last_day > 0:
        time_enter_sw_utc = sc_time_utc[sw_indices[0]]
    else:
        time_enter_sw_utc = None

    # If there are no observations in the last day
    # but are SW obs on the first day, can get the time
    # of the exit as the last entry in sw_indices
    if n_sw_in_last_day == 0 and n_sw_in_first_day > 0:
        time_leave_sw_utc = sc_time_utc[sw_indices[-1]]
    else:
        time_leave_sw_utc = None

    return n_obs_in_sw, time_enter_sw_utc, time_leave_sw_utc


def sw_entry_exit_string(entry_utc, exit_utc, time_in_sw):

    line = "".join(
        (
            "In SW from ",
            entry_utc.strftime("%Y %b %d %H:%M:%S"),
            " to ",
            exit_utc.strftime("%Y %b %d %H:%M:%S"),
            ", duration of SW coverage: ",
            str(time_in_sw),
            ", N Carrington rotations: ",
            str(time_in_sw.days / n_carrington_days),
        )
    )

    return line


def find_sw_coverage(
    data_directory,
    data_source,
    start_search_time,
    end_search_time,
    load_spice_kernels=None,
    print_sw_enter_exit_times=None,
    print_sw_month_scans=None,
    download_spice_if_updated=True
):

    """Retrieve the durations that the MAVEN spacecraft is in the solar wind using
    the Trotignon et al. [2003] fits and the spacecraft position over time as
    retrieved from
     - MAVEN ancillary IDL sav files: Load month-by-month (NOTE: EXTREMELY slow)
     - MAVEN SEP ancillary IDL sav files: Load day-by-day
     - Spice kernels: Make hour-minute resolution times to retrieve."""

    # Clean/convert dates provided
    start_search_date, end_search_date, n_days = helper.sanitize_date_inputs(
        start_date=start_search_time, end_date=end_search_time)

    if data_source == "spice" and load_spice_kernels:
        k = spice.load_kernels(
            data_directory,
            start_date=start_search_time, end_date=end_search_time,
            download_if_not_available=download_spice_if_updated,
            verbose=None)
        # print(spice.currently_loaded_kernels())

    # If the end search time lacks the timezone info,
    # assign to UTC
    if end_search_time.tzinfo is None:
        end_search_time = end_search_time.replace(
            tzinfo=dt.timezone.utc)
    if start_search_time.tzinfo is None:
        start_search_time = start_search_time.replace(
            tzinfo=dt.timezone.utc)

    first_day_of_month_rrule = rrule.rrule(
        rrule.MONTHLY, dtstart=start_search_time,
        until=end_search_time, bymonthday=1)
    last_day_of_month_rrule = rrule.rrule(
        rrule.MONTHLY, dtstart=start_search_time,
        until=end_search_time, bymonthday=-1)

    first_day_each_month = first_day_of_month_rrule.between(
        start_search_time, end_search_time, inc=True)
    last_day_each_month = last_day_of_month_rrule.between(
        start_search_time, end_search_time, inc=True)

    if first_day_each_month[0] > start_search_time:
        first_day_each_month.insert(0, start_search_time)
    if last_day_each_month[-1] < end_search_time:
        last_day_each_month.append(end_search_time)

    # print(first_day_each_month)
    # print(last_day_each_month)
    # input()

    enter_sw_times_utc = []
    exit_sw_times_utc = []

    for day_i, day_f in zip(first_day_each_month, last_day_each_month):
        # Go month by month to find SW coverage
        # Three possibilities:
        # - whole month in SW
        # - whole month not in SW
        # - S/C goes to/from having SW coverage to/from not having SW coverage

        # Get the number of observations in the solar wind
        # on the first and last day of the month (n_obs_in_sw),
        # and the UTC time when daily solar wind coverage stopped
        # (time_leave_sw_utc) or started (time_enter_sw_utc).
        # These are set to None if N/A.
        if data_source == "maven_anci":
            # Use monthly spacecraft ephemeris files that have
            # spacecraft position in MSO to see when the SW coverage
            # starts/stops:
            n_obs_in_sw, time_enter_sw_utc, time_leave_sw_utc =\
                monthly_sw_coverage(data_directory, day_i, day_f)
        else:
            # Otherwise do a binary search on daily files:
            n_obs_in_sw, time_enter_sw_utc, time_leave_sw_utc =\
                search_sw_coverage(data_directory, day_i, day_f, data_source)

        # If neither first day or last day are in SW, skip month
        # If both are, mark as in SW,
        # If one is in Sw but the other is not,
        # use the month dataset to pinpoint the
        # turnover point.
        Nsw_i, Nsw_f = n_obs_in_sw
        # print(n_obs_in_sw, time_enter_sw_utc, time_leave_sw_utc)

        if print_sw_month_scans:
            daterange_str = "{}-{}".format(
                day_i.strftime("%Y %b %d"), day_f.strftime("%d"))
            if Nsw_i > 0 and Nsw_f > 0:
                print("{}: SW coverage in month, next...".format(
                    daterange_str))
            elif Nsw_i == 0 and Nsw_f == 0:
                print("{}: No SW coverage in month, next...".format(
                    daterange_str))
            elif Nsw_i == 0 and Nsw_f > 0:
                enter_str = time_enter_sw_utc.strftime("%Y-%m-%d %H:%M")
                # print(enter_str)
                print("{}: SW coverage gained! Time of coverage"
                      " return: {}".format(daterange_str, enter_str))
            elif Nsw_i > 0 and Nsw_f == 0:
                leave_str = time_leave_sw_utc.strftime("%Y-%m-%d %H:%M")
                # print(leave_str)
                print("{}: SW coverage lost! Time of coverage"
                      " loss: {}".format(daterange_str, leave_str))

        # If already have continuous SW coverage on first day
        # searched, assume has entered that day:
        if Nsw_i > 0 and day_i == first_day_each_month[0]:
            enter_sw_times_utc.append(start_search_time)

        # If the time of contiguous SW coverage ends:
        if time_leave_sw_utc:
            exit_sw_times_utc.append(time_leave_sw_utc)
            last_entry_time = enter_sw_times_utc[-1]
            # print(last_entry_time, time_leave_sw_utc)
            time_in_sw = time_leave_sw_utc - last_entry_time
            if print_sw_enter_exit_times:
                print(sw_entry_exit_string(
                    last_entry_time, time_leave_sw_utc, time_in_sw))

        # If the time of contiguous SW coverage starts:
        if time_enter_sw_utc:
            enter_sw_times_utc.append(time_enter_sw_utc)

        # # If no prev exit times recorded yet:
        # prior_Nsw_i, prior_Nsw_f = prior_n_obs_sw
        # if (Nsw_i > 0 and Nsw_f > 0) and (prior_Nsw_i == 0 and prior_Nsw_f == 0):
        #     if print_sw_month_scans:
        #         print("Entered SW on month-boundary!")
        #     time_enter_sw_utc = current_time

        # if (prior_Nsw_i > 0 and prior_Nsw_f > 0) and (Nsw_i == 0 and Nsw_f == 0):
        #     print("Exited SW on month-boundary!")
        #     time_leave_sw_utc = current_time - dt.timedelta(days=1)
        #     if print_sw_enter_exit_times:
        #         time_in_sw = time_leave_sw_utc - time_enter_sw_utc
        #         print(sw_entry_exit_string(time_enter_sw_utc, time_leave_sw_utc, time_in_sw))
        #     enter_sw_times_utc.append(time_enter_sw_utc)
        #     exit_sw_times_utc.append(time_leave_sw_utc)

        # prior_n_obs_sw = n_obs_in_sw

    return enter_sw_times_utc, exit_sw_times_utc


def make_sw_coverage_file(save_dir, filename, sw_entry_utc, sw_exit_utc):

    header = "Enter SW (YYYY-MM-DD HH:mm), Exit SW (YYYY-MM-DD HH:mm)\n"

    with open(os.path.join(save_dir, filename), "w") as fh:
        fh.write(header)
        for enter_i, exit_i in zip(sw_entry_utc, sw_exit_utc):
            enter_i_str = enter_i.strftime("%Y-%m-%d %H:%M")
            exit_i_str = exit_i.strftime("%Y-%m-%d %H:%M") + "\n"
            fh.write(", ".join((enter_i_str, exit_i_str)))


def read_sw_coverage_file(file):
    sw_entry_utc = []
    sw_exit_utc = []

    with open(file, "r") as fh:
        for line_num, line in enumerate(fh):
            if line_num > 0:
                entry_i, exit_i = line.split(", ")
                exit_i = exit_i.split("\n")[0]
                # print(entry_i, exit_i)
                sw_entry_utc.append(
                    dt.datetime.strptime(entry_i, "%Y-%m-%d %H:%M"))
                sw_exit_utc.append(
                    dt.datetime.strptime(exit_i, "%Y-%m-%d %H:%M"))

    return sw_entry_utc, sw_exit_utc


def make_sw_pass_ranges(sw_coverage_datafile):

    enter_sw_times_utc, exit_sw_times_utc = read_sw_coverage_file(sw_coverage_datafile)
    # print(enter_sw_times_utc, exit_sw_times_utc)

    range_enter_i = []
    range_exit_i = []

    # Need to exclude times the s/c was in safemode -- doesn't count
    safe_mode = [dt.datetime(2022, 2, 22), dt.datetime(2022, 4, 23)]

    for sw_enter_i, sw_exit_i in zip(enter_sw_times_utc, exit_sw_times_utc):
        sw_enter_j = dt.datetime(
            sw_enter_i.year, sw_enter_i.month, sw_enter_i.day)
        sw_exit_j = dt.datetime(
            sw_exit_i.year, sw_exit_i.month, sw_exit_i.day)

        # print(sw_enter_j, sw_exit_j)

        # Safe mode check and exclude:

        # if safe mode ended before the solar wind entry
        # and if safe mode began before solar wind exit
        if (safe_mode[0] < sw_enter_j and safe_mode[0] < sw_exit_j) and (
            safe_mode[1] < sw_enter_j and safe_mode[1] < sw_exit_j
        ):
            # print("Excluded - safe mode before")
            pass
        if (sw_enter_j < safe_mode[0] and sw_exit_j < safe_mode[0]) and (
            sw_enter_j < safe_mode[1] and sw_exit_j < safe_mode[1]
        ):
            # print("Excluded - safe mode after")
            pass
        else:
            # print("safe mode")

            if safe_mode[0] < sw_enter_j and sw_exit_j < safe_mode[1]:
                # If SW coverage began after safe mode start and ended before safe mode ended, skip.
                # print("skip")
                continue
            if sw_enter_j < safe_mode[0] and sw_exit_j < safe_mode[1]:
                # If SW coverage began before safe mode start and ended after safe mode ended,
                # end at the beginning of safe mode
                # print("exit sw coverage at beginning of safe mode")
                sw_exit_j = safe_mode[0]

            if safe_mode[0] < sw_enter_j and safe_mode[1] < sw_exit_j:
                # If SW coverage began after safe mode start and ended after safe mode ended,
                # start indexing from the end of safe mode
                # print("enter sw coverage at end of safe mode")
                # print(safe_mode)
                sw_enter_j = safe_mode[1]
                sw_enter_i = sw_enter_j

            if sw_enter_j < safe_mode[0] and safe_mode[1] < sw_exit_j:

                # print("split sw coverage at end of safe mode")
                # print(safe_mode)

                range_enter_i.append(sw_enter_j)
                range_exit_i.append(safe_mode[0])
                range_enter_i.append(safe_mode[1])
                range_exit_i.append(sw_exit_j)
                continue

        n_days = (sw_exit_i - sw_enter_i).days

        # print(n_days)
        if n_days > 200:
            halfway_time = sw_enter_j + dt.timedelta(days=int(n_days / 2))
            halfway_time = dt.datetime(halfway_time.year, halfway_time.month, halfway_time.day)
            range_enter_i.append(sw_enter_j)
            range_exit_i.append(halfway_time)
            range_enter_i.append(halfway_time)
            range_exit_i.append(sw_exit_j)
        else:
            range_enter_i.append(sw_enter_j)
            range_exit_i.append(sw_exit_j)

        # print(sw_enter_j)

    return range_enter_i, range_exit_i


def make_Carrington_sw_pass_ranges(sw_coverage_datafile):

    enter_sw_times_utc, exit_sw_times_utc = read_sw_coverage_file(sw_coverage_datafile)
    # print(enter_sw_times_utc, exit_sw_times_utc)

    range_enter_i = []
    range_exit_i = []

    for sw_enter_i, sw_exit_i in zip(enter_sw_times_utc, exit_sw_times_utc):
        sw_enter_j = dt.datetime(sw_enter_i.year, sw_enter_i.month, sw_enter_i.day)
        sw_exit_j = dt.datetime(sw_exit_i.year, sw_exit_i.month, sw_exit_i.day)
        n_days = (sw_exit_i - sw_enter_i).days
        n_carrington_cycles = n_days / n_carrington_days
        for i in range(int(np.floor(n_carrington_cycles))):
            sw_exit_j = sw_enter_j + dt.timedelta(days=i * n_carrington_cycles)
            range_enter_i.append(sw_enter_j)
            range_exit_i.append(sw_exit_j)
            sw_enter_j = sw_exit_j

    return range_enter_i, range_exit_i


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data_directory",
        help="Directory containing MAVEN data.",
        required=True
    )
    parser.add_argument(
        "--start_date",
        help="Start date of SW coverage search (YYYY-MM-DD).",
        type=lambda s: dt.datetime.strptime(s, "%Y-%m-%d"),
        default=parse(specification.Mars_orbit_insertion) + dt.timedelta(days=1),
    )
    parser.add_argument(
        "--end_date",
        help="Stop date of comparison timeframe (YYYY-MM-DD).",
        type=lambda s: dt.datetime.strptime(s, "%Y-%m-%d"),
        default=dt.datetime.now() - dt.timedelta(days=1),
    )

    parser.add_argument(
        "--make_file",
        help="Make a file containing all the entry/exit"
             " times of s/c in solar wind.",
        action="store_true",
    )

    parser.add_argument(
        "--data_source",
        help="Read data source ('spice' or 'sep_anci' or 'maven_anci').",
        default="spice",
    )

    parser.add_argument(
        "--download",
        help="Download NAIF files to data_directory if not found.",
        action="store_true",
    )

    parser.add_argument(
        "--show_month_scans",
        help="Show the results of month-to-month"
             " scans (useful for debugging).",
        action="store_true",
    )

    args = parser.parse_args()
    data_directory = args.data_directory
    data_source = args.data_source

    start_search_time = args.start_date
    end_search_time = args.end_date

    # start_search_time = dt.datetime(2014, 12, 1)
    # start_search_time = dt.datetime(2021, 1, 1)
    # start_search_time = dt.datetime(2022, 1, 1)
    # end_search_time = dt.datetime(2015, 4, 1)
    # end_search_time = dt.datetime(2021, 8, 1)
    # end_search_time = dt.datetime(2016, 6, 1)
    # end_search_time = dt.datetime(2021, 8, 25)
    # end_search_time = dt.datetime(2022, 11, 11)

    enter_sw_times_utc, exit_sw_times_utc = find_sw_coverage(
        data_directory,
        data_source,
        start_search_time,
        end_search_time,
        load_spice_kernels=(data_source == "spice"),
        print_sw_enter_exit_times=True,
        print_sw_month_scans=args.show_month_scans,
        download_spice_if_updated=args.download
    )

    if args.make_file:
        start_search_time_str = start_search_time.strftime("%Y%m%d")
        end_search_time_str = end_search_time.strftime("%Y%m%d")
        filename = "maven_sw_coverage_{}to{}.dat".format(
            start_search_time_str, end_search_time_str)
        make_sw_coverage_file(
            data_directory, filename, enter_sw_times_utc, exit_sw_times_utc)
        sw_e, sw_exit = read_sw_coverage_file(
            os.path.join(data_directory, filename))
