import datetime as dt
import re
import os
from threading import Lock as ThreadLock
import sys
import time

from dateutil.parser import parse as parsedt
import posixpath
import requests
from bs4 import BeautifulSoup

from . import helper
from . import specification

# Requires requests, bs4, html5lib


class AuthError(Exception):
    pass


class NotFound(Exception):
    pass


def verify_url(url, session=None, username='', password='',
               verbose=False):

    '''Verify if URL is accessible, and raises AuthError
    if we lack username/password or NotFound if 404'''

    # auth = (username, password)
    if session is None:
        req = requests
    else:
        req = session
    # print(session.auth)
    # print(auth)
    # response = req.get(url, auth=auth)
    response = req.head(url, allow_redirects=True)
    status = response.status_code

    if verbose:
        print("Root URL: ", url, ", status code:", status)

    # input()

    if status < 300:
        return

    if status == 401:
        # 401 Unauthorized
        if not username:
            if session is not None:
                session.close()
            raise AuthError(
                "Need authentication to access '{url}'."
                " Supply username and password.".format(
                    url=url))
        else:
            if session is not None:
                session.close()
            raise AuthError(
                "Username '{}'' and password '' failed to"
                " access '{}'. Check your credentials.".format(url))

    elif status == 404:
        # 404 Not Found
        raise NotFound("URL does not exist: {}".format(url))
    else:
        print("Unknown response code: {}".format(status))


def fix_filename_version(filename, v_number=None, r_number=None):

    '''This function determines if the supplied version number
    and/or revision number is provided (i.e. not None), and if so,
    replaces the section of the string with a zero padded string
    that matches the provided string or integer.

    Useful for setting a fixed version # and/or revision #.'''

    if v_number is not None:
        if "v[0-9][0-9][0-9]" in filename:
            v_substr = "v[0-9][0-9][0-9]"
            v_sigfig = 3
        elif "v[0-9][0-9]" in filename:
            v_substr = "v[0-9][0-9]"
            v_sigfig = 2
        else:
            raise ValueError("Unknown v_number format.")
        v = 'v{}'.format(str(v_number).rjust(v_sigfig, '0'))
        filename = filename.replace(v_substr, v)

    if r_number is not None:
        if "r[0-9][0-9]" in filename:
            v_substr = "r[0-9][0-9]"
            v_sigfig = 2
        else:
            raise ValueError("Unknown r_number format.")
        v = 'r{}'.format(str(r_number).rjust(v_sigfig, '0'))
        filename = filename.replace(v_substr, v)

    return filename


def sdc_retrieve(instrument, destination_dir='',
                 mirrored_local_source_tree="ssl_sprg",
                 start_date=None, end_date=None, orbnum=None, n_days=None,
                 orb_to_t_func=None,
                 level="l2", dataset_name="", ext="", res="", coord="",
                 source="ssl_sprg", username="", password="",
                 v_number=None, r_number=None,
                 verbose=None,
                 prompt_for_download=False):

    '''
    Routine to download MAVEN data from lasp.colorado.edu/maven/sdc
    or sprg.ssl.berkeley.edu/data.

    instrument: string, name of MAVEN instrument with data to be downloaded
    destination_dir: string, path where data will be saved. (required)
    source: string, name of remote to download data from. (required)
    mirrored_local_source_tree: string, if 'ssl_sprg' or 'lasp_sdc',
        will save data to a location within destination_dir that
        mirrors the SPRG or SDC directory tree. if '' or None,
        will not mirror a tree and will just save all the data
        to the provided destination_dir.
    start_date, end_date, orbnum, n_days: string or datetime or int objects,
        defines the span of time to be retrieved.
    orb_to_t_func: function, necessary to determine the file date if
        provided an orbit number.
    username: optional, string, username to access remote if accessing
        embargoed data
    password: optional, string, password to access remote if accessing
        embargoed data
    v_number, r_number: int or str, supply the version and revision number
        to retrieve if already known. (defaults None)
    verbose: Boolean, if True will make print statements useful
        for debugging
    '''

    # Check if the destination directory is provided
    if not destination_dir:
        raise FileNotFoundError(
            "Must provide a folder (destination_dir) where data will"
            " be saved.")
    if not os.path.exists(destination_dir):
        raise FileNotFoundError(
            "Destination directory does not exist, "
            "check input / see if drive connected.")

    # If we already have our local data directory
    # structured to mirror a remote directory tree
    # (e.g. sprg.ssl.berkeley.edu or lasp.colorado.edu),
    # we want to save to a location consistent with that.
    # By default, assumes the ssl_sprg structure. If set to '',
    # assume all files saved direct to destination_dir.
    if mirrored_local_source_tree:
        local_source_folder_root = specification.remote_base_directory(
            mirrored_local_source_tree)

        # Check if the destination dir + source folder exists:
        local_root_dir = os.path.join(
            destination_dir, *local_source_folder_root)
        if not os.path.exists(local_root_dir):
            make_dir_yn = input(
                "No dir available at {}. Do you"
                " want to create one? (Y/N) ".format(local_root_dir))
            if make_dir_yn.lower() == 'n':
                sys.exit("Unable to save data without directory.")

        destination_dir = local_root_dir

    # Get the TLA for the instrument
    tla = instrument[:3].lower()

    # if dataset requested is only a string and not none,
    # make into tuple:
    if isinstance(dataset_name, str):
        dataset_name = (dataset_name,)

    # Iterate through dataset names to get unique
    # paths to batch datasets saved in same dir:
    unique_file_paths = {}
    for d_name in dataset_name:
        # print(d_name, ",")
        # See if the requested dataset exists
        specification.check_if_dataset_exist(
            instrument, level=level, ext=ext,
            dataset=d_name, coord=coord, res=res)

        # And see if it is available on this remote
        specification.check_if_dataset_on_remote(
            tla, source, auth_present=(username != ''),
            level=level, dataset=d_name,
            ext=ext, start_utc=start_date)

        # Get path and filename:
        filedir_tuple = specification.path(
            tla, level, ext=ext, res=res, dataset_name=d_name)
        filename_regex = specification.filename(
            tla, level=level, dataset_name=d_name, ext=ext,
            coord=coord, res=res)
        # print("/".join(filedir_tuple))

        # If supply the version & revision, replace in line:
        filename_regex = fix_filename_version(
            filename_regex, v_number=v_number, r_number=r_number)

        if filedir_tuple in unique_file_paths:
            unique_file_paths[filedir_tuple].append(filename_regex)
        else:
            unique_file_paths[filedir_tuple] = [filename_regex]
    # print(unique_file_paths)
    # input()

    if verbose:
        for path_i in unique_file_paths:
            print("In ", "/".join(path_i))
            print("Search for filenames: ", unique_file_paths[path_i])

    # Make datetime range
    dt_range = helper.daterange(
        start_date=start_date, end_date=end_date, n_days=n_days,
        orbnum=orbnum, orb_to_t_func=orb_to_t_func)

    # Select URL and the subdomain for
    # a given data "source":
    source_folder_root = specification.remote_base_directory(source)

    if "lasp" in source:
        remote_domain = "https://lasp.colorado.edu"
    elif "ssl" in source:
        remote_domain = "http://sprg.ssl.berkeley.edu/data"
    # elif "pds" in source:
    #     remote_domain = "https://pds-ppi.igpp.ucla.edu/"\
    #         "ditdos/download?id=pds://PPI"

    # Iterate through each webpage:
    session = requests.session()
    # print(username, password)
    if username is not None:
        session.auth = (username, password)
    if verbose:
        print("Session opened, iterating through yyyy/mm...")

    # Check if url / source is reachable:
    url_root = "/".join((remote_domain, *source_folder_root))
    verify_url(url_root, session=session, verbose=verbose)

    # Organize the file paths by yyyy/mm since all datasets in SDC/SPRG
    # filed by that in the root dir:
    for filedir_tuple in unique_file_paths:

        # If file dir tuple has yyyy in it, want to see if
        # the URL up to that point exists:
        yyyy_str = '{yyyy}'
        mm_str = '{mm}'
        if yyyy_str not in filedir_tuple:
            yyyy_index = len(filedir_tuple)
        else:
            yyyy_index = filedir_tuple.index(yyyy_str)

        # check if url / source up to the yyyy info is reachable:
        base_url_i = "/".join((url_root, *filedir_tuple[:yyyy_index]))
        verify_url(base_url_i, session=session, verbose=verbose)

        # Check to see if folders for the requested years are on remote:
        if yyyy_str in filedir_tuple:
            # Make the url up to the year:
            yyyy_url = "/".join((url_root, *filedir_tuple[:(yyyy_index + 1)]))

            # If year in url, get all unique years in the requested dt_range:
            uniq_years = set([i.strftime("%Y") for i in dt_range])
            # Iterate over to see if the URLs for each exists:
            avail_years = []
            for yyyy_i in uniq_years:
                yyyy_url_i = yyyy_url.format(yyyy=yyyy_i)
                # print(yyyy_url_i)
                try:
                    verify_url(yyyy_url_i, session=session, verbose=verbose)
                    avail_years.append(yyyy_i)
                except NotFound:
                    pass

            # Restrict the dt range to available years:
            dt_range = [i for i in dt_range if i.strftime("%Y") in avail_years]
            if verbose:
                if len(dt_range) == 0:
                    print("No data available in requested years.")
                else:
                    print("Data available between ",
                          dt_range[0].strftime("%Y-%m-%d"),
                          " and ", dt_range[-1].strftime("%Y-%m-%d"))

        # See if the yyyy/mm folders exist on remote,
        # and restrict the searched time range to when files are available.
        if mm_str in filedir_tuple:
            # Make the url up to the year and month:
            yyyy_mm_url = "/".join((url_root, *filedir_tuple))

            # If year in url, get all unique years in the requested dt_range:
            uniq_year_months = set(
                [(i.strftime("%Y"), i.strftime("%m")) for i in dt_range])
            # print(uniq_year_months)

            avail_yearmonth = []
            for (yyyy_i, mm_i) in uniq_year_months:
                yyyy_mm_url_i = yyyy_mm_url.format(yyyy=yyyy_i, mm=mm_i)
                # print(yyyy_mm_url_i)
                try:
                    verify_url(yyyy_mm_url_i, session=session, verbose=verbose)
                    avail_yearmonth.append((yyyy_i, mm_i))
                except NotFound:
                    pass

            # Restrict the dt range to available years:
            dt_range = [i for i in dt_range if
                        (i.strftime("%Y"), i.strftime("%m"))
                        in avail_yearmonth]
            if verbose:
                if len(dt_range) == 0:
                    print("No data available in requested period.")
                else:
                    print("Data available between ",
                          dt_range[0].strftime("%Y-%m-%d"),
                          " and ",
                          dt_range[-1].strftime("%Y-%m-%d"))

        # If no datetimes remain after the removal, no files to retrieve,
        # abort here.
        if not dt_range:
            print("No data available {} to {}.".format(start_date, end_date))
            session.close()
            return

        # Make the format string for the dataset url:
        remote_url_fstring = "/".join((url_root, *filedir_tuple))

        # print(remote_url_fstring)
        # input()

        destination_dir = os.path.join(local_root_dir, *filedir_tuple)

        # Get the unique year-month combos and track the days in
        # each. The keys of this dict will be used to
        # retrieve the HTML soup of the directory, and the values
        # will be the filenames.

        folder_filename_dict = {}

        # Get regexes for files:
        filename_regex = unique_file_paths[filedir_tuple]

        # print(filedir_tuple)
        # print(filename_regex)
        # input()

        if verbose:
            print("Determine unique year-month keys for"
                  " seeking datafiles corresponding to requested time...")

        for dt_i in dt_range:
            dt_fstring_i = helper.dt_fstring(dt_i)
            yyyy_mm_i = (dt_fstring_i['yyyy'], dt_fstring_i['mm'])

            # Get filename for specific day info.
            filename_i =\
                [regex_i.format(**dt_fstring_i) for regex_i in filename_regex]
            # print(filename_i)
            # input()

            # Get URL and destination folder for year/month.
            # Will index by this, since this is the root downloaded
            # from and will account for same root files (e.g. pfp l0 svy)
            remote_url_i = remote_url_fstring.format(
                yyyy=yyyy_mm_i[0], mm=yyyy_mm_i[1])
            # print(remote_url_i)

            # Add dictionary entry for the remote URL if not created,
            # otherwise add the day to existing entry.
            if remote_url_i in folder_filename_dict:
                # Append to filenames:
                folder_filename_dict[remote_url_i]['filename'] +=\
                    filename_i
            else:
                destination_dir_i = destination_dir.format(
                    yyyy=yyyy_mm_i[0], mm=yyyy_mm_i[1])

                folder_filename_dict[remote_url_i] =\
                    {'filename': filename_i,
                     'destination': destination_dir_i
                     }

        if verbose:
            print("Done.")

        # print(folder_filename_dict)
        # input()
        for url_i in folder_filename_dict:

            if verbose:
                print("URL: {}".format(url_i))
                # print('{}/{}:'.format(*yyyy_mm_i))

            info_i = folder_filename_dict[url_i]
            filenames_i = info_i['filename']
            destination_i = info_i['destination']

            # Get the HTML of the URL, formatted
            # as a "BeautifulSoup object":
            # SLOWEST STEP
            try:
                html_soup_i = html_retrieve(url_i, session=session)
            except requests.exceptions.HTTPError:
                if verbose:
                    print("Files in {} don't exist online yet.".format(url_i))
                break
            # print(html_soup_i)

            # Next, confirm if directory exists, and if not,
            # make it (would have exited if remote didn't have dir.)
            if not os.path.exists(destination_i):
                os.makedirs(destination_i)

            # Get list of local files present in the directory:
            local_files_i = os.listdir(destination_i)

            # For each filename(s) per day in the year/month
            for filename_ij in filenames_i:
                if verbose:
                    print('search for: ', filename_ij)

                # Find any local copy that matches:
                local_copy_ij = [i for i in local_files_i if
                                 re.match(re.compile(filename_ij), i)]

                # Get last modified time of local file and convert
                # to UTC:
                local_modtime_posix =\
                    [os.path.getmtime(os.path.join(destination_i, j))
                     for j in local_copy_ij]
                local_modtime_utc = [dt.datetime.fromtimestamp(
                    j, tz=dt.timezone.utc) for j in local_modtime_posix]

                if len(local_copy_ij) != 0:
                    # Get only the most recently modified file
                    local_copy_ij, local_modtime_utc = newest_file(
                        local_copy_ij, local_modtime_utc,
                        multiple_files_per_day=True)
                else:
                    local_copy_ij = ''
                    local_modtime_utc = 'N/A'

                if verbose:
                    print('local filename: ', local_copy_ij)
                    print('local modtime: ', local_modtime_utc)
                # input()

                # Retrieve most recently modified file(s)
                remote_filename_ij, remote_filename_modtime_dt =\
                    newest_file_from_html_index(
                        html_soup_i, filename_ij,
                        multiple_files_per_day=True)
                if verbose:
                    print("at url: ", url_i)
                    print('remote filename: ', remote_filename_ij)
                    print('remote modtime: ', remote_filename_modtime_dt)
                # input()

                # Get all newest files on remote that are not present locally:
                only_on_remote = [i for i in remote_filename_ij
                                  if i not in local_copy_ij]
                if verbose:
                    print('files only on remote: ', only_on_remote)

                # Download the file(s) if they are not available locally:
                if only_on_remote:
                    if verbose:
                        print("Remote file newer than local, downloading:")
                    if prompt_for_download:
                        input("hit enter to continue")

                    for only_on_remote_i in only_on_remote:
                        only_on_remote_url_i = "/".join((url_i, only_on_remote_i))
                        only_on_remote_destination_i = os.path.join(
                            destination_i, only_on_remote_i)

                        download_file(
                            only_on_remote_url_i,
                            only_on_remote_destination_i,
                            session=session, chunk_size=1048576,
                            verbose=verbose)
                    if verbose:
                        print("Done.")

    session.close()
    if verbose:
        print("Session closed.")


def html_retrieve(url, auth=None, session=None):
    '''Returns a Beautiful Soup object containing
    the HTML of a given url'''

    if session is None:
        req = requests
    else:
        req = session

    http_data_folder = posixpath.join(url)
    r = req.get(http_data_folder, auth=auth)

    # Raise error if unauthorized / page not exist
    r.raise_for_status()

    # Convert the HTML into content
    soup = BeautifulSoup(r.content, 'html5lib')

    return soup


def extract_date_modified(a):
    '''Returns neighboring datemodified
    tag for a given filename as listed on:
    - SSL sprg & NAIF (same HTML structure)
    - LASP SDC'''

    # print('a tag: ', a)
    a_next = a.next_sibling
    # print("next a tag: ", a_next)

    # if no next sibling, use find_next instead
    # (necessary for LASP sdc)
    if not a_next:
        a_datemodified_str = a.find_next().string.strip()
        size_str = a.find_next().find_next().string.strip()

    else:
        a_next = a_next.strip()
        date_modified_size_str = a_next.split("  ")
        # print('date modified (str): ', date_modified_size_str)
        size_str = date_modified_size_str[1].strip()
        a_datemodified_str = date_modified_size_str[0]

    date_modified_dt = parsedt(a_datemodified_str).replace(
        tzinfo=dt.timezone.utc)
    # print('date modified (dt): ', date_modified_dt)
    # print('size :', size_str)
    # input()

    return date_modified_dt


def newest_file_from_html_index(html_soup, filename, verbose=None,
                                ck_check=None, ck_end_dt=None,
                                multiple_files_per_day=None):

    '''Function to return the most recently updated file
    matching a provided 'filename' (can be regex) in an
    html soup object with hrefs that describe each filename.

    html_soup: BeautifulSoup object
    filename: string, can be refex

    '''

    if verbose:
        print("Searching for: ", filename)

    if isinstance(filename, str):
        filename = (filename,)

    newest = []
    last_modified = []

    for i, filename_i in enumerate(filename):
        regex = re.compile(filename_i)
        matching_file_a_tags = html_soup.find_all('a', href=regex)
        # print(filename_i)
        # print(regex)
        # print(matching_file_a_tags)
        # # print(html_soup)
        # input()
        # print(matching_file_a_tags)
        if not matching_file_a_tags:
            newest.append('')
            last_modified.append(None)
            continue

        date_modified =\
            [extract_date_modified(a) for a in matching_file_a_tags]
        # print(date_modified)

        filenames = [a.string for a in matching_file_a_tags]
        # print(filenames)

        newest_file_name_i, most_recent_modified_i = newest_file(
            filenames, date_modified,
            multiple_files_per_day=multiple_files_per_day)
        # print(newest_file_name_i, most_recent_modified_i)

        newest += newest_file_name_i
        last_modified += most_recent_modified_i

        # If we're looking at ck files, we don't want to retrieve
        # the other file_fmts (dailys and quicks) if we have
        # the rel (final) file, which is searched first.

        # Otherwise, we want to restrict the search to files
        # AFTER the last rel file ends.
        if (i == 0 and len(newest_file_name_i) != 0) and ck_check:
            last_ck_name = newest_file_name_i[-1]
            last_yymmdd = dt.datetime.strptime(
                last_ck_name.split("_")[4], "%y%m%d")
            if last_yymmdd >= ck_end_dt:
                break
            # else:
            #     n_days = (ck_end_dt - last_yymmdd).days
            #     daily_dt = [last_yymmdd + dt.timedelta(days=i) for i in range(n_days)]
            #     daily_filfmt = ""
            #     filename[1] = new_red_search

            # "mvn_{platform}_red_" + yymmdd + "_v[0-9][0-9].bc"

    last_modified = [i for i in last_modified if i]
    newest = [i for i in newest if i]
    # print(last_modified, newest)
    # input()

    return newest, last_modified


def newest_file(filenames, date_modified, multiple_files_per_day=None):

    '''Get the newest file for a set of files
    and a set of last-date-modified's.

    Will return the filename corresponding to the most recent
    modification time.

    If set multiple_file_prefix to True, will return
    a list of filenames and times that correspond
    to each unique prefix. This is necessary for
    NGIMS, ACC, and IUVS data
    due to the unpredictability in their names
    (e.g. NGIMS has an arbitrary script name reference
    in the name that is unique to each orbit, e.g.
    '53728' in 'mvn_ngi_l2_ion-abund-53728_20231201T151008_v08_r01.csv'.)

    We want the newest version of each unique file.

    filenames: list of strings, NO REGEX.
    date_modified: dt or posix time
    '''

    most_recent_filenames = []
    most_recent_date_modified = []

    if multiple_files_per_day:

        # When multiple files are present per day e.g. per orbit,
        # they typically have the same leading text up to the
        # substring "_v[0-9][0-9]([0-9])". So we can split on
        # "_v" and get the first part:
        filename_startswith = [i.split('_v')[0] for i in filenames]
        # print(filename_startswith)
        # Select only the unique ones:
        unique_filename_startswith = list(set(filename_startswith))
        # print(unique_filename_startswith)

        for uniq_filename_i in unique_filename_startswith:
            # Get index of filenames that start with that:
            subset_i = [filenames.index(i) for i in filenames
                        if i.startswith(uniq_filename_i)]
            # print(subset_i)
            date_modified_i = [date_modified[i] for i in subset_i]

            # Get most recently modified
            most_recent_date_modified_i = max(date_modified_i)
            most_recent_index_i = date_modified_i.index(
                most_recent_date_modified_i)
            most_recent_filename_i = filenames[subset_i[most_recent_index_i]]

            most_recent_filenames.append(most_recent_filename_i)
            most_recent_date_modified.append(most_recent_date_modified_i)

    else:

        # Grab the most recently modified version if there's multiple versions
        most_recent_date_modified_i = max(date_modified)
        file_index = date_modified.index(most_recent_date_modified_i)

        most_recent_filenames.append(filenames[file_index])
        most_recent_date_modified.append(most_recent_date_modified_i)

    return most_recent_filenames, most_recent_date_modified


def remote_file_updated(remote_file_url, local_file, verbose=None, session=None):

    '''Compare two files to see if the remote has
    been modified after the local.
    remote_file_url: url to remote file
    local_file: path to a local file
    '''

    if session is None:
        req = requests
    else:
        req = session


    # Get last modified time of local file and convert
    # to UTC:
    local_modtime_posix = os.path.getmtime(local_file)
    local_modtime_utc = dt.datetime.fromtimestamp(
        local_modtime_posix, tz=dt.timezone.utc)

    # Likewise for the remote file (assumes Last-Modified
    # is assigned)
    r = req.head(remote_file_url)

    if r.status_code != 200:
        return FileNotFoundError(
            "No file available at '{url}', check "
            "URL and try again.".format(url=remote_file_url))

    remote_modtime_utc = parsedt(r.headers['Last-Modified'])

    # print(remote_modtime_utc, local_modtime_utc)

    if remote_modtime_utc > local_modtime_utc:
        return True
    else:
        if verbose:
            print("No changes to remote file '{remote_file}' "
                  " (last modified {remote_time_utc}) since "
                  "local file '{local_file}' (last modified"
                  " on {local_time_utc}).".format(
                    remote_file=remote_file_url, local_file=local_file,
                    local_time_utc=local_modtime_utc,
                    remote_time_utc=remote_modtime_utc))
        return False


def download_file(file_url, local_filename, session=None,
                  chunk_size=1048576, auth=None, verbose=None):

    '''Function to download a file given the URL (file_url),
    a place to save it (local_filename), and maybe a session
    (which shares auth tokens). Will explode if the status
    code is anything other than 200.

    file_url: string, url of file to be downloaded.
    local_file_name: string ending in file name.
    chunk_size: size of downloaded bytes, default set to 1 MB
    auth: if haven't defined a session, contains username/password
    verbose: '''

    if session is None:
        req = requests
    else:
        req = session

    # Check if directory exist, make if not:
    containing_dir, filename = os.path.split(local_filename)
    if not os.path.exists(containing_dir):
        os.makedirs(containing_dir)

    # Check if file present on remote, exit if not found:
    response = req.head(file_url)
    if response.status_code != 200:
        return FileNotFoundError(
            "No file available at '{url}', check "
            "URL and try again.".format(url=file_url))

    # Download with default chunk size of 1 MB
    size_kb = int(int(response.headers["Content-Length"])/1024)
    # if verbose:
    print("Downloading '{}''...".format(file_url))
    progress_bar = Spinner(size_kb, 'kb')
    tot_i = 0

    try:
        with req.get(file_url, stream=True, auth=auth) as r:
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=chunk_size):

                    # Write:
                    f.write(chunk)

                    # if verbose:
                    # Spinner
                    tot_i += chunk_size/1024
                    if tot_i > size_kb:
                        tot_i = size_kb
                    progress_bar.increment(int(tot_i))
    except KeyboardInterrupt:
        os.remove(local_filename)
        sys.exit()

    print("\nDone.")

    return


class Spinner(object):
    """
    [Written by Autumn Rose Jolitz]
    A Generator is an excellent metaphor for conveying data/state,
    but it always needs to be started. Allocate a spinner generator that
    infinitely cycles throught `states` and zero pads the numbers to have
    constant line length.

    This allows us to "rewrite the line" using the `\r` terminal escape code
    which places the cursor back at the start of the line allowing
    overwriting it.
    """

    def __init__(self, limit, unit_str):
        self.lock = ThreadLock()
        self.last_spin = 0

        def make_generator():
            places = len(str(limit))
            states = ("|", "/", "-", "\\")
            index = 0
            current = 0
            fmt = ("\r{} {:0%sd} / %d " + unit_str) % (places, limit)
            while True:
                index += 1
                index %= 4
                count = yield
                if count is not None:
                    current = count
                sys.stdout.write(fmt.format(states[index], current))
                sys.stdout.flush()

        self.gen = make_generator()
        next(self.gen)  # Whoosh

    def increment(self, value=None):
        with self.lock:
            now = time.time()
            if value is not None:
                self.gen.send(value)
                self.last_spin = now
            elif now - self.last_spin > 0.1:
                next(self.gen)
                self.last_spin = now
