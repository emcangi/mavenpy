import sys

import numpy as np
from matplotlib import pyplot as plt

from mavenpy import anc, file_path


# If you don't use IDL or have that set,
# replace this with the data directory you have
# or want data saved in.

# Search for the IDL data directory
spedas_data_directory = file_path.get_IDL_data_dir()
spedas_data_directory = None

if spedas_data_directory is None:
    # If spedas env variable not found, prompt for kernel directory:
    data_directory = input(
        "No SPEDAS data directory found. Enter desired location of "
        "spice kernels: ")

    # Do you already have spice kernels downloaded
    # in this folder the way SPEDAS saves them?
    # (e.g. in /misc/spice/naif/MAVEN/kernels)?
    # It will download according to this organization,
    # mirroring NAIF / SSL.
    mirror_spedas_dir_tree = input(
        "Mirror SPEDAS kernel distribution? (y/n): ").lower().strip() == "y"

else:
    _ = input("Found SPEDAS data directory ('{}'), use this folder?"
              " (y/n or control-C to escape): ".format(spedas_data_directory))

    # If found but don't want to use, inquire again:
    if "n" in _.strip().lower():
        data_directory = input(
            "Will not use that directory. Enter desired location of "
            "spice kernels: ")
        mirror_spedas_dir_tree = input(
            "Mirror SPEDAS kernel distribution? (y/n): ").lower().strip() == "y"
    else:
        data_directory = spedas_data_directory
        mirror_spedas_dir_tree = True


# Ask if download additional:
download = input(
    "Do you want to download kernels if not found "
    "locally? (y/n) ").lower().strip() == "y"

# By default, will ask before downloading any new files:
prompt_for_download = input(
    "Ask before download data? (y/n) ").lower().strip() == "y"

# Can set verbose to True if want to
# see process of loading:
verbose = False


# Set a time range of data to look at
# start_date = '2015 3 1'
# end_date = '2015 3 8'
start_date = '2023 12 7'
end_date = '2023 12 11'

# Load spice kernels needed for MAVEN
# analysis:
# Import spice module from mavenpy
# Then load kernels
from mavenpy import spice
k = spice.load_kernels(
    data_directory,
    start_date=start_date, end_date=end_date,
    download_if_not_available=download,
    mirror_spedas_dir_tree=mirror_spedas_dir_tree,
    prompt_for_download=prompt_for_download)
print("Loaded kernels.")

print("List current kernels:")
print(spice.currently_loaded_kernels())
# input()

sc_time_utc, sc_time_unx, x, y, z = spice.load_MAVEN_position(
    start_date, end_date=end_date,
    n_sample_points=400)

alt = np.sqrt(x**2 + y**2 + z**2)

fig, ax = plt.subplots(nrows=3, sharex=True)
fig_2, ax_2 = plt.subplots()

# Retrieve the spacecraft ephemeris:
eph_sc = anc.read_spacecraft_ephemeris(
    data_directory, 'mso', start_date=start_date,
    end_date=end_date, verbose=verbose,
    download_if_not_available=download,
    restrict_to_timerange=True)
sc_alt = np.sqrt(eph_sc['x']**2 + eph_sc['y']**2 + eph_sc['z']**2)

ax[0].plot(sc_time_utc, x + 100)
ax[1].plot(sc_time_utc, y + 100)
ax[2].plot(sc_time_utc, z + 100)
ax_2.plot(sc_time_utc, alt)

ax[0].plot(eph_sc['time_utc'], eph_sc['x'])
ax[1].plot(eph_sc['time_utc'], eph_sc['y'])
ax[2].plot(eph_sc['time_utc'], eph_sc['z'])
ax_2.plot(eph_sc['time_utc'], sc_alt)
ax_2.set_ylabel('Alt. km')
# plt.show()


# Retrieve orbit ephemeris for this time period
eph = anc.read_orbit_ephemeris(
    data_directory, start_date=start_date, end_date=end_date,
    download_if_not_available=download)
print("Entries in the ephemeris: ", eph.keys())

# Calculate the orbit number for a time in the middle
time_utc_i = '2015 3 3 14:00'

orb_num_i = anc.orbit_num(time_utc=time_utc_i, ephemeris=eph)

print("Time UTC ", time_utc_i, " corresponds to ", orb_num_i)

# Can also go backwards to get a time!
time_alt_i = anc.orbit_num(orbit_num=orb_num_i, ephemeris=eph)

print("Maps back to time UTC ", time_alt_i)
# orb_num = anc.orbit_num(time_unix=eph_sc['t'], ephemeris=eph)


# Create a secondary axis with orbit number.

anc.add_orbit_axis(ax[0], ephemeris=eph, label='Orb. #')
anc.add_orbit_axis(ax_2, ephemeris=eph, label='Orb. #')
# input()

# plt.plot(orb_num, sc_alt)
plt.show()
