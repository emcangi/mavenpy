import os
import time

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import cdflib

from mavenpy import swia, helper


def time_uncompress(dim, first_index, compressed, n_iter=100, iterate='time'):

    if iterate == 'time':
        iterate_over_time = True
        iterate_over_index = False
    elif iterate == "index":
        iterate_over_time = False
        iterate_over_index = True

    # Counts over all axes
    init_time = time.time()
    for i in range(n_iter):
        new_counts = swia.uncompress(
            dim, first_index, compressed,
            iterate_over_time=iterate_over_time,
            iterate_over_index=iterate_over_index)
    dur = time.time() - init_time
    dur = np.around(dur, 2)

    print("Counts indexed over "
          "{} (w/ None) (x{}): {} s, Dim: {}".format(
            iterate, n_iter, dur, new_counts.shape))


def plot_beam_spread(time_index, energy, theta, phi, counts):

    counts_e = np.sum(counts, axis=(1, 2))

    fig, ax = plt.subplots(ncols=3, figsize=(12, 5))

    spec = counts_e[time_index, :]

    ax[0].loglog(energy, spec)

    # plt.show()

    # Energy of max flux
    en_index = np.argmax(spec)
    en_max = energy[en_index]
    en_max_round = np.around(en_max, decimals=0)

    above_e = np.where(energy > 1.5*en_max)[0][0]
    print(above_e)

    en_1_index = above_e + np.argmax(spec[above_e:])
    en_1 = energy[en_1_index]
    en_1_round = np.around(en_1, decimals=0)
    # print(en_max_round, en_1_round)
    # plt.show()

    ax[0].set_title("T = {}".format(time_unix[time_index]))
    ax[0].axvline(en_max, color='r')
    ax[0].axvline(en_1, color='k')
    avg_fov = counts[time_index, :, :, en_index]
    # avg_fov = np.sum(
    #     counts[time_index:(time_index + 10), :, :, en_index], axis=0)
    ax[1].set_title("E = {} eV".format(en_max_round))
    p = ax[1].pcolormesh(phi, theta[:, en_index], avg_fov.T, norm=LogNorm())
    plt.colorbar(p)
    # plt.show()
    ax[1].scatter(180, 0, color='r', marker='x')

    ax[2].set_title("E = {} eV".format(en_1_round))
    avg_fov = counts[time_index, :, :, en_1_index]
    p = ax[2].pcolormesh(phi, theta[:, en_1_index], avg_fov.T, norm=LogNorm())
    plt.colorbar(p)
    ax[2].scatter(180, 0, color='r', marker='x')


def multiplot(epoch, theta, phi, energy, counts):

    counts_after_overthetaphi = np.sum(counts, axis=(1, 2))
    counts_after_overenphi = np.sum(counts, axis=(1, 3))
    counts_after_overthetaen = np.sum(counts, axis=(2, 3))

    fig_0, ax_0 = plt.subplots(nrows=3, sharex=True)
    ax_0[0].pcolormesh(
        epoch, theta_fine[:, 0], counts_after_overenphi.T, norm=LogNorm())
    ax_0[1].pcolormesh(
        epoch, energy_fine, counts_after_overthetaphi.T, norm=LogNorm())
    ax_0[1].set_yscale('log')
    ax_0[2].pcolormesh(
        epoch, phi_fine, counts_after_overthetaen.T, norm=LogNorm())

    # time_index = np.searchsorted(time_unix, 1664613911)
    # # time_index = np.searchsorted(time_unix, 1664615000)
    # time_index = np.searchsorted(time_unix, 1664613500)
    # # time_index = np.searchsorted(time_unix, 1664613300)
    # time_index = np.searchsorted(time_unix, 1664613376)

    # time_index = np.searchsorted(time_unix, 1664587240)

    # time_i = time_unix[time_index]

    # for ax_i in ax_0:
    #     ax_i.axvline(time_i)
    #     ax_i.set_xlim([time_i - 200, time_i + 200])


n_iter = 10


swi_dir = "/Users/rjolitz/Downloads/test_maven_data/"

swif_filename = os.path.join(swi_dir, 'mvn_swi_l2_finearc3d_20221001_v02_r01.cdf')
swic_filename = os.path.join(swi_dir, 'mvn_swi_l2_coarsearc3d_20221001_v02_r01.cdf')

# print(dir(swia))

fields = ("time_unix", "epoch", "energy_fine",
          "atten_state",
          "theta_fine", "theta_atten_fine", "phi_fine",
          "estep_first", "dstep_first",
          "diff_en_fluxes", "counts", "de_over_e_fine")

# swif = swia.read(swif_filename, "f", sum_over_fov=True, preprocess=True)
swif = swia.read(swif_filename, "f", fields=fields, sum_over_fov=False,
                 preprocess=True, uncompress_fine=True, lib='cdflib')
# swic = swia.read(swic_filename, "c")

time_unix = swif['time_unix'][0]
epoch_utc = swif['epoch'][0]
estep_first = swif['estep_first'][0]
dstep_first = swif['dstep_first'][0]
energy_fine = swif['energy_fine'][0]
theta_fine = swif['theta_fine'][0]
phi_fine = swif['phi_fine'][0]
phi_first = np.zeros(shape=time_unix.shape).astype('int')

counts = swif['counts'][0]  # N_time x N_phi x N_theta x N_en

print(len(energy_fine), len(theta_fine), len(phi_fine), counts.shape)

time_index = np.searchsorted(time_unix, 1664613911)
# time_index = np.searchsorted(time_unix, 1664587240)
plot_beam_spread(time_index, energy_fine, theta_fine, phi_fine, counts)

# input()

multiplot(epoch_utc, theta_fine, phi_fine, energy_fine, counts)

plt.show()

# Counts over all axes
# None (phi) x theta x energies
time_uncompress(
    (len(phi_fine), len(theta_fine), len(energy_fine)),
    ("None", dstep_first, estep_first),
    counts, n_iter=n_iter, iterate='time')
time_uncompress(
    (len(phi_fine), len(theta_fine), len(energy_fine)),
    ("None", dstep_first, estep_first),
    counts, n_iter=n_iter, iterate='index')
time_uncompress(
    (len(phi_fine), len(theta_fine), len(energy_fine)),
    (phi_first, dstep_first, estep_first),
    counts, n_iter=n_iter, iterate='time')
time_uncompress(
    (len(phi_fine), len(theta_fine), len(energy_fine)),
    (phi_first, dstep_first, estep_first),
    counts, n_iter=n_iter, iterate='index')
input()

# fig, ax = plt.subplots(nrows=3, sharex=True, sharey=True)

# Counts - N_time x N_en
counts_overthphi = np.sum(counts, axis=(1, 2))
time_uncompress(
    (len(energy_fine),),
    (estep_first,),
    counts_overthphi, n_iter=n_iter, iterate='time')
time_uncompress(
    (len(energy_fine),),
    (estep_first,),
    counts_overthphi, n_iter=n_iter, iterate='index')

# init_time = time.time()
# for i in range(n_iter):
#     new_counts = uncompress(
#         (len(energy_fine),), (estep_first,), counts_overthphi,
#         iterate_over_index=True)
# print("Sum(Counts [theta, phi]) iterated over "
#       "energy (x{}): {} s".format(n_iter, time.time() - init_time))


# ax[0].pcolormesh(time_unix, energy_fine, new_counts.T, norm=LogNorm())
# ax[0].set_yscale('log')
# # plt.show()



# init_time = time.time()
# for i in range(n_iter):
#     new_counts = uncompress(
#         (len(phi_fine), len(theta_fine), len(energy_fine)),
#         (phi_first, dstep_first, estep_first), counts,
#         iterate_over_time=True)
# print("Counts indexed over "
#       "time (x{}): {} s".format(n_iter, time.time() - init_time))




init_time = time.time()
for i in range(n_iter):
    new_counts = uncompress(
        (len(phi_fine), len(theta_fine), len(energy_fine)),
        ("None", dstep_first, estep_first), counts,
        iterate_over_index=True)
print("Counts indexed over "
      "index (w/ None) (x{}): {} s".format(n_iter, time.time() - init_time))


init_time = time.time()
for i in range(n_iter):
    new_counts = uncompress(
        (len(phi_fine), len(theta_fine), len(energy_fine)),
        (phi_first, dstep_first, estep_first), counts,
        iterate_over_time=True)
print("Counts indexed over "
      "time (w/ phi index) (x{}): {} s".format(n_iter, time.time() - init_time))

init_time = time.time()
for i in range(n_iter):
    new_counts = uncompress(
        (len(phi_fine), len(theta_fine), len(energy_fine)),
        (phi_first, dstep_first, estep_first), counts,
        iterate_over_index=True)
print("Counts indexed over "
      "index (w/ phi index) (x{}): {} s".format(n_iter, time.time() - init_time))



print(new_counts.shape)
input()


# Counts over phi
counts_overphi = np.sum(counts, axis=1)
init_time = time.time()
for i in range(n_iter):
    new_counts = uncompress(
        (len(theta_fine), len(energy_fine)),
        (dstep_first, estep_first), counts_overphi,
        iterate_over_index=True)
print("Sum(Counts [phi]) iterated over "
      "time (x{}): {} s".format(n_iter, time.time() - init_time))
print(new_counts.shape)
counts_after_overtheta = np.sum(new_counts, axis=1)
counts_after_overen = np.sum(new_counts, axis=2)

ax[2].pcolormesh(
    time_unix, energy_fine, counts_after_overtheta.T, norm=LogNorm())
ax[2].set_yscale('log')

fig_0, ax_0 = plt.subplots(nrows=2, sharex=True)
ax_0[0].pcolormesh(
    time_unix, theta_fine[:, 0], counts_after_overen.T, norm=LogNorm())
ax_0[1].pcolormesh(
    time_unix, energy_fine, counts_after_overtheta.T, norm=LogNorm())
ax_0[1].set_yscale('log')
plt.show()


# Counts - N_time x N_en
init_time = time.time()
for i in range(n_iter):
    new_counts = uncompress(
        (len(energy_fine),), (estep_first,), counts_overthphi,
        iterate_over_index=True)
print("Sum(Counts [theta, phi]) iterated over "
      "energy (x{}): {} s".format(n_iter, time.time() - init_time))

ax[0].pcolormesh(time_unix, energy_fine, new_counts.T, norm=LogNorm())
ax[0].set_yscale('log')
# plt.show()
# input()

init_time = time.time()
for i in range(n_iter):
    new_counts = uncompress(
        (len(energy_fine),), (estep_first,), counts_overthphi,
        iterate_over_time=True)
print("Sum(Counts [theta, phi]) iterated over "
      "time (x{}): {} s".format(n_iter, time.time() - init_time))
ax[1].pcolormesh(time_unix, energy_fine, new_counts.T, norm=LogNorm())
ax[1].set_yscale('log')




