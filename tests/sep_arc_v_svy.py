import os

import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scipy.io import readsav
from mavenpy.helper import numpy_UNX_to_UTC

test_dir = "/Users/rjolitz/Downloads/test_maven_data"

l1_full_name = "mvn_sep_l1_20141201_1day.sav"
# l1_full_name = "mvn_sep_l1_20160205_1day.sav"
l1_full_name = "mvn_sep_l1_20170912_1day.sav"
l1_full_name = "mvn_sep_l1_20220216_1day.sav"
# # l1_full_name = "mvn_sep_l1_20230101_1day.sav"
# l1_full_name = "mvn_sep_l1_20230722_v004.sav"
# l1_full_name = "mvn_sep_l1_20230723_v005.sav"
l1_full_name = "mvn_sep_l1_20150326_1day.sav"

# Notes on Level 1 info:
# For compare "s[1,2]_[svy,arc]"
# There are many columns, but some are ALWAYS constant
# and others are redundant
# - F0, DSEQCNTR, and VALID fixed at 0
#   (Verified on four files: 2014/12/1, 2017/9/12, 2016/2/5, 2023/1/1)
# - MAPID fixed at 9 since Mars orbit insertion
# - SENSOR fixed at 1 or 2 depending on which detector checked
# - CCODE and SUBSEC same in each detector, altho ccode can be
# different in svy v arc (2017-09-12)
# - SEQ_CNTR varies over time and svy v arc v 1 and 2
# - Delta_time = Duration (except DUR = 0 where delta_time=nan)
# - RATE = COUNTS_TOTAL/DURATION
# So want to look at: SUBSEC, SEQ_CNTR, CCODE, ATT, DURATION,
# COUNTS_TOTAL, DURATION, CFACTOR

# l1_5min_name = "mvn_sep_l1_20160401_5min.sav"
# l1_5min = readsav(os.path.join(test_dir, l1_5min_name))

# print("L1 5min")
# print(l1_5min['s1_svy'].dtype)
# print(l1_5min['s1_arc'].dtype)

l1_all = readsav(os.path.join(test_dir, l1_full_name))

print("L1 all")
print(l1_all.keys())
for ds in l1_all:
    print("{} dtype:".format(ds), l1_all[ds].dtype)

l1_svy1 = l1_all["s1_svy"]
l1_svy2 = l1_all["s2_svy"]
l1_arc2 = l1_all["s2_arc"]
l1_arc1 = l1_all["s1_arc"]

names = l1_svy1.dtype.names
print(names)

skip = ["TIME", "TRANGE", "MET", "ET", "DATA",
        "F0", "VALID", "DSEQCNTR", "MAPID", "SUBSEC"]
plot_names = [i for i in names if i not in skip]

datasets = ['s1_svy', 's1_arc', 's2_svy', 's2_arc']
colors = ['orange', 'blue', 'red', 'purple']

fig, ax = plt.subplots(nrows=len(plot_names), sharex=True)
ax[0].set_title("s[1,2]_[svy, arc] from {}".format(l1_full_name))
for i, name_i in enumerate(plot_names):

    for j, ds in enumerate(datasets):
        if i == 0:
            label = ds
        else:
            label = ""

        print(name_i, l1_all[ds][name_i][:10])
        ax[i].plot(l1_all[ds]["time"], l1_all[ds][name_i],
                   color=colors[j], label=label)
        ax[i].set_ylabel(name_i)
ax[0].legend(ncols=len(datasets))
# plt.show()


# Noise array plot

skip = ["TIME", "MET", "ET", "DATA", "MAPID",
        "CFACTOR", "VALID", "DURATION", "CCODE", "NOISE_RES"]
datasets = ["s1_nse", "s2_nse"]

names = l1_all["s1_nse"].dtype.names
print("s1_nse NAMES")
print(names)
plot_names = [i for i in names if i not in skip]

fig, ax = plt.subplots(nrows=len(plot_names), sharex=True)
ax[0].set_title("s[1,2]_nse from {}".format(l1_full_name))
for i, name_i in enumerate(plot_names):

    for j, ds in enumerate(datasets):
        if i == 0:
            label = ds
        else:
            label = ""

        data_i = l1_all[ds][name_i]

        if data_i[0].shape:
            ncols = data_i[0].shape[0]
            print(name_i, data_i[:10])
            data_i = np.stack(data_i)
            for i_j in range(ncols):
                print(i_j)
                ax[i].plot(l1_all[ds]["time"], data_i[:, i_j])

        else:
            print(name_i, data_i[:10])
            ax[i].plot(l1_all[ds]["time"], data_i,
                       color=colors[j], label=label)
        ax[i].set_ylabel(name_i)



ax[0].legend(ncols=len(datasets))

plt.show()


# datasets = ['s1_svy', 's1_arc', 's2_svy', 's2_arc']
datasets = ['s1_svy', 's2_svy']
fig, ax = plt.subplots(nrows=len(datasets), sharex=True, sharey=True)
ax[0].set_title(l1_full_name)

bins_id = np.linspace(0, 255, 256)

for i, s in enumerate(datasets):
    time_i = l1_all[s]['time']

    duration_i = l1_all[s]['duration']
    data_i = np.array(
        [j.astype('float') for j in l1_all[s]['data']])

    print(s, ", shape: ", data_i.shape)
    print('Min: ', np.nanmax(data_i), ", Max: ", np.nanmin(data_i))

    non_nan = (~np.isnan(time_i))

    time_i = time_i[non_nan]
    duration_i = duration_i[non_nan]
    data_i = data_i[non_nan, :]

    count_rate_i = data_i/duration_i[:, np.newaxis]

    epoch_i = numpy_UNX_to_UTC(time_i)
    ax[i].pcolormesh(epoch_i, bins_id, count_rate_i.T, norm=LogNorm())
    ax[i].set_ylabel(s)

plt.show()
