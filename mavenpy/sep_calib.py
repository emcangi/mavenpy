import itertools

import numpy as np

# Calibration info:
# The detector pattern referenced in Level 2 files:
l2_coincidence_groupings = ["X", "O", "T", "OT", "F", "FO", "FT", "FTO", "Mixed"]
n_pattern = len(l2_coincidence_groupings)
telescope = ["A", "B"]
n_telescope = len(telescope)
geom = np.array([0.18, 0.0018])
efficiency = 1.

# The height of the shaped pulse output from
# the ADC is directly proportional to the deposited energy.
# This relationship was determined by a calibration
# with the 59.54 keV gamma emission line from Americium-241
# It is reported in ADC units per keV in the SIS, but is
# internally written in SPEDAS mvn_sep_lut2map as
# keV per ADC
# Note that it is fixed to 41 for all coincidence channels

keV_in_adc =\
    {'1A': {"O": 43.77, "T": 38.49, "F": 41.13, "X": 1,
            "FT": 41., "OT": 41., "FTO": 41.},
     '1B': {"O": 41.97, "T": 40.29, "F": 42.28, "X": 1,
            "FT": 41., "OT": 41., "FTO": 41.},
     '2A': {"O": 40.25, "T": 44.08, "F": 43.90, "X": 1,
            "FT": 41., "OT": 41., "FTO": 41.},
     '2B': {"O": 43.22, "T": 43.97, "F": 41.96, "X": 1,
            "FT": 41., "OT": 41., "FTO": 41.}}

# Map refers to
#  mapid 8 -> mapname='Flight2'
#  used from 2014-03-17/22:45 UTC until MOI: 2014-09-22/19:40 UTC
#       (S/C in hybernation since 2014-07-17 UTC until MOI)
#  mapid 9 -> mapname='Flight3'
#  used from 2014-09-22/19:40 UTC until present
mapids = [8, 9]

# Map 9 bin width info
# (same as MAP_ADC_DELTA in Level 2 raw counts)
# (which is the standard for MAVENs mission)
# The raw counts are in 256 bins,
# 128 per telescope, first A and then B
# The bins are currently ordered as O (30 bins),
# then T (12 bins), then F (30 bins), then
# OT (19 bins) and FT (19 bins), and finally FTO (18 bins)
map9_BW_O = [6, 1, 1, 1, 1, 2, 2, 3, 4, 5, 7, 9, 12, 16, 21, 27,
             36, 47, 63, 82, 108, 143, 188, 248, 327, 430, 566,
             746, 993, 1]
map9_BW_F = map9_BW_O
map9_BW_T = [8, 2, 4, 7, 12, 21, 37, 63, 300, 906, 996, 1740]
map9_BW_OT = [8, 2, 4, 7, 12, 21, 37, 27, 36, 47, 63, 190, 331, 575,
              429, 567, 746, 984, 10]
map9_BW_OT = [i*2 for i in map9_BW_OT]
map9_BW_FT = map9_BW_OT
map9_BW_FTO = [8, 2, 4, 7, 12, 21, 16, 21, 27, 36, 110, 190, 331, 248,
               327, 429, 567, 1740]
map9_BW_FTO = [i*4 for i in map9_BW_FTO]

maps = {9: {'telescope_order': ("A", "B"),
            'detector_pattern_order': ("O", "T", "F", "OT", "FT", "FTO"),
            'bin_length': (30, 12, 30, 19, 19, 18),
            "bin_width": {"O": map9_BW_O, "F": map9_BW_F, "T": map9_BW_T,
                          "OT": map9_BW_OT, "FT": map9_BW_FT,
                          "FTO": map9_BW_FTO}}
        }


particle_to_dp = {"ion": "O", "elec": "F"}
p_ld_to_telescope = {"ion_f": "B", "ion_r": "A", "elec_f": "A", "elec_r": "B"}

# Note +1 is added on the index slice
# since Python handles intervals as [i, f)
# rather than [i, f] as in IDL
bin_map_index_slice =\
    {8: {"ion": slice(3, 30 + 1), "elec": slice(3, 17 + 1)},
     9: {"ion": slice(2, None), "elec": slice(2, 16 + 1)}}


def norm_flux_factor(sensor_num, look_direction, particle_name, mapnum=9):

    '''Returns the factor used to calculate the flux for a given
    SEP sensor (#), look director (f or r), and particle name (ion or elec)
    from counts (#/s). mapnum set to default (9).'''

    # Get in lower case format:
    particle_name = (particle_name.lower())[:4]
    look_direction = look_direction.lower()
    sensor_num = str(sensor_num)

    # Get the detector element for the requested particle:
    detector = particle_to_dp[particle_name]

    # Get the telescope name:
    telescope =\
        p_ld_to_telescope["{}_{}".format(particle_name, look_direction)]

    # Energy per ADC for sensor/telescope/detector
    keV_per_adc_i = keV_in_adc["{}{}".format(sensor_num, telescope)][detector]

    # Get the binwidth and subset of bins that apply:
    slice_i = bin_map_index_slice[mapnum][particle_name]
    binwidth_i = maps[mapnum]["bin_width"][detector][slice_i]

    # Calculate the actual energy bin width
    ewidth = 59.5/keV_per_adc_i * np.array(binwidth_i)

    # returns keV cm2 ster

    return ewidth * geom[0]


def dataname_to_subset(data_name='', telescope='', detector_pattern='',
                       t_id='', dp_id='', mapnum=9,
                       telescope_order='', detector_pattern_order='',
                       n_bins=''):

    '''Return the range of indices to access a given
    dataname or telescope + detector_pattern from a map read from
    Level 2 raw or Level 1.
    data_name: string, e.g. 'A-FTO' '''

    if not detector_pattern and not data_name:
        raise ValueError(
            "Need detector_pattern (e.g. 'FTO') & telescope (e.g. 'A')"
            " OR data_name (e.g. 'A-FTO')")
    elif data_name and not data_name:
        potential_delimiters = ["_", " ", "-"]
        for del_i in potential_delimiters:
            if del_i in data_name:
                telescope, detector_pattern = data_name.split(del_i)
                break

    # Get telescope/detector pattern index in map
    if not t_id:
        if not telescope_order:
            telescope_order = maps[mapnum]["telescope_order"]
        t_id = telescope_order.index(telescope)

    if not dp_id:
        if not detector_pattern_order:
            detector_pattern_order = maps[mapnum]["detector_pattern_order"]
        dp_id = detector_pattern_order.index(detector_pattern)

    # Get the n_bins for each t_dp
    if not n_bins:
        n_bins = maps[mapnum]["bin_length"]

    # Grab current location for n_bins of requested
    # and sum all prior to get starting index
    n_bins_t_dp_i = n_bins[dp_id]
    prior_len_dp_i = sum(n_bins[:dp_id])

    lower_lim = t_id * 128 + prior_len_dp_i
    upper_lim = lower_lim + n_bins_t_dp_i

    return lower_lim, upper_lim


def map_to_bins(binname, output='map_array', mapnum=9,
                sensors=('1', '2'), telescopes=("A", "B"),
                detector_pattern=("O", "T", "F", "OT", "FT", "FTO"),
                preceding_str=''):

    # Get map information so can address the ADC array
    BW = maps[mapnum]["bin_width"]
    t_order = maps[mapnum]["telescope_order"]
    dp_order = maps[mapnum]["detector_pattern_order"]
    n_bins = maps[mapnum]["bin_length"]

    # If retrieving the ADC, same for each detector,
    # so only need to iterate over one.
    if binname == "ADC":
        sensors = ("1",)

    # ADC info
    bin_dict = {}

    if output == "map_array":
        names = ["MAP_FTO", "MAP_TID"]
        dtypes = [">i8", ">i8"]

        if binname == "ADC":
            names += ['MAP_ADC_LOW', 'MAP_ADC_HIGH', 'MAP_ADC_AVG',
                      'MAP_ADC_DELTA']
            dtypes += ['>i8', ">i8", 'f', 'f']

        elif binname == "energy":
            names += ["MAP_NRG_MEAS_AVG", "MAP_NRG_MEAS_DELTA"]
            dtypes += ['f', 'f']

        for name_i, dtype_i in zip(names, dtypes):
            bin_dict[name_i] = np.zeros(shape=(256), dtype=dtype_i)

    elif output == "dict":
        if binname == "ADC":
            if not preceding_str:
                label = "{telescope}-{dp_pattern}_{{name}}"
            else:
                label = "{precede}_{telescope}-{dp_pattern}_{{name}}"
        elif binname == "energy":
            if not preceding_str:
                label = "{s}_{telescope}-{dp_pattern}_{{}}"
            else:
                label = "{precede}{s}_{telescope}-{dp_pattern}_{{}}"

    for s, t, dp in itertools.product(sensors, telescopes, detector_pattern):

        # Retrieve telescope / detector pattern index
        t_id = t_order.index(t)
        dp_id = dp_order.index(dp)

        # select indices to subset
        i, f = dataname_to_subset(
            telescope=t, detector_pattern=dp, t_id=t_id, dp_id=dp_id,
            n_bins=n_bins)

        # print(s, t, dp, i, f)

        # Get bin widths
        BW_i = np.array(BW[dp])

        # ADC low/high
        adc_high_i = np.cumsum(BW_i)
        adc_low_i = adc_high_i - BW_i
        adc_avg_i = (adc_high_i + adc_low_i)/2

        # Calculate the energy bins
        if binname == "energy":
            keV_per_adc_s_t = keV_in_adc["{s}{t}".format(s=s, t=t)]
            keV_per_adc_i = keV_per_adc_s_t[dp]
            nrg_meas_avg_i = 59.5/keV_per_adc_i * adc_avg_i
            nrg_meas_delta_i = 59.5/keV_per_adc_i * BW_i

            # Identify where ADC > 4096
            # can cause overflow, so energy bin corrected accordingly.
            overflow = np.where(adc_high_i >= 4096)[0]
            # print(list(adc_high_i))
            # print(overflow)
            # input()
            if len(overflow) != 0:
                overflow = overflow[-1]
                overflow_fudge = 0.3
                nrg_meas_delta_i[overflow] = nrg_meas_avg_i[overflow] * overflow_fudge
                nrg_meas_avg_i[overflow] += nrg_meas_delta_i[overflow] / 2

        # Write the ADC or energy bins to array

        if output == "dict":
            label_i = label.format(
                precede=preceding_str, s=s, telescope=t, dp_pattern=dp)

            if binname == "ADC":
                bin_dict[label_i.format("ADC_avg")] = adc_avg_i
                bin_dict[label_i.format("ADC_low")] = adc_low_i
                bin_dict[label_i.format("ADC_high")] = adc_high_i
                bin_dict[label_i.format("ADC_delta")] = BW_i

            elif binname == "energy":
                bin_dict[label_i.format("energy")] = nrg_meas_avg_i
                bin_dict[label_i.format("denergy")] = nrg_meas_delta_i

        elif output == "map_array":
            bin_dict["MAP_FTO"][i:f] = l2_coincidence_groupings.index(dp)
            bin_dict["MAP_TID"][i:f] = t_id

            if binname == "ADC":
                bin_dict["MAP_ADC_LOW"][i:f] = adc_low_i
                bin_dict["MAP_ADC_HIGH"][i:f] = adc_high_i
                bin_dict["MAP_ADC_AVG"][i:f] = adc_avg_i
                bin_dict["MAP_ADC_DELTA"][i:f] = BW_i

            elif binname == "energy":
                bin_dict["MAP_NRG_MEAS_AVG"][i:f] = nrg_meas_avg_i
                bin_dict["MAP_NRG_MEAS_DELTA"][i:f] = nrg_meas_delta_i

    return bin_dict


def map_data_to_dict(data, telescopes, detector_pattern, data_name,
                     mapnum=9, preceding_str=''):

    '''Deconvolve the 256 bin by N_time array
    into a dictionary with each subset labeled by the telescope
    and detector pattern.'''
    # data: N_time x 256 bin

    counts = {}
    t_order = maps[mapnum]["telescope_order"]
    dp_order = maps[mapnum]["detector_pattern_order"]
    n_bins = maps[mapnum]["bin_length"]

    if not preceding_str:
        label = "{telescope}-{detector_pattern}_{data_name}"
    else:
        label = "{precede}_{telescope}-{detector_pattern}_{data_name}"

    for t, dp in itertools.product(telescopes, detector_pattern):

        # Retrieve telescope / detector pattern index
        t_id = t_order.index(t)
        dp_id = dp_order.index(dp)

        # select indices to subset
        i, f = dataname_to_subset(
            telescope=t, detector_pattern=dp, t_id=t_id, dp_id=dp_id,
            n_bins=n_bins)

        label_i = label.format(
            precede=preceding_str, telescope=t, detector_pattern=dp,
            data_name=data_name)

        counts[label_i] = data[:, i:f]

    return counts


def raw_to_calibrated(raw_data_dict, mapid=9, look_directions=("F", "R"),
                      particle=("elec", "ion"), unit="rates"):

    calibrated_data = {}

    # Retain non-flux/rates arrays
    for n in raw_data_dict:
        if ("flux" not in n and "rate" not in n) and "energy" not in n:
            calibrated_data[n] = raw_data_dict[n]

    # Recover when the attenuator actuates
    atten = raw_data_dict["attenuator_state"]
    atten_actuate = np.where(np.ediff1d(atten) != 0)[0]

    for (dir_i, p_i) in itertools.product(look_directions, particle):

        calib_data_name_i = "{}_{}".format(dir_i, p_i)

        # Get the detector element for the requested particle:
        detector = particle_to_dp[particle_name]

        # Get the telescope name:
        telescope =\
            p_ld_to_telescope["{}_{}".format(look_direction, look_direction)]

        # Get the indices used for the requested particle and map
        index_slice = bin_map_index_slice[mapid][p_i]

        # print(dir_i, p_i, detector_pattern_i, telescope_i)
        raw_data_name_i = "{}-{}".format(telescope_i, detector_pattern_i)

        raw_en_i = raw_data_dict["{}_energy".format(raw_data_name_i)]
        raw_de_i = raw_data_dict["{}_denergy".format(raw_data_name_i)]
        raw_data_i = raw_data_dict["{}_{}".format(raw_data_name_i, unit)]

        # The electronic noise threshold is approximately 11 keV
        # and the amount of energy is digitized with a resolution
        # of 1.36-1.54 keV
        # print(raw_data_name_i, raw_en_i.shape,
        #       raw_de_i.shape, raw_data_i.shape)
        # input()
        calib_en_i = raw_en_i[index_slice] + 10.2
        calib_de_i = raw_de_i[index_slice]
        calib_data_i = raw_data_i[:, index_slice]

        # Next, false counts can be created when the
        # attenuator is opened/closed.
        # Detect when atten changes state and NaN the flux
        # before/after actuation.
        calib_data_i[atten_actuate, :] = np.nan
        calib_data_i[atten_actuate + 1, :] = np.nan
        # calib_data_i[atten_actuate + 2, :] = np.nan

        calibrated_data["{}_energy".format(calib_data_name_i)] = calib_en_i
        calibrated_data["{}_denergy".format(calib_data_name_i)] = calib_de_i
        calibrated_data["{}_{}".format(calib_data_name_i, unit)] = calib_data_i

    return calibrated_data

