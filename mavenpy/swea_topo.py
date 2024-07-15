import numpy as np

from scipy.io import readsav

from .helper import numpy_UNX_to_UTC


#########################################
#      Magnetic topology                #
#########################################

# Dimensionality of indices give specific topology:
topo_dict = {
    "Unknown": 0,
    "Closed - Dayside": 1,
    "Closed - Day-Night": 2,
    "Closed - Night, Trapped": 3,
    "Closed - Night, Void": 4,
    "Open - Day": 5,
    "Open - Night": 6,
    "Draped": 7,
}


def swe_topo_matrix():

    # Topology matrix
    # Create a 3 x 3 x 3 x 3 x 3 x 3 x 4 matrix that
    # will return the specific topology for a particular combination
    # of parameters. In order of the dimension, the dimension axes are:
    # 0: Shape parameter of away-from-planet directed electrons (upward):
    #    0 - photoelectrons traveling away from planet
    # 1: Shape parameter of toward-planet electrons (downward)
    #    0 - photoelectrons traveling toward planet
    # 2: Flux above or below void threshold
    #    0 - electron void
    # 3: PAD parameter of away-from-planet 100-300 eV electrons (upward):
    #    1 - loss cone
    # 4: PAD parameter of toward-planet 100-300 eV electrons (downward):
    #    1 - loss cone
    # 5: PAD behavior: 0 - loss cone, 1 - isotropic, 2 - upward beamed
    # (If any is this is 2, NaN - no information.)
    # jaway, jtoward, jvoid, score_up, score_down, lc_iso_upward_beam

    topomtrx = np.zeros(shape=(3, 3, 3, 3, 3, 4)).astype("int")

    # Neutralize NaN coded indices
    topomtrx[2, :, :, :, :, :] = topo_dict["Unknown"]
    topomtrx[:2, 2, :, :, :, :] = topo_dict["Unknown"]
    topomtrx[:2, :2, 2, :, :, :] = topo_dict["Unknown"]
    topomtrx[:2, :2, :2, 2, :, :] = topo_dict["Unknown"]
    topomtrx[:2, :2, :2, :2, 2, :] = topo_dict["Unknown"]
    topomtrx[:2, :2, :2, :2, :2, 2] = topo_dict["Unknown"]

    # (up/down)ward photoelectrons + no (down/up)ward photoelectrons
    # + not a void + (up/down)ward losscone + no (downup)ward losscone
    topomtrx[1, 0, 1, 1, 0, :] = topo_dict["Closed - Day-Night"]
    topomtrx[0, 1, 1, 0, 1, :] = topo_dict["Closed - Day-Night"]

    # Open - Day - 5:
    # Down photo electrons + down losscone + up sw electrons ("rare")
    topomtrx[1, 0, 1, 0, 1, :] = topo_dict["Open - Day"]
    topomtrx[1, 0, 1, 2, 1, :] = topo_dict["Open - Day"]

    # Only toward-planet photo-electrons, no losscones by PAD or
    # flux ratio:
    topomtrx[1, 0, 1, 0, 0, 1:3] = topo_dict["Open - Day"]

    # Only toward-planet photo-electrons, no losscones by PAD score,
    # but losscone by flux ratio.
    topomtrx[1, 0, 1, 0, 0, 0] = topo_dict["Closed - Day-Night"]

    # Up photo electrons + no down photo electrons, no down losscone
    topomtrx[0, 1, 1, 1, 0, :] = topo_dict["Open - Day"]
    topomtrx[0, 1, 1, 0, 0, :] = topo_dict["Open - Day"]

    # Open - Night - 6
    # No photo electrons + up onesided losscone
    topomtrx[1, 1, 1, 1, 0, :] = topo_dict["Open - Night"]

    # Draped - 7
    # No photo electrons + down onesided losscone
    topomtrx[1, 1, 1, 0, 1, :] = topo_dict["Draped"]
    # No photo electrons + no losscones
    topomtrx[1, 1, 1, 0, 0, :] = topo_dict["Draped"]

    # Closed - Night, Trapped -- 3:
    # Regardless of photoelectrons, if losscones on both ends by PAD
    # score, mark as trapped.
    topomtrx[:, :, 1, 1, 1, :] = topo_dict["Closed - Night, Trapped"]

    # If there are photoelectrons in both directions and not a void,
    # should be a closed field line connected to dayside on both
    # ends.
    topomtrx[0, 0, 1, :, :, :] = topo_dict["Closed - Dayside"]

    # If photoelectrons in both directions but upward beamed
    # or a losscone:
    topomtrx[0, 0, 1, :, :, 2] = topo_dict["Closed - Day-Night"]
    topomtrx[0, 0, 1, :, :, 0] = topo_dict["Closed - Day-Night"]

    # Closed - Night, Void -- 4:
    # If void, ignore everything else, mark as void.
    topomtrx[:, :, 0, :, :, :] = topo_dict["Closed - Night, Void"]

    return topomtrx


def read_swea_l3_swia_regid(data_directory, observation_datetime):
    # Use to filter out sheath / solar wind regions
    swiid_filename = instrument_data_file(
        data_directory, "swea", observation_datetime, level="3", swe_topo_data="swia_regid"
    )
    swiid = readsav(swiid_filename)["regid"]
    swiid_time = swiid["time"]
    swi_id = swiid["id"]

    return {"time": (swiid_time, "unx"), "swi_id": (swi_id, "0 - sw, 1 - sheath, 2/3")}


def read_swea_l3_topo_parameters(data_directory, observation_datetime):

    """Return SWEA energy fluxes or counts and energy bins for recorded spectra
    in units of 'eflux' or 'counts'.
    Adapted from mvn_swia_load_l2_data.pro.

    data_directory: string, path where MAVEN data is saved.
    observation_datetime: datetime object, time of requested dataset

    Returns dictionary of time, energy, and spectra corresponding to requested dataset.
    """

    # Open required data - shape, padscore, and swia id parameters.

    shape_filename = instrument_data_file(
        data_directory, "swea", observation_datetime, level="3", swe_topo_data="shape"
    )
    shape = readsav(shape_filename)["strday"]
    # print("retrieve shape")

    padscore_filename = instrument_data_file(
        data_directory, "swea", observation_datetime, level="3", swe_topo_data="padscore"
    )
    padscore = readsav(padscore_filename)["padtopo"]
    # print("retrieve pad score")

    # The shape parameter is a 3 x 3 matrix for each time slice.
    # The first dimension describes the described electron population
    # (away/toward/trapped). The second dimension describes the
    # relative pitch angle range, which is:
    # Index 0 : Relative PA 0-30:
    # 0-30 for towards, 150-180 for away, 30-150 for trapped
    # Index 1 : Relative PA 0-45:
    # 0-45 for towards, 135-180 for away, 45-135 for trapped
    # Index 2 : Relative PA 0-60:
    # 0-60 for towards, 120-180 for away, 60-120 for trapped
    # By default, S. Xu's routine uses Index 0 (0-30), so
    # we will use that as default here.
    relative_pa_shape_index = 0

    shape_parameter = np.stack(shape["shape"], axis=0)
    shape_away = shape_parameter[:, relative_pa_shape_index, 0]
    shape_toward = shape_parameter[:, relative_pa_shape_index, 1]

    flux_at_40eV = shape["f40"]

    # plt.figure()
    # plt.plot(numpy_UNX_to_UTC(shape["t"]), shape_away, color='g')
    # plt.plot(numpy_UNX_to_UTC(shape["t"]), shape_toward, color='r')
    # plt.ylim([0, 3])

    shape_threshold = 1
    void_flux_threshold = 1e5

    # For away/towards shape parameters:
    # If shape parameter is below the shape threshold (usually 1),
    # there is photoelectrons - set to 0.
    # Else no photoelectrons - set to 1.
    # Else shape is nan - set to 2
    # For trapped: Determine if flux is below a threshold,
    # making it a void - set to 0 if void.
    # Else not a void - set to 1.
    # Else flux is NaN - set to 2.

    # Shape scoring
    jaway = np.where(shape_away < shape_threshold, 0, 1)
    jtoward = np.where(shape_toward < shape_threshold, 0, 1)
    jvoid = np.where(flux_at_40eV < void_flux_threshold, 0, 1)

    # jaway = np.where(np.floor(shape_away/shape_threshold) < 1, 0, 1)
    # jtoward = np.where(np.floor(shape_away/shape_threshold) < 1, 0, 1)

    jaway[np.isnan(shape_away)] = 2
    jtoward[np.isnan(shape_toward)] = 2
    jvoid[np.isnan(flux_at_40eV)] = 2

    # plt.figure()
    # plt.plot(numpy_UNX_to_UTC(shape["t"]), jtoward, color='r')
    # plt.plot(numpy_UNX_to_UTC(shape["t"]), jaway, color='b')
    # plt.show()

    # Draped v Open-to-night
    # Flux ratio is a 3 x 2 matrix for each time slice.
    # Flux ratio = (flux away from Mars)/(flux toward Mars)
    #  for flux in relative pitch angle range & energy range.
    # The first dimension describes the pitch angle range
    # (see above), the second dimension describes the energy range
    # (low (35 - 60 eV) - 0, high (100-300 eV) - 1).
    # Flux ratio < 0.20: Flux away < 20% of flux toward.
    #   -- set to 0: loss cone
    # Flux ratio < 5: 20% of Flux away < Flux toward.
    #   -- if loss cone and set to 0: loss cone
    #   -- if not loss cone and set to 0: isotropic
    #   -- if loss cone and set to 1: isotropic
    #   -- if not loss cone and set to 1: upward beamed
    energy_range_index = 0

    flux_ratio = np.stack(shape["fratio_a2t"], axis=0)
    ratio_pa_en = flux_ratio[:, relative_pa_shape_index, energy_range_index]

    # plt.figure()
    # plt.plot(numpy_UNX_to_UTC(shape["t"]), ratio_pa_en, color='k')
    # plt.ylim([0.1, 30])
    # plt.yscale('log')

    loss_cone_flux_ratio_threshold = 0.2
    beam_flux_ratio_threshold = 5

    lc_presence = np.where(ratio_pa_en < loss_cone_flux_ratio_threshold, 0, 1)
    beam_presence = np.where(ratio_pa_en < beam_flux_ratio_threshold, 0, 1)
    lc_iso_upward_beam = lc_presence + beam_presence
    lc_iso_upward_beam[np.isnan(ratio_pa_en)] = 3

    # PAD scoring
    # Compare parallel to perpendicular flux.
    # if parallel flux is 3sigma below perpendicular flux, mark
    # as a loss cone.
    # Z = (x - μ)/σ
    # Upward PAD - (Z < -3) - 1 if loss cone, 0 if not.
    # Downward PAD - (Z < -3) - 1 if loss cone, 0 if not.

    loss_cone_threshold = 3

    zscore_down = padscore["zscoredown"]
    zscore_up = padscore["zscoreup"]

    # plt.figure()
    # plt.plot(numpy_UNX_to_UTC(padscore["time"]), zscore_up, color='b')
    # plt.plot(numpy_UNX_to_UTC(padscore["time"]), zscore_down, color='r')
    # plt.ylim([-100, -3])

    up_lc = np.where(zscore_up < -loss_cone_threshold, 1, 0)
    up_lc[np.isnan(zscore_up)] = 2
    down_lc = np.where(zscore_down < -loss_cone_threshold, 1, 0)
    down_lc[np.isnan(zscore_down)] = 2

    # Can replicate the behavior of tplot interp.pro by using numpy interp and flooring:
    shape_time = shape["t"]
    padsc_time = padscore["time"]
    up_lc_shape_coordinate = np.floor(np.interp(shape_time, padsc_time, up_lc)).astype("int")
    down_lc_shape_coordinate = np.floor(np.interp(shape_time, padsc_time, down_lc)).astype("int")

    # plt.show()

    topo_matrx = swe_topo_matrix()
    topo = topo_matrx[
        jaway, jtoward, jvoid, up_lc_shape_coordinate, down_lc_shape_coordinate, lc_iso_upward_beam
    ]

    # SWIA ID match
    # Use to filter out sheath / solar wind regions
    swiid_filename = instrument_data_file(
        data_directory, "swea", observation_datetime, level="3", swe_topo_data="swia_regid"
    )
    swiid = readsav(swiid_filename)["regid"]
    # print("retrieve swi id")

    swiid_time = swiid["time"]
    swiid_id = swiid["id"]

    id_indices = np.searchsorted(swiid_time, shape_time)
    id_indices[id_indices == swiid_time.size] = swiid_time.size - 1
    # id_indices = np.floor(np.interp(shape_time, swiid_time, swiid_id)).astype("int")
    shape_swiid = swiid_id[id_indices]

    sheath_or_sw = np.where((shape_swiid == 1) | (shape_swiid == 2))[0]
    topo[sheath_or_sw] = topo_dict["Draped"]

    # fig, ax = plt.subplots(nrows=5, sharex=True)
    # ax[0].plot(numpy_UNX_to_UTC(padscore["time"]), down_lc, color='b')
    # ax[0].plot(numpy_UNX_to_UTC(padscore["time"]), up_lc, color='r')
    # ax[1].plot(numpy_UNX_to_UTC(shape["t"]), jtoward, color='r')
    # ax[1].plot(numpy_UNX_to_UTC(shape["t"]), jaway, color='b')
    # ax[2].plot(numpy_UNX_to_UTC(shape["t"]), down_lc_shape_coordinate, color='b')
    # ax[2].plot(numpy_UNX_to_UTC(shape["t"]), up_lc_shape_coordinate, color='r')
    # ax[3].plot(numpy_UNX_to_UTC(shape["t"]), lc_iso_upward_beam, color='k')
    # ax[4].plot(numpy_UNX_to_UTC(shape["t"]), shape_swiid, color='k')
    # ax[4].plot(numpy_UNX_to_UTC(swiid_time), swiid_id, color='k', linestyle='--')
    # ax[0].set_xlim([dt.datetime(2017, 9, 15, 21, 57),
    #                 dt.datetime(2017, 9, 15, 22, 3)])
    # ax[0].set_xlim([dt.datetime(2017, 9, 15),
    #                 dt.datetime(2017, 9, 15, 0, 5)])

    return {"time": (shape_time, "unx"), "topo": (topo, "unitless")}
