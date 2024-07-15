import os
import argparse
import datetime as dt
from dateutil.parser import parse
from collections.abc import Iterable

from scipy.io import readsav
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge


from mavenpy import spice, mars_shape_conics, coordinates, helper, anc

Rm = mars_shape_conics.Mars_radius


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data_directory",
        help="Directory containing MAVEN data.",
        required=True
    )
    parser.add_argument(
        "--start",
        help="Start of MAVEN plot timeframe (YYYY-MM-DD HH:MM:SS).",
        type=str
    )
    parser.add_argument(
        "--end",
        help="End of MAVEN plot timeframe (YYYY-MM-DD HH:MM:SS).",
        type=str,
    )

    parser.add_argument(
        "--orbit_number",
        help="MAVEN orbit number.",
        type=float
    )

    parser.add_argument(
        "--download",
        help="Download NAIF files to data_directory if not found.",
        action="store_true",
    )

    parser.add_argument(
        "--plot_b",
        help="Add contours for Br.",
        action="store_true",
    )

    parser.add_argument(
        "--plot_path_color",
        help="Color the path by progressing time.",
        action="store_true",
    )

    args = parser.parse_args()

    default_xyzlim = [-3, 3]
    N = 500
    axis_color = '0.75'
    # axis_color = '0.925'

    # Clean the date inputs:
    data_directory = args.data_directory
    # start_date, n_days, end_date = helper.sanitize_date_inputs(
    #     start_date=start, end_date=end)

    if args.orbit_number:
        eph = anc.read_orbit_ephemeris(
            data_directory,
            start_date='2014 12 1', end_date='2024 7 1',
            download_if_not_available=True)
        onum = [args.orbit_number - 0.5, args.orbit_number + 0.5]
        start, end = anc.orbit_num(ephemeris=eph, orbit_num=onum)
    else:
        start = parse(args.start)
        end = parse(args.end)

    print(start, end)

    k = spice.load_kernels(
        data_directory,
        start_date=start, end_date=end,
        download_if_not_available=args.download,
        verbose=None, spk_ext='bsp')
    # print(spice.currently_loaded_kernels())

    delta_t = (end - start).total_seconds()
    sc_time_utc =\
        [start + dt.timedelta(seconds=(delta_t*i/N)) for i in range(N)]
    x_mso, y_mso, z_mso = spice.MAVEN_position(sc_time_utc, frame="MAVEN_MSO")
    x_geo, y_geo, z_geo = spice.MAVEN_position(sc_time_utc, frame="IAU_MARS")

    if args.plot_path_color:
        color_i = cm.viridis(np.linspace(0, 1, N))
        lonalt_color_i = cm.Reds(np.linspace(0, 1, N))
        lonlat_color_i = cm.Blues(np.linspace(0, 1, N))
    else:
        color_i = 'g'
        lonalt_color_i = 'r'
        lonlat_color_i = 'b'

    alt, lat, lon = coordinates.cartesian_to_geographic(
        x_geo, y_geo, z_geo)

    # mso_spline = mars_shape_conics.cartesian_spline(
    #     sc_time_unx, sc_x_mso, sc_y_mso, sc_z_mso)

    # sc_mso_short = sc_mso[start_time_index:end_time_index]
    # rho_mso_scan = np.sqrt(sc_mso_scan[1]**2 + sc_mso_scan[2]**2)/Rm
    x_Rm = x_mso / Rm
    y_Rm = y_mso / Rm
    z_Rm = z_mso / Rm

    # Conics:
    x_bs, y_bs, z_bs = mars_shape_conics.bow_shock()
    x_mpb, y_mpb, z_mpb = mars_shape_conics.MPB()

    fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(10, 6.5))
    path = ((x_Rm, y_Rm, z_Rm), (x_Rm, z_Rm, y_Rm), (y_Rm, z_Rm, -x_Rm))
    bs_line = ((x_bs, y_bs), (x_bs, z_bs), None)
    mpb_line = ((x_mpb, y_mpb), (x_mpb, z_mpb), None)
    plot_x_label = ("MSO X, Rm", "MSO X, Rm", "MSO Y, Rm")
    plot_y_label = ("MSO Y, Rm", "MSO Z, Rm", "MSO Z, Rm")

    for i, ax_i in enumerate(ax[0, :]):
        if i == 2:
            mars_circle = Wedge((0, 0), 1, 0, 360, fc="w", ec="salmon")
            ax_i.add_patch(mars_circle)
        else:
            mars_circle = Wedge((0, 0), 1, -90, 90, fc="w", ec="salmon")
            shadow_hatch = Wedge(
                (0, 0), 1, 90, 90 + 180, fc="w", hatch="//////", ec='salmon')
            ax_i.add_patch(mars_circle)
            ax_i.add_patch(shadow_hatch)

        # Bow shock:
        bs_i = bs_line[i]
        if bs_i is not None:
            plot_x_bs, plot_y_bs = bs_i[0], bs_i[1]
            ax_i.plot(plot_x_bs, plot_y_bs, linestyle="--", color="k")

        # MPB:
        mpb_i = mpb_line[i]
        if mpb_i is not None:
            plot_x_mpb, plot_y_mpb = mpb_i[0], mpb_i[1]
            ax_i.plot(plot_x_mpb, plot_y_mpb, linestyle=":", color="k")

        # path
        plot_x, plot_y, plot_z = path[i]
        index = ((np.sqrt(plot_x**2 + plot_y**2) < 1) & (plot_z > 0))
        # print(index)

        ax_i.scatter(plot_x[~index], plot_y[~index], c=color_i[~index],
                     marker='.', s=4)
        ax_i.scatter(plot_x[index], plot_y[index], c=color_i[index],
                     marker='.', s=0.1)

        ax_i.set_aspect("equal")
        ax_i.set_xlim(default_xyzlim)
        ax_i.set_ylim(default_xyzlim)
        ax_i.set_xlabel(plot_x_label[i])
        ax_i.set_ylabel(plot_y_label[i])
        ax_i.set_facecolor(axis_color)

    # Split the lower row into two plots:
    gs = ax[1, 0].get_gridspec()
    gs01 = gs[1, :].subgridspec(1, 2)
    for ax_i in ax[1, :]:
        ax_i.remove()
    # ax_l = fig.add_subplot(gs[1:, 1:])
    ax_l = fig.add_subplot(gs01[:, 0])
    ax_r = fig.add_subplot(gs01[:, 1])
    ax_r2 = ax_r.twinx()

    # The rho plot:
    rho_bs = np.sqrt(z_bs ** 2 + y_bs ** 2)
    rho_mpb = np.sqrt(z_mpb ** 2 + y_mpb ** 2)
    rho_Rm = np.sqrt(z_Rm**2 + y_Rm**2)

    mars_circle = Wedge((0, 0), 1, -90, 90, fc="w", ec="salmon")
    shadow_hatch = Wedge(
        (0, 0), 1, 90, 90 + 180, fc="w", hatch="////", ec='salmon')
    ax_l.add_patch(mars_circle)
    ax_l.add_patch(shadow_hatch)
    ax_l.plot(x_bs, rho_bs, linestyle="--", color="k")
    ax_l.plot(x_mpb, rho_mpb, linestyle=":", color="k")
    ax_l.scatter(x_Rm, rho_Rm, c=color_i, marker='.', s=2)
    ax_l.set_aspect('equal')
    ax_l.set_xlim(default_xyzlim)
    ax_l.set_ylim([0, default_xyzlim[1]])
    ax_l.set_xlabel("MSO X, Rm")
    ax_l.set_ylabel("ρ [√($Y^2$ + $Z^2$)], Rm")
    ax_l.set_facecolor(axis_color)

    # Lon / Lat

    # Longitude wrapping from 0 to 360 appears as horizontal bars.
    # We can eliminate those by plotting segments where wrapping occurs.
    wrapped_indices = np.where(np.abs(np.ediff1d(lon)) > 180)[0]
    wrapped_indices = np.append(wrapped_indices, len(lon) - 1)

    # input()
    if len(wrapped_indices) == 0:
        ax_r.scatter(lon, lat, color=lonlat_color_i, s=3, marker='.', zorder=2)
        ax_r2.scatter(lon, alt, color=lonalt_color_i, s=3, marker='.', zorder=2)
    else:
        init_index = 0
        for wrap_index in wrapped_indices:
            lon_i = lon[init_index:(wrap_index + 1)]
            lat_i = lat[init_index:(wrap_index + 1)]
            alt_i = alt[init_index:(wrap_index + 1)]

            if isinstance(lonlat_color_i, Iterable):
                lonlat_color_j = lonlat_color_i[init_index:(wrap_index + 1)]
            else:
                lonlat_color_j = lonlat_color_i

            if isinstance(lonlat_color_i, Iterable):
                lonalt_color_j = lonalt_color_i[init_index:(wrap_index + 1)]
            else:
                lonalt_color_j = lonalt_color_i

            ax_r.scatter(lon_i, lat_i, color=lonlat_color_j, s=3, marker='.', zorder=2)
            ax_r2.scatter(lon_i, alt_i, color=lonalt_color_j, s=3, marker='.', zorder=2)
            init_index = wrap_index + 1

    ax_r.set_aspect('equal')
    ax_r.set_xlim([0, 360])
    ax_r.set_yticks([-90 + 30 * i for i in range(7)])
    ax_r.set_ylim([-90, 90])
    ax_r.set_xticks([60 * i for i in range(7)])
    ax_r.set_xlabel("East Longitude, deg.")
    ax_r.set_ylabel("Latitude, deg.", color='b')
    ax_r2.set_ylabel("Altitude, km", color='r')
    ax_r.set_facecolor(axis_color)

    fig.tight_layout()

    if args.plot_b:
        plot_b_file = "/Users/rjolitz/Downloads/Morschhauser_spc_dlat0.25_delon0.25_dalt5.sav"
        b_struct = readsav(plot_b_file)['morschhauser']

        r = b_struct['radius'][0]
        lon = b_struct['longitude'][0]
        lat = b_struct['latitude'][0]
        b = b_struct['b'][0]

        altitude_index = helper.find_closest_index(r, 400 + 3390)
        b_subset = b[:, :, altitude_index, :]
        bx = b_subset[:, :, 0]
        by = b_subset[:, :, 1]
        bz = b_subset[:, :, 2]
        # print(bx.shape, by.shape, bz.shape)

        theta = (np.radians(90.0 - lat))[:, np.newaxis]
        phi = (np.radians(lon))[np.newaxis, :]

        br, bt, bp = coordinates.cartesian_to_spherical_vector(
            bx, by, bz, theta, phi)

        ax_r.contourf(lon, lat, br, vmin=-40, vmax=40, cmap=cm.RdBu, zorder=1)
        # plt.colorbar(label='Br, nT')

    plt.show()
