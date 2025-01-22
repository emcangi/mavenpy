import os
import collections.abc
import re


import numpy as np

from scipy.io import readsav


from spacepy import pycdf
import cdflib

# from .file_path import instrument_data_file, instrument_data_directory


def read_cdf(cdf_file_path, field_names='', lib='cdflib', show_info=False):

    data_dict = {}

    if lib == 'spacepy':
        data_cdf = pycdf.CDF(cdf_file_path)

        if not field_names:
            field_names = [i for i in data_cdf.keys()]

        for field_name_i in field_names:
            data_i = data_cdf[field_name_i][...]
            # if field_name_i == "epoch":
            #     time_unx = [i for i in data_cdf if "unix" in i or "unx" in i][0]
            #     print(data_i[:10])
            #     input()
            #     data_i = numpy_UNX_to_UTC(data_cdf[time_unx][...])
            data_dict[field_name_i] = data_i

        data_cdf.close()

    elif lib == 'cdflib':
        data_cdf = cdflib.CDF(cdf_file_path)

        if show_info:
            print(data_cdf.cdf_info())
        if not field_names:
            field_names = data_cdf.cdf_info().zVariables

        for field_name_i in field_names:
            data_i = data_cdf.varget(field_name_i)

            # print(field_name_i, type(data_i))

            if field_name_i == "epoch":
                # print("epoch")
                # print(data_i[:10])
                # NOTE: cdfepoch lacks fractional second precision

                data_i = cdflib.cdfepoch.to_datetime(data_i)

                # NaTs fixed as of 2025/01/15 in cdflib 1.3.3.
                # # This is a patch around fillvals in epoch
                # # ruining the rest of the array for heretofor
                # # unknown reasons.
                # # data_i = np.array(data_i)
                # nonsense_index = np.where(data_i < 0)[0]
                # print(nonsense_index)
                # if nonsense_index.size == 0:
                #     data_i = cdflib.cdfepoch.to_datetime(data_i)
                # else:
                #     full_epoch = np.empty(data_i.shape, dtype='datetime64[s]')
                #     not_nonsense_index = np.where(data_i > 0)[0]
                #     not_nonsense_epoch = cdflib.cdfepoch.to_datetime(
                #         data_i[not_nonsense_index])
                #     full_epoch[not_nonsense_index] = not_nonsense_epoch
                #     full_epoch[nonsense_index] = np.datetime64('1900-01-01')
                #     data_i = full_epoch

                # Fixes errant conversion into list, which breaks append:
                data_i = np.array(data_i)
                # # print(data_i[28:34])
                # # input()

            # print(field_name_i, type(data_i))

            data_dict[field_name_i] = data_i
        # print('end varget')
        # input()

    return data_dict


def read_tplot(tplot_file_path, return_plot_parameters=True, verbose=None):

    '''Reads a .tplot file (which is an IDL sav file that contains
    a list of tplot variables, each with parameters used to make
    the plot in IDL). Outputs the datasets (in data_dict)
    and the plot_parameters (if requested, which it will by default).'''

    # Open the tplot file, which is a stylized
    # IDL sav file:
    data_sav = readsav(tplot_file_path)

    # Get the dataset info and # of datasets
    data = data_sav["dq"]
    n_datasets = len(data)

    # Set up dict to contain data and plot limits:
    data_dict = {}

    # If requested, make structure to hold the info for
    # plotting the variables. This is useful for replicating
    # the actual plots, e.g. tohban.
    if return_plot_parameters:
        # Get the plot options:
        tplot_options = data_sav["tv"][0]['options']
        tplot_options_columns = tplot_options.dtype.names
        # print(tplot_options_columns)
        # print(tplot_options['varnames'])
        # print("VARNAMES" in tplot_options_columns)
        # input()

        plot_info_dict = {}

        if "VARNAMES" in tplot_options_columns:
            # Get the plotted names in order of retrieved.
            plot_name_order =\
                [i.decode() for i in tplot_options['varnames'][0]]
            plot_info_dict["tplot_order"] = plot_name_order

    # Iterate through the referenced datasets:
    for i in range(n_datasets):
        struct_i = data[i]

        # Get the name:
        name_i = struct_i["name"].decode()
        if verbose:
            print(name_i)

        # Get the dtype (=1 if has data, =3 if reference to other datasets:)
        dtype_i = struct_i["dtype"]

        # Data-containing struct:
        data_i = struct_i["dh"]

        if dtype_i == 3:
            # This dtype is when the tplot variable
            # is a list of other tplot names,
            # so return the reference names:
            data_dict[name_i] = [j.decode() for j in data_i]

        else:
            # Select the first element of the data structure,
            # since IDL sav always has it as a tuple:
            data_i = data_i[0]

            # In tplot variables, the x coordinate is always a Posix time.
            time_unx = data_i['x']
            data_dict_i = {"time_unix": time_unx}

            # Get the column names as in the data struct:
            data_column_names = data_i.dtype.names

            if "V" in data_column_names:
                data_dict_i['y'] = data_i["V"]
                data_dict_i['z'] = data_i["Y"]
                if verbose:
                    print("Y dim: ", data_dict_i['y'].shape, ", Z dim:",
                          data_dict_i['z'].shape)
            elif "Y" in data_column_names:
                data_dict_i['y'] = data_i["Y"]
                if verbose:
                    print("Y dim: ", data_dict_i['y'].shape)

            data_dict[name_i] = data_dict_i

        # input()
        if return_plot_parameters:
            # Get the limit & dlimit structures:
            plot_info_dict_i = parse_tplot_plot_parameters(
                struct_i, verbose=verbose)
            plot_info_dict[name_i] = plot_info_dict_i

    if return_plot_parameters:
        return data_dict, plot_info_dict
    else:
        return data_dict, 0


def parse_tplot_plot_parameters(data, verbose=None):

    # Grab the limit and dlimit and merge the structures:
    limit = data["lh"]
    dlimit = data["dl"]

    # Merge limit and dlimit bc we don't care which is used,
    # and names/ranges/etc are sometimes in one or the other.
    dlimit_names = (dlimit.dtype.names or ())
    limit_names = (limit.dtype.names or ())
    all_lim = {}
    for name in limit_names:
        limit_i = limit[name][0]
        if isinstance(limit_i, np.ndarray):
            if limit_i.size == 1:
                limit_i = limit_i[0]
        if isinstance(limit_i, bytes):
            limit_i = limit_i.decode()
        all_lim[name] = limit_i
        # print(name, limit_i, type(limit_i))

    for name in dlimit_names:
        dlimit_i = dlimit[name][0]
        if isinstance(dlimit_i, np.ndarray):
            if dlimit_i.size == 1:
                dlimit_i = dlimit_i[0]
        if isinstance(dlimit_i, bytes):
            dlimit_i = dlimit_i.decode()
        # print(name, dlimit_i, type(dlimit_i))
        all_lim[name] = dlimit_i
    all_lim_names = [i for i in all_lim]
    # print(all_lim_names)
    # input()

    # for i in all_lim:
    #     print(i, ":", all_lim[i])

    if "SPEC" in all_lim_names:
        axes = ("y", "z")
    else:
        axes = ("y",)

    # Final plot info dictionary:
    plot_info_dict_i = {}

    for axis_name in axes:
        axis = axis_name.upper()

        title_i = "{}TITLE".format(axis)
        subtitle_i = "{}SUBTITLE".format(axis)
        range_name = "{}RANGE".format(axis)
        log_toggle = "{}LOG".format(axis)

        # if "SPEC" in spec_toggle and "ZTITLE" not in all_lim_names:
        #     zname_i = "<>"

        if title_i in all_lim_names:
            title_str = all_lim[title_i]
        else:
            title_str = ""

        if subtitle_i in all_lim_names:
            # print("subtitle present")
            title_str += '\n'
            limit_substr = ""
            # print("subtitle in limit struct")
            limit_substr = all_lim[subtitle_i]
            # print(limit_struct[subtitle_i][0])
            # print(limit_substr)
            title_str += limit_substr

        # print(ytitle)
        # Replace !C with newline character
        title_str = "\n".join(title_str.split("!C"))

        title_str = title_str.replace("e!E-!N", "Elec")

        # Replace !N and !E with the on/off carat for the tecplot
        if "!N" in title_str:
            title_str = re.sub("!N", "}",  re.sub("!E", "^{", title_str))
            title_str = "".join(("$", title_str, "$"))

        if verbose:
            print(axis_name, ":", title_str)

        plot_info_dict_i["{}label".format(axis_name)] = title_str

        # If there's a range provided but the start=end,
        # act as if None:
        if range_name in all_lim_names:
            range_i = all_lim[range_name]
            if range_i[0] != range_i[1]:
                plot_info_dict_i["{}range".format(axis_name)] = range_i
            if verbose:
                print("{}range: ".format(axis_name), range_i)

        if log_toggle in all_lim_names:
            plot_info_dict_i["{}scale".format(axis_name)] = 'log'
        else:
            plot_info_dict_i["{}scale".format(axis_name)] = 'linear'

    if "CONSTANT" in all_lim_names:
        plot_info_dict_i["constant"] = all_lim["CONSTANT"]

    if "AXIS" in all_lim_names:
        axis_dl = all_lim["AXIS"]

        second_ax = {}
        dtype_names = axis_dl.dtype.names

        for sec_ax_i in ("YTITLE", "YRANGE"):

            if sec_ax_i in dtype_names:
                sec_yabel = axis_dl[sec_ax_i]

                if "TITLE" in sec_ax_i:
                    key_i = "ylabel"
                if "RANGE" in sec_ax_i:
                    key_i = "yrange"

                if isinstance(sec_yabel, np.ndarray):
                    if sec_yabel.size == 1:
                        sec_yabel = sec_yabel[0]
                if isinstance(sec_yabel, bytes):
                    sec_yabel = sec_yabel.decode()

                second_ax[key_i] = sec_yabel
        # print(axis_dl)
        # input()
        plot_info_dict_i["second_axis"] = second_ax
    if verbose:
        print(plot_info_dict_i)
    # input()

    return plot_info_dict_i


def read_sav(sav_file_path, field_names='', struct_name='',
             flatten_struct=True):

    '''Reads an IDL sav file using scipy.io.readsav
    and iterates through the elements to remove tuple indexing,
    allowing easier access of the data.

    sav_file_path: location of IDL sav file
    field_names: the names of fields in the IDL structure,
        if already known.
    struct_names: names of the saved structures, if
        already known and if there are more than one.
    flatten_struct: By default, will flatten all structures
        into one final dict with the fields pre-pended by
        the struct name. If set to False, will instead save
        structs as separate dicts in the final dict.'''

    # Open IDL sav file, which is read into
    # a dictionary containing structure names,
    # each structure name corresponding to a
    # record array.
    data_sav = readsav(sav_file_path)

    # Retrieve struct name if not provided
    if not struct_name:
        available_struc_names = [i for i in data_sav]
        if not struct_name:
            # If left as default null string, set to the list of
            # available struct names
            struct_name_arr = available_struc_names

    elif isinstance(struct_name, str):
        # Make a list of the struct name
        struct_name_arr = [struct_name]
    elif isinstance(struct_name, tuple) or\
            isinstance(struct_name, list):
        struct_name_arr = struct_name
    else:
        file_name = os.path.split(sav_file_path)[1]
        raise NameError(
            "Unrecognized format of structure name for retrieval"
            " from IDL sav file {filename}, provide string or"
            " list/tuple, else default to pulling all structs.".format(
                filename=file_name))

    # Want to avoid overwriting field names
    # as read per structure if not previously defined
    retrieve_field_names_per_str = not field_names

    # print(struct_name_arr)

    data_dict = {}

    for struct_name_i in struct_name_arr:
        # Fill data sav
        data_sav_i = data_sav[struct_name_i]
        # if len(data_sav) != 1:

        # print(struct_name_i)

        # If the structure is not a record array or record,
        # save the contents directly and don't iterate.
        if not isinstance(data_sav_i, (np.recarray, np.record)):
            # print("Skipping record flattening on {}".format(struct_name_i))
            isbytstr, data_sav_i = check_convert_byte_str_arr(data_sav_i)
            data_dict[struct_name_i] = data_sav_i
            continue

        # Retrieve fields
        if retrieve_field_names_per_str:
            field_names = [i for i in data_sav_i.dtype.names]
        # print(field_names)

        # Fill data dict:
        # If you are pulling more than one structure from
        # an IDL sav file, the variable names inside will
        # be pre-pended by the structure name in order
        # to avoid overwriting datasets with the same name
        # e.g. SEP Level 1 data has a column called "DATA"
        # in six structures: s[1,2]-[nse,svy,arc].

        if len(struct_name_arr) != 1 and flatten_struct:
            precede_str_i = struct_name_i
        else:
            precede_str_i = ''

        data_dict_i = reformat_sav(
            data_sav_i, field_names, precede_str=precede_str_i)

        if flatten_struct:
            data_dict.update(data_dict_i)
        else:
            data_dict[struct_name_i] = data_dict_i
        data_dict_i = {}
        data_sav_i = {}

    return data_dict


def check_convert_byte_str_arr(data_i):

    # byte array handling
    # print(data_i)
    # print(type(data_i))

    # if isinstance(data_i,  (int, str, float, bytes)):
    if not isinstance(data_i, collections.abc.Sized) or\
            isinstance(data_i, (bytes, str)):
        # Leave alone if its a lone element
        # Although format it if a byte string

        # print("data is byt/str/float/int")
        if isinstance(data_i, bytes):
            data_formatted = data_i.decode()
        else:
            data_formatted = data_i

        return True, data_formatted

    #
    data_i_0 = data_i[0]

    if isinstance(data_i_0, (bytes, str)):
        # print("element is byt/str/float/int")

        if isinstance(data_i_0, bytes):
            data_formatted = [i.decode() for i in data_i]

        data_formatted = data_i

        if len(data_formatted) == 1:
            data_formatted = data_formatted[0]

        return True, data_formatted
    # print("No change")

    return False, data_i


def reformat_sav(old_data_dict, names, iter=0, precede_str=''):

    new_data_dict = {}

    # print(names)
    # print(old_data_dict.dtype)

    for og_name in names:
        # always lower case the fieldnames
        new_name = og_name.lower()

        if precede_str:
            if precede_str != "time":
                new_name = "{}_{}".format(precede_str, new_name)

        data_i = old_data_dict[og_name]

        # print(old_data_dict.dtype)
        # print(og_name, "->", new_name)
        # print(data_i.shape, data_i.dtype, data_i.dtype.names)

        # Flatten the data element
        if len(data_i) == 1:
            data_i = flatten_rec(data_i)
            # print(data_i)
            # print(type(data_i))
            # print("Flattened shape, dtype, dtypenames: ",
            #       data_i.shape, data_i.dtype, data_i.dtype.names)

        # byte array handling
        isbytstr, data_i = check_convert_byte_str_arr(data_i)

        if isbytstr:
            new_data_dict[new_name] = data_i
            continue

        # print(data_i)

        # If any named, read record array accordingly:
        if data_i.dtype.names:
            names = data_i.dtype.names

            subset_dict = reformat_sav(
                data_i, names, precede_str=og_name.lower())

            for i in subset_dict:
                new_data_dict[i] = subset_dict[i]

            # print(subset_dict.keys())
            # input()
        else:

            if data_i[0].size != 1:
                data_i = np.stack(data_i)
                # print("0th element contains arr, stacked")
                # print("New shape, dtype, dtypenames:",
                #       data_i.shape, data_i.dtype, data_i.dtype.names)

            new_data_dict[new_name] = data_i

    return new_data_dict


def flatten_rec(data):
    '''The IDL sav reader produces a record array
    that can have a bunch of zero indexes. This routine
    flattens so only the rec array / tuple / record
    has N elements.
    '''

    # Subscript at 0 until you get the multielement
    # length
    while True:
        # print(len(data), data.dtype.names)
        # input()
        if len(data) == 1:
            data = data[0]
        else:
            break
    # input()
    return data
