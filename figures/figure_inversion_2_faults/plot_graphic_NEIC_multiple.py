# -*- coding: utf-8 -*-
"""The routines here allow to plot the solution model FFM modelling, as well
as moment rate function, and waveform fits.
"""

import argparse
import collections
import errno
import glob
import json
import os
import pathlib
from datetime import datetime
from shutil import move
from typing import Any, List, Optional, Tuple, Union

import cartopy.crs as ccrs  # type: ignore
import cartopy.feature as cf  # type: ignore
import cartopy.io.shapereader as shpreader  # type: ignore
import matplotlib  # type: ignore
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import pygmt  # type: ignore
from matplotlib import ticker  # type: ignore
from matplotlib import colormaps, colors, gridspec  # type: ignore
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap  # type: ignore
from matplotlib.image import AxesImage  # type: ignore
from obspy.imaging.beachball import beach, beachball  # type: ignore
from pyproj import Geod  # type: ignore
from pyrocko import moment_tensor as pmt  # type: ignore
from pyrocko import plot  # type: ignore
from pyrocko.plot import beachball  # type: ignore
from scipy.interpolate import griddata  # type: ignore

#
# local modules
#
import ffm.fault_plane as pf
import ffm.plane_management as pl_mng
import ffm.seismic_tensor as tensor
import ffm.shakemap_tools as shakemap
import ffm.velocity_models as mv
from ffm import get_outputs, load_ffm_model
from ffm.plot_maps_NEIC import plot_map, set_map_cartopy
from ffm.static2fsp import static_to_fsp
from ffm.static2srf import static_to_srf

# from ffm.waveform_plots_NEIC import plot_waveform_fits
from ffm.waveform_plots_NEIC import add_metadata
from typing import List, Literal, Optional, Union

plt.rc("axes", titlesize=14)
plt.rc("axes", labelsize=12)
plt.rc("xtick", labelsize=12)
plt.rc("ytick", labelsize=12)
plt.rc("font", size=12)


def plot_waveforms(
    axes: List[plt.Axes],
    times: List[Union[list, np.ndarray]],
    waveforms: List[Union[list, np.ndarray]],
    weights: List[float],
    type_str: Optional[str] = None,
    comp: Optional[List[str]] = None,
    color: str = "blue",
    custom: Optional[Literal["fill", "syn"]] = None,
    yrel_annotation_max: float = 0.6,
) -> List[plt.Axes]:
    """Plot timeseries waveform data

    :param axes: The axes to plot the data on
    :type axes: List[plt.Axes]
    :param times: The time
    :type times: List[np.ndarray]
    :param waveforms: The waveform values
    :type waveforms: List[np.ndarray]
    :param weights: The waveform weight
    :type weights: List[float]
    :param type_str: The data type, defaults to None
    :type type_str: Optional[str], optional
    :param comp: The component, defaults to None
    :type comp: Optional[List[str]], optional
    :param color: The color of the plot line, defaults to "blue"
    :type color: str, optional
    :param custom: The custom option (fill), defaults to None
    :type custom: Optional[Literal[&quot;fill&quot;, &quot;syn&quot;]], optional
    :return: The updated axes
    :rtype: List[plt.Axes]
    """
    nPlot = 0
    for ax, time, waveform, weight in zip(axes, times, waveforms, weights):
        nPlot += 1
        if weight == 0.0:
            ax.plot(time, waveform, color=color, linewidth=2, linestyle="dashed")
        else:
            ax.plot(time, waveform, color=color, linewidth=2 * weight)
        min_time, max_time = ax.get_xlim()
        min_time = np.minimum(np.min(time), min_time)
        max_time = np.maximum(np.max(time), max_time)
        if custom == "fill":
            min_val, max_val = ax.get_ylim()
            min_val = np.minimum(np.min(waveform), min_val)
            max_val = np.maximum(np.max(waveform), max_val)
            ax.set_ylim((-(max(abs(min_val), max_val)), max(abs(min_val), max_val)))
            min_val, max_val = ax.get_ylim()
            ax.vlines(0, min_val, max_val, "k", lw=1)
            if type_str == "body":
                ax.text(
                    np.max(time),
                    yrel_annotation_max * max_val,
                    "{:0.1f}".format(max(abs(min_val), max_val)),
                    ha="right",
                    va="center",
                )
                ax.hlines(0, -20, np.max(time), "k", lw=1)
                ax.set_xlim((-20, np.max(time)))
            elif type_str == "surf":
                ax.text(
                    np.max(time),
                    yrel_annotation_max * max_val,
                    "{:0.3f}".format(max(abs(min_val), max_val)),
                    ha="right",
                    va="center",
                )
                ax.hlines(0, -350, np.max(time), "k", lw=1)
                ax.set_xlim((-350, np.max(time)))
            elif type_str == "cgps" or type_str == "strong":
                min_wval = np.min(waveform)
                max_wval = np.max(waveform)
                if max_wval > abs(min_wval):
                    ax.text(
                        np.max(time),
                        yrel_annotation_max * max_val,
                        "{:0.2f}".format(max_wval),
                        ha="right",
                        va="center",
                    )
                else:
                    ax.text(
                        np.max(time),
                        yrel_annotation_max * max_val,
                        "{:0.2f}".format(min_wval),
                        ha="right",
                        va="center",
                    )
                ax.hlines(0, -15, np.max(time), "k", lw=1)
                ax.set_xlim((-15, np.max(time)))
            min_time, max_time = ax.get_xlim()
            if type_str == "body" and comp == "BHZ":
                ax.text(
                    1.1 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "P",
                    ha="right",
                    va="bottom",
                )
            if type_str == "body" and comp == "SH":
                ax.text(
                    1.1 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "SH",
                    ha="right",
                    va="bottom",
                )
            if type_str == "surf" and comp == "BHZ":
                ax.text(
                    1.2 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "Z",
                    ha="right",
                    va="bottom",
                )
            if type_str == "surf" and comp == "SH":
                ax.text(
                    1.2 * min_time,
                    0.2 * max(abs(min_val), max_val),
                    "T",
                    ha="right",
                    va="bottom",
                )
        if custom == "syn":
            max_val = np.maximum(abs(min(waveform)), max(waveform))
            tmin, tmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            if type_str == "body":
                ax.text(
                    tmax,
                    yrel_annotation_max * ymin,
                    "{:0.1f}".format(max_val),
                    ha="right",
                    va="center",
                    color=color,
                )
            elif type_str == "surf":
                ax.text(
                    tmax,
                    yrel_annotation_max * ymin,
                    "{:0.3f}".format(max_val),
                    ha="right",
                    va="center",
                    color=color,
                )
            elif type_str == "cgps" or type_str == "strong":
                ax.text(
                    tmax,
                    yrel_annotation_max * ymin,
                    "{:0.2f}".format(max_val),
                    ha="right",
                    va="center",
                    color=color,
                )
        if type_str == "body":
            ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5, min_n_ticks=5))
            ax.yaxis.set_major_locator(ticker.NullLocator())
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        elif type_str == "surf":
            ax.xaxis.set_major_locator(ticker.MultipleLocator(500))
            ax.yaxis.set_major_locator(ticker.NullLocator())
        elif type_str == "cgps" or type_str == "strong":
            ax.xaxis.set_major_locator(ticker.MultipleLocator(40))
            ax.yaxis.get_major_locator().set_params(integer=True)  # type:ignore
        if nPlot > len(weights) - 3:
            ax.set_xlabel("Time After OT (s)")
        ax.grid(axis="x", which="both", linestyle="dotted", color="0.5")
    return axes


def plot_waveform_fits(
    files: List[dict],
    components: list,
    type_str: str,
    start_margin: int = 10,
    plot_directory: Union[pathlib.Path, str] = pathlib.Path(),
    additional_files: List[dict] = [],
):
    """Plot fit of observed to synthetic data for selected channels

    :param files: Waveform files to plot
    :type files: List[dict]
    :param components: Components (channels) of data selected for plotting
    :type components: list
    :param type_str: Data type of given waveform files
    :type type_str: str
    :param start_margin: Start margin of data for plotting, defaults to 10
    :type start_margin: int, optional
    :param plot_directory: Directory where plots should be saved, defaults to pathlib.Path()
    :type plot_directory: Union[pathlib.Path, str], optional
    """
    plot_directory = pathlib.Path(plot_directory)
    if type_str == "body" or type_str == "surf":
        files = [file for file in files if file["component"] in components]
        print("Creating Waveform Fit Plot: " + str(type_str) + " " + str(components[0]))
    if type_str == "cgps" or type_str == "strong":
        files = [file for file in files]
        print("Creating Waveform Fit Plot: " + str(type_str))
    if not len(files):
        print(
            f"No files found matching components ({components}). "
            "Not plotting waveform fits."
        )
        return
    files = sorted(files, key=lambda k: (k["azimuth"], k["component"]))
    if additional_files:
        if type_str == "body" or type_str == "surf":
            additional_files = [
                file for file in additional_files if file["component"] in components
            ]

        additional_files = sorted(
            additional_files, key=lambda k: (k["azimuth"], k["component"])
        )
        additional_syn_waveforms = [file["synthetic"] for file in additional_files]

    sampling = [file["dt"] for file in files]
    names = [file["name"] for file in files]
    azimuths = [file["azimuth"] for file in files]
    distances = [file["distance"] for file in files]
    weights = [file["trace_weight"] for file in files]
    obs_waveforms = [file["observed"] for file in files]
    syn_waveforms = [file["synthetic"] for file in files]
    comp = [file["component"] for file in files]
    zipped = zip(sampling, syn_waveforms)
    syn_times = [dt * np.arange(0, len(synthetic)) for dt, synthetic in zipped]
    start_waveform: List[float] = []
    for file in files:
        dt = file["dt"]
        nstart = file["start_signal"]
        margin = int(start_margin / dt)
        margin = min(nstart, margin)
        start_waveform = start_waveform + [margin]
    obs_times = [
        dt * np.arange(-start, len(observed) - start)
        for dt, start, observed in zip(sampling, start_waveform, obs_waveforms)
    ]

    if additional_files:
        shift_additional_file = [
            f1["start_signal"] - f2["start_signal"]
            for f1, f2 in zip(files, additional_files)
        ]
        # syn_times need to be udpate to account for potential different shift match
        additional_syn_times = [
            arr - dt * shift_additional_file[k] for k, arr in enumerate(syn_times)
        ]

    numrows_phase = len(files) // 3 + 1
    fig, axes = plt.subplots(
        max(3, numrows_phase), 3, figsize=(20, int(2.2 * numrows_phase))
    )
    axes2: List[Axes] = axes.ravel()  # type: ignore
    for ax in axes2[len(files) :]:
        ax.axis("off")

    if type_str == "body" or type_str == "surf":
        axes2 = plot_waveforms(
            list(axes2),
            obs_times,
            obs_waveforms,
            weights,
            type_str=type_str,
            comp=comp[0],
            color="black",
            custom="fill",
        )
        axes2 = plot_waveforms(
            list(axes2),
            syn_times,
            syn_waveforms,
            weights,
            type_str=type_str,
            comp=comp[0],
            color="red",
            custom="syn",
            yrel_annotation_max=0.5,
        )
        if additional_files:
            axes2 = plot_waveforms(
                list(axes2),
                additional_syn_times,
                additional_syn_waveforms,
                weights,
                type_str=type_str,
                comp=comp[0],
                color="blue",
                custom="syn",
                yrel_annotation_max=0.8,
            )
        dict = {
            "weights": weights,
            "azimuths": azimuths,
            "names": names,
            "distances": distances,
            "type_str": type_str,
            "comps": comp[0],
        }
    if type_str == "cgps" or type_str == "strong":
        axes2 = plot_waveforms(
            list(axes2),
            obs_times,
            obs_waveforms,
            weights,
            type_str=type_str,
            comp=comp,
            color="black",
            custom="fill",
        )
        axes2 = plot_waveforms(
            axes2,
            syn_times,
            syn_waveforms,
            weights,
            type_str=type_str,
            comp=comp,
            color="red",
            custom="syn",
            yrel_annotation_max=0.5,
        )
        if additional_files:
            axes2 = plot_waveforms(
                axes2,
                additional_syn_times,
                additional_syn_waveforms,
                weights,
                type_str=type_str,
                comp=comp,
                color="blue",
                custom="syn",
                yrel_annotation_max=0.8,
            )

        dict = {
            "weights": weights,
            "azimuths": azimuths,
            "names": names,
            "distances": distances,
            "type_str": type_str,
            "comps": comp,
        }
    axes2 = add_metadata(list(axes2), **dict)

    if type_str == "body":
        if "BHZ" in components:
            plot_name = "P_body_waves"
        if "SH" in components:
            plot_name = "SH_body_waves"

    if type_str == "surf":
        if "BHZ" in components:
            plot_name = "Rayleigh_surf_waves"
        if "SH" in components:
            plot_name = "Love_surf_waves"

    if type_str == "cgps":
        plot_name = "cGPS_waves"

    if type_str == "strong":
        plot_name = "strong_motion_waves"

    plt.savefig(plot_directory / (plot_name + ".png"), dpi=300, bbox_inches='tight')
    plt.savefig(plot_directory / (plot_name + ".ps"))
    plt.close()
    return


def retrieve_addition_traces(directory, data_type, stations=None):
    if data_type == "body":
        config_file = "tele_waves.json"
        syn_file = "synthetics_body.txt"
    elif data_type == "surf":
        config_file = "surf_waves.json"
        syn_file = "synthetics_surf.txt"
    elif data_type == "strong":
        config_file = "strong_motion_waves.json"
        syn_file = "synthetics_strong.txt"
    additional_traces_info = []
    with open(directory / config_file) as t:
        additional_traces_info = json.load(t)
    additional_traces_info = get_outputs.get_data_dict(
        additional_traces_info,
        syn_file=syn_file,
        directory=directory,
    )
    if stations:
        additional_traces_info = [
            file for file in additional_traces_info if file["name"] in stations
        ]
    return additional_traces_info


def plot_misfit(
    used_data_type: List[str],
    directories: Union[pathlib.Path, str] = pathlib.Path(),
    stations=None,
):
    """Plot misfit of observed and synthetic data

    :param used_data_type: The list of data types used
    :type used_data_type: List[str]
    :param directory: The location of data files, defaults to pathlib.Path()
    :type directory: Union[pathlib.Path, str], optional
    :raises FileNotFoundError: When a data type's json file is not found
    """
    directories = [pathlib.Path(directory) for directory in directories]
    directory = directories[0]

    if "body" in used_data_type:
        if not os.path.isfile(directory / "tele_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "tele_waves.json"
            )
        with open(directory / "tele_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_body.txt", directory=directory
        )
        if stations:
            traces_info = [file for file in traces_info if file["name"] in stations]

        additional_traces_info = (
            retrieve_addition_traces(directories[1], "body", stations)
            if len(directories)>1
            else []
        )

        values = [["BHZ"], ["SH"]]
        for components in values:
            plot_waveform_fits(
                traces_info,
                components,
                "body",
                plot_directory=directory,
                additional_files=additional_traces_info,
            )

    if "surf" in used_data_type:
        if not os.path.isfile(directory / "surf_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "surf_waves.json"
            )
        with open(directory / "surf_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_surf.txt", margin=0, directory=directory
        )
        if stations:
            traces_info = [file for file in traces_info if file["name"] in stations]

        additional_traces_info = (
            retrieve_addition_traces(directories[1], "surf", stations)
            if len(directories)>1
            else []
        )
        values = [["BHZ"], ["SH"]]
        for components in values:
            plot_waveform_fits(
                traces_info,
                components,
                "surf",
                plot_directory=directory,
                additional_files=additional_traces_info,
            )
    if "strong" in used_data_type:
        if not os.path.isfile(directory / "strong_motion_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "strong_motion_waves.json"
            )
        with open(directory / "strong_motion_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_strong.txt", directory=directory
        )
        if stations:
            traces_info = [file for file in traces_info if file["name"] in stations]
        additional_traces_info = (
            retrieve_addition_traces(directories[1], "strong", stations)
            if len(directories)
            else []
        )
        values = [["HLZ", "HNZ"], ["HLE", "HNE"], ["HLN", "HNN"]]
        plot_waveform_fits(
            traces_info,
            values,
            "strong",
            start_margin=10,
            plot_directory=directory,
            additional_files=additional_traces_info,
        )
    if "cgps" in used_data_type:
        if not os.path.isfile(directory / "cgps_waves.json"):
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), "cgps_waves.json"
            )
        with open(directory / "cgps_waves.json") as t:
            traces_info = json.load(t)
        traces_info = get_outputs.get_data_dict(
            traces_info, syn_file="synthetics_cgps.txt", directory=directory
        )
        values = [["LXZ", "LHZ", "LYZ"], ["LXE", "LHE", "LYE"], ["LXN", "LHN", "LYN"]]
        plot_waveform_fits(
            traces_info, values, "cgps", start_margin=10, plot_directory=directory
        )
    return


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="compare multiple WASP inversion results"
    )

    parser.add_argument(
        "input_folders",
        help="inversion folder separated by ','",
    )
    parser.add_argument(
        "waveform_types",
        help="waveform types (e.g. body,surf,strong) , separated by ','",
    )

    parser.add_argument(
        "--stations_body",
        nargs=1,
        help="stations for body waves. n values separated by ','",
    )
    parser.add_argument(
        "--stations_surf",
        nargs=1,
        help="stations for surface waves. n values separated by ','",
    )
    parser.add_argument(
        "--stations_strong",
        nargs=1,
        help="stations for strong waves. n values separated by ','",
    )
    args = parser.parse_args()

    directories = args.input_folders.split(",")
    waveform_types = args.waveform_types.split(",")

    if "body" in waveform_types:
        stations = None
        if args.stations_body:
            stations = args.stations_body[0].split(",")
        plot_misfit(["body"], directories, stations)

    if "surf" in waveform_types:
        stations = None
        if args.stations_surf:
            stations = args.stations_surf[0].split(",")
        plot_misfit(["surf"], directories, stations)

    if "strong" in waveform_types:
        stations = None
        if args.stations_strong:
            stations = args.stations_strong[0].split(",")
        plot_misfit(["strong"], directories, stations)
