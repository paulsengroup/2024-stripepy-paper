# Copyright (C) 2024 Andrea Raffo <andrea.raffo@ibv.uio.no>
#
# SPDX-License-Identifier: MIT


import itertools
import os
import pathlib
from typing import Dict, List, Tuple, Union

import bioframe as bf
import h5py
import hictkpy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib.ticker import EngFormatter, ScalarFormatter
from numpy.typing import NDArray

# Colors
fruit_punch = sns.blend_palette(["white", "red"], as_cmap=True)
colors = ["#f4a261", "#e9c46a", "#2a9d8f", "#ae6c96"]


def _raise_invalid_bin_type_except(f: hictkpy.File):
    raise RuntimeError(f"Only files with a uniform bin size are supported, found \"{f.attributes()['bin-type']}\".")


def cmap_loading(path: os.PathLike, resolution: int):
    try:
        if not isinstance(resolution, int):
            raise TypeError("resolution must be an integer.")

        if resolution <= 0:
            raise ValueError("resolution must be greater than zero.")

        if hictkpy.is_scool_file(path):
            raise RuntimeError(".scool files are not currently supported.")
        if hictkpy.is_cooler(path):
            f = hictkpy.File(path)
            if f.resolution() == 0:
                _raise_invalid_bin_type_except(f)
            if f.resolution() != resolution:
                raise RuntimeError(f"expected {resolution} resolution, found {f.resolution()}.")
        else:
            f = hictkpy.MultiResFile(path)[resolution]
    except RuntimeError as e:
        raise RuntimeError(f'error opening file "{path}"') from e

    if f.attributes().get("bin-type", "fixed") != "fixed":
        _raise_invalid_bin_type_except(f)

    # Retrieve metadata:
    chr_starts = [0]  # left ends of each chromosome  (in matrix coordinates)
    chr_ends = []  # right ends of each chromosome (in matrix coordinates)
    chr_sizes = []  # integer bp lengths, one per chromosome
    for bp_length in f.chromosomes().values():
        chr_sizes.append(bp_length)
        chr_ends.append(chr_starts[-1] + int(np.ceil(bp_length / resolution)))
        chr_starts.append(chr_ends[-1])
        # print(f"{chr_starts[-2]}-{chr_ends[-1]}-{chr_sizes[-1]}")
    chr_starts.pop(-1)

    return f, chr_starts, chr_ends, chr_sizes


def chromosomes_to_study(chromosomes, length_in_bp, min_size_allowed):
    # Extract the list of chromosomes:
    chr_ids = list(range(len(chromosomes)))
    c_pairs = list(zip(chr_ids, chromosomes))

    # Remove overly-short chromosomes:
    surviving_indices = [i for i, e in enumerate(length_in_bp) if e > min_size_allowed]
    deleted_indices = [i for i, e in enumerate(length_in_bp) if e <= min_size_allowed]
    if len(deleted_indices) > 0:
        print(
            f"ATT: The following chromosomes are discarded because shorter than MIN_SIZE_CHROMOSOME = "
            f"{min_size_allowed} bp: {[chromosomes[i] for i in deleted_indices]}"
        )
        c_pairs = [c_pairs[i] for i in surviving_indices]

        # If there is no chromosome left, exit:
        if len(c_pairs) == 0:
            raise ValueError(f"\nNo chromosome is long enough... decrease the parameter MIN_SIZE_CHROMOSOME")

    return c_pairs


def retrieve_ground_truth(
    path: pathlib.Path,
    file_name: str,
    chromosome: str,
    resolution: int,
    probability_threshold: float,
    number_of_bins: int,
) -> Dict[str, Dict[str, NDArray]]:
    """
    This function loads ground truth (GT) information and stores it in a dictionary.

    Parameters
    ----------
    path: pathlib.Path
        path to the ground truth file
    file_name: str
        name of the bed file containing coordinates and occupancy of each extrusion barrier
    chromosome: str
        chromosome name
    resolution: int
        resolution of the contact map
    probability_threshold: float
        cut-off value for the occupancy
    number_of_bins: int
        number of bins for current chromosome at given resolution -- used to define a classification vector

    Returns
    -------
    Dict[str, Dict[str, NDArray]]
        a dictionary for the ground truth annotations containing:

         * location of the anchor points
         * classification vectors w.r.t. the genomic bins (0: bin does not host anchor, 1: bin hosts anchor)
    """

    # Ground truth for the chromosome under study:
    df = bf.read_table(path / f"{file_name}", schema="bed6")
    df_sub = bf.select(df, chromosome).copy()
    df_sub["anchor"] = (df_sub["end"] + df_sub["start"]) / 2

    # Lower ground truth stripes:
    LT_stripes = df_sub[(df_sub.score > probability_threshold) & (df_sub.strand == "+")]

    # Upper ground truth stripes:
    UT_stripes = df_sub[(df_sub.score > probability_threshold) & (df_sub.strand == "-")]

    # Retrieve ground truth anchors
    GT = {"anchors": {}, "clas_vec": {}}
    GT["anchors"]["LT"] = np.unique(np.round(LT_stripes.anchor.values / resolution)).astype(int)
    GT["anchors"]["UT"] = np.unique(np.round(UT_stripes.anchor.values / resolution)).astype(int)

    # Classification vectors:
    GT["clas_vec"]["LT"] = np.isin(np.arange(number_of_bins), GT["anchors"]["LT"]).astype(int)
    GT["clas_vec"]["UT"] = np.isin(np.arange(number_of_bins), GT["anchors"]["UT"]).astype(int)

    return GT


def retrieve_stripepy(
    path: pathlib.Path, chromosome: str, number_of_bins: int, threshold: float
) -> Tuple[Union[Dict[str, NDArray], None]]:
    """
    This function loads predictions produced by StripePy and stores them in a dictionary.

    Parameters
    ----------
    path: pathlib.Path
        path to the file produced by StripePy
    chromosome: str
        chromosome name
    number_of_bins: int
        number of bins for current chromosome at given resolution -- used to define a classification vector
    threshold: float
        cut-off value for the relative change parameter

    Returns
    -------
    Tuple[Union[Dict[str, NDArray], None]]
        four dictionaries containing:

         * coordinates of the horizontal domains
         * coordinates of the vertical domains
         * predicted classification vectors w.r.t. the genomic bins (0: bin does not host seed, 1: bin hosts seed)
         * seed sites
    """

    # Open hdf5 file:
    hf = h5py.File(str(path / "results.hdf5"), "r")

    # If the current chromosome is in the HDF5 file...:
    if chromosome in hf:

        # Initialize dictionaries for relevant information:
        HIoIs = {"LT": None, "UT": None}
        VIoIs = {"LT": None, "UT": None}
        seeds = {"LT": None, "UT": None}
        bio_descriptors = {"LT": None, "UT": None}
        geo_descriptors = {"LT": None, "UT": None}
        candidates2keep = {"LT": None, "UT": None}
        clas_vec = {"LT": None, "UT": None}

        for w in ["LT", "UT"]:
            bio_descriptors[w] = np.array(hf[f"{chromosome}/stripes/{w}/bio-descriptors"])
            geo_descriptors[w] = np.array(hf[f"{chromosome}/stripes/{w}/geo-descriptors"])

            candidates2keep[w] = np.where(bio_descriptors[w][:, 2] > threshold)[0]
            HIoIs[w] = geo_descriptors[w][candidates2keep[w], 2:4].astype(int)
            VIoIs[w] = geo_descriptors[w][candidates2keep[w], 4:6].astype(int)
            seeds[w] = geo_descriptors[w][candidates2keep[w], 0].astype(int)
            clas_vec[w] = np.where(np.isin(range(number_of_bins), seeds[w]), 1, 0)

    else:
        # Group does not exist
        print("Group does not exist")
        HIoIs, VIoIs, clas_vec, seeds, top_pers = None, None, None, None, None

    return HIoIs, VIoIs, clas_vec, seeds


def retrieve_chromosight(
    path: pathlib.Path, chromosome: str, number_of_bins: int, resolution: int
) -> Tuple[Dict[str, NDArray]]:
    """
    This function loads predictions produced by Chromosight and stores them in a dictionary.

    Parameters
    ----------
    path: pathlib.Path
        path to the file produced by Chromosight
    chromosome: str
        chromosome name
    number_of_bins: int
        number of bins for current chromosome at given resolution -- used to define a classification vector
    resolution: int
        resolution of the contact map

    Returns
    -------
    Tuple[Dict[str: NDArray]]
        four dictionaries containing:

         * coordinates of the horizontal domains
         * coordinates of the vertical domains
         * predicted classification vectors w.r.t. the genomic bins (0: bin does not host anchor, 1: bin hosts anchor)
         * predicted anchor sites
    """

    # Load predictions -- lower-triangular:
    Ldf = pd.read_csv(str(path / "left" / "output.tsv"), sep="\t")
    Ldf_chr = Ldf[Ldf["chrom1"] == chromosome]

    # Load predictions -- upper-triangular:
    Udf = pd.read_csv(str(path / "right" / "output.tsv"), sep="\t")
    Udf_chr = Udf[Udf["chrom1"] == chromosome]

    # Gather lower- and upper-triangular candidates:
    # NB: Chromosight does not estimate the width of a stripe, but constrained it to 1 bin.
    HIoIs = {
        "LT": (Ldf_chr[["start1", "start1"]].values / resolution).astype(int),
        "UT": (Udf_chr[["start2", "start2"]].values / resolution).astype(int),
    }
    VIoIs = {
        "LT": (Ldf_chr[["start2", "start2"]].values / resolution).astype(int),
        "UT": (Udf_chr[["start1", "start1"]].values / resolution).astype(int),
    }

    # Gather anchors:
    anchors = {"LT": np.array([x[0] for x in HIoIs["LT"]]), "UT": np.array([x[0] for x in HIoIs["UT"]])}

    # Classification vectors (predicted):
    clas_vec = {
        "LT": np.where(np.isin(range(number_of_bins), anchors["LT"]), 1, 0),
        "UT": np.where(np.isin(range(number_of_bins), anchors["UT"]), 1, 0),
    }

    return HIoIs, VIoIs, clas_vec, anchors


def retrieve_stripecaller(
    path: pathlib.Path, chromosome: str, number_of_bins: int, resolution: int
) -> Tuple[Dict[str, NDArray]]:
    """
    This function loads predictions produced by StripeCaller and stores them in a dictionary.

    Parameters
    ----------
    path: pathlib.Path
        path to the file produced by StripeCaller
    chromosome: str
        chromosome name
    number_of_bins: int
        number of bins for current chromosome at given resolution -- used to define a classification vector
    resolution: int
        resolution of the contact map

    Returns
    -------
    Tuple[Dict[str: NDArray]]
        four dictionaries containing:

         * coordinates of the horizontal domains
         * coordinates of the vertical domains
         * predicted classification vectors w.r.t. the genomic bins (0: bin does not host anchor, 1: bin hosts anchor)
         * predicted anchor sites
    """

    # Load predictions:
    df = bf.read_table(str(path / "output.bedpe"), schema="bed6")
    df_chr = df[df["chrom"] == chromosome]
    X1 = df_chr["start"].values.tolist()
    X2 = df_chr["end"].values.tolist()
    Y1 = df_chr["score"].values.tolist()
    Y2 = df_chr["strand"].values.tolist()

    # Gather lower- and upper-triangular candidates:
    # NB: StripeCaller does not estimate the width of a stripe, but constrained it to 1 bin.
    L_ids = [i for i, (x1, x2, y1, y2) in enumerate(zip(X1, X2, Y1, Y2)) if x2 - x1 < y2 - y1]
    U_ids = [i for i, (x1, x2, y1, y2) in enumerate(zip(X1, X2, Y1, Y2)) if y2 - y1 < x2 - x1]
    HIoIs = {
        "LT": np.array([[s / resolution, s / resolution] for idx, s in enumerate(X1) if idx in L_ids]).astype(int),
        "UT": np.array([[s / resolution, s / resolution] for idx, s in enumerate(Y1) if idx in U_ids]).astype(int),
    }
    VIoIs = {
        "LT": np.array(
            [[s / resolution, e / resolution] for idx, (s, e) in enumerate(zip(Y1, Y2)) if idx in L_ids]
        ).astype(int),
        "UT": np.array(
            [[s / resolution, e / resolution] for idx, (s, e) in enumerate(zip(X1, X2)) if idx in U_ids]
        ).astype(int),
    }

    # Gather anchors:
    anchors = {
        "LT": np.array([x[0] for x in HIoIs["LT"]]).astype(int),
        "UT": np.array([x[0] for x in HIoIs["UT"]]).astype(int),
    }

    # Classification vectors (predicted):
    clas_vec = {
        "LT": np.array([1 if i in anchors["LT"] else 0 for i in range(number_of_bins)]).astype(int),
        "UT": np.array([1 if i in anchors["UT"] else 0 for i in range(number_of_bins)]).astype(int),
    }

    return HIoIs, VIoIs, clas_vec, anchors


def retrieve_stripenn(
    path: pathlib.Path, chromosome: str, number_of_bins: int, resolution: int, filter_: bool = False
) -> Tuple[Dict[str, NDArray]]:
    """
    This function loads predictions produced by Stripenn and stores them in a dictionary.

    Parameters
    ----------
    path: pathlib.Path
        path to the file produced by Stripenn
    chromosome: str
        chromosome name
    number_of_bins: int
        number of bins for current chromosome at given resolution -- used to define a classification vector
    resolution: int
        resolution of the contact map
    filter_: bool
        if True, filtered predictions are loaded; if False, unfiltered predictions are loaded (default: False)

    Returns
    -------
    Tuple[Dict[str: NDArray]]
        four dictionaries containing:

         * coordinates of the horizontal domains
         * coordinates of the vertical domains
         * predicted classification vectors w.r.t. the genomic bins (0: bin does not host anchor, 1: bin hosts anchor)
         * predicted anchor sites
    """

    # Retrieve dataframe:
    if filter_ is False:
        df = pd.read_csv(str(path / "result_unfiltered.tsv"), sep="\t")
    else:
        df = pd.read_csv(str(path / "result_filtered.tsv"), sep="\t")

    # Filter rows for the specified chromosome and candidate types (lower- and upper-triangular candidates):
    L_df_chr = df[(df["chr"] == chromosome) & (df["pos1"] == df["pos3"])]
    U_df_chr = df[(df["chr"] == chromosome) & (df["pos2"] == df["pos4"])]

    # Gather lower-triangular candidates:
    HIoIs = {
        "LT": (L_df_chr[["pos1", "pos2"]].values / resolution).astype(int),
        "UT": (U_df_chr[["pos1", "pos2"]].values / resolution).astype(int),
    }
    VIoIs = {
        "LT": (L_df_chr[["pos3", "pos4"]].values / resolution).astype(int),
        "UT": (U_df_chr[["pos3", "pos4"]].values / resolution).astype(int),
    }

    # Gather anchors:
    anchors = {
        "LT": np.array([int((x[0] + x[1]) / 2) for x in HIoIs["LT"]]).astype(int),
        "UT": np.array([int((x[0] + x[1]) / 2) for x in HIoIs["UT"]]).astype(int),
    }

    # Classification vectors (predicted):
    clas_vec = {
        "LT": np.array([1 if i in anchors["LT"] else 0 for i in range(number_of_bins)]).astype(int),
        "UT": np.array([1 if i in anchors["UT"] else 0 for i in range(number_of_bins)]).astype(int),
    }

    return HIoIs, VIoIs, clas_vec, anchors


def format_ticks(ax, x=True, y=True, rotate=True):

    if x:
        ax.xaxis.set_major_formatter(EngFormatter("b"))
        ax.xaxis.tick_bottom()
    else:
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.tick_bottom()

    if y:
        ax.yaxis.set_major_formatter(EngFormatter("b"))
    else:
        ax.yaxis.set_major_formatter(ScalarFormatter())

    if rotate:
        ax.tick_params(axis="x", rotation=45)


def StripeBench_csv_counts(
    results: pd.DataFrame,
    resolutions: List[int],
    contact_densities: List[int],
    noises: List[int],
    output_path: pathlib.Path,
):
    """
    This function produces Supplementary Table 1 in the manuscript.

    Parameters
    ----------
    results: pd.DataFrame
        dataframe containing measures and other quantities for the callers under analysis (Stripebench benchmark)
    resolutions: List[int]
        list of resolutions in the Stripebench benchmark
    contact_densities: List[int]
        list of contact densities in the Stripebench benchmark
    noises: List[int]
        list of noise levels in the Stripebench benchmark
    output_path: pathlib.Path
        path to output folder
    """

    # List of lists, that will contain each row of the dataframe:
    rows = []

    # Loop over resolutions:
    for resolution in resolutions:
        n_res = -1

        # Loop over contact densities:
        for contact in contact_densities:
            n_res += 1

            # Initialize list of entries for current row:
            row = []

            # Enter resolution (in kb) and contact density:
            row.append(f"{round(resolution / 1000)}kb")
            row.append(f"{round(contact)}")

            # Loop over levels of noise:
            for noise in noises:

                # Slice the dataframe:
                sliced_results = results[
                    (results["Resolution"] == resolution)
                    & (results["Noise"] == noise)
                    & (results["Contact Density"] == contact)
                ]

                for m in ["stripepy", "chromosight", "stripecaller", "stripenn"]:
                    num_stripes = np.sum(
                        sliced_results[sliced_results["Method"] == m]["classification_vector"].values[0]
                    )
                    row.append(num_stripes)

            rows.append(row)

    # Store results in a pandas dataframe:
    df = pd.DataFrame(
        rows,
        columns=[
            "rho",
            "delta",
            "M1",
            "M2",
            "M3",
            "M4",
            "M1",
            "M2",
            "M3",
            "M4",
            "M1",
            "M2",
            "M3",
            "M4",
            "M1",
            "M2",
            "M3",
            "M4",
        ],
    )

    df.to_csv(str(output_path / "tables" / "table1.csv"), index=False)
    print("Done.")


def StripeBench_csv_measures(
    results: pd.DataFrame,
    resolutions: List[int],
    contact_densities: List[int],
    noises: List[int],
    output_path: pathlib.Path,
):
    """
    This function produces Supplementary Tables 2-5 in the manuscript.

    Parameters
    ----------
    results: pd.DataFrame
        dataframe containing measures and other quantities for the callers under analysis (Stripebench benchmark)
    resolutions: List[int]
        list of resolutions in the Stripebench benchmark
    contact_densities: List[int]
        list of contact densities in the Stripebench benchmark
    noises: List[int]
        list of noise levels in the Stripebench benchmark
    output_path: pathlib.Path
        path to output folder
    """

    # Loop over resolutions:
    for n_res, resolution in enumerate(resolutions):

        # List of lists, that will contain each row of the dataframe:
        rows = []

        # Loop over levels of noise:
        for noise in noises:

            # Slice the dataframe:
            sliced_results = results[(results["Resolution"] == resolution) & (results["Noise"] == noise)]

            # Print:
            for num_meas, m in enumerate(
                ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]
            ):

                # Initialize list of entries for current row:
                row = []

                # Enter noise level and measure acronym:
                row.append(f"{round(noise / 1000)}k")
                row.append(m)

                for contact in contact_densities:
                    for method in ["stripepy", "chromosight", "stripecaller", "stripenn"]:
                        row.append(
                            f"{sliced_results[
                                          (sliced_results["Method"] == method) &
                                          (sliced_results["Contact Density"] == contact)][m].values[0] * 100:>6.2f}"
                        )

                rows.append(row)

        # Store results in a pandas dataframe:
        df = pd.DataFrame(
            rows,
            columns=[
                "sigma",
                "measure",
                "M1",
                "M2",
                "M3",
                "M4",
                "M1",
                "M2",
                "M3",
                "M4",
                "M1",
                "M2",
                "M3",
                "M4",
                "M1",
                "M2",
                "M3",
                "M4",
            ],
        )
        df.to_csv(str(output_path / "tables" / f"table{n_res+2}.csv"), index=False)
    print("Done.")


def marginal_plots(
    results: pd.DataFrame,
    resolutions: List[int],
    contact_densities: List[int],
    noises: List[int],
    output_path: pathlib.Path,
):
    """
    This function produces Figure 4 and Extended Data Figures 1-3 in the manuscript.

    Parameters
    ----------
    results: pd.DataFrame
        dataframe containing measures and other quantities for the callers under analysis (Stripebench benchmark)
    resolutions: List[int]
        list of resolutions in the Stripebench benchmark
    contact_densities: List[int]
        list of contact densities in the Stripebench benchmark
    noises: List[int]
        list of noise levels in the Stripebench benchmark
    output_path: pathlib.Path
        path to output folder
    """

    excluded_methods = ["stripenn"]
    results = results[~results["Method"].isin(excluded_methods)]

    # 1) CHANGE IN RESOLUTION

    # Storing the medians for plots!
    Q2s_by_res = dict()
    for meas in ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]:
        Q2s_by_res[meas] = []

    # Looping in resolution:
    fig, axes = plt.subplots(12, 4, figsize=(6, 11))
    for n_res, resolution in enumerate(resolutions):

        # Extract rows referring to current resolution:
        sliced_results = results.loc[results["Resolution"] == resolution]

        # Medians for current resolution:
        Q2s_this_res = dict()

        for n_meas, clas_meas_name in enumerate(
            ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]
        ):

            # Medians for current resolution and current measure:
            Q2s_this_res[clas_meas_name] = []
            for method in ["stripepy", "chromosight", "stripecaller"]:
                Q2s_this_res[clas_meas_name].append(
                    np.median(sliced_results.loc[sliced_results["Method"] == method][clas_meas_name].values)
                )
            Q2s_by_res[clas_meas_name].append(Q2s_this_res[clas_meas_name])

            # Box Plots
            ax = axes[n_meas, n_res]
            sns.boxplot(
                y=clas_meas_name,
                x="Method",
                data=sliced_results,
                ax=ax,
                palette=colors[0:3],
                hue="Method",
                orient="v",
                width=0.5,
                fliersize=3,
            )
            ax.set(xlabel=None, ylabel=None)
            ax.tick_params(labelbottom=False, bottom=False)
            ax.yaxis.set_tick_params(labelsize=7)
            ax.axes.get_xaxis().set_visible(False)
            if n_res == 0:
                ax.set_ylabel(clas_meas_name, fontsize=12)
            if n_meas == 0:
                ax.set_title(r"$\rho$" f" = {int(resolution / 1000)}kb", fontsize=12)
                ax.set_ylim((-0.02, 0.52))
                minor_ticks = np.linspace(0.00, 0.50, 11)
                major_ticks = np.linspace(0.00, 0.50, 6)
            elif n_meas == 1:
                ax.set_ylim((0.90, 1.005))
                minor_ticks = np.linspace(0.90, 1.00, 11)
                major_ticks = np.linspace(0.90, 1.00, 6)
            elif n_meas == 2:
                ax.set_ylim((-0.02, 1.02))
                minor_ticks = np.linspace(0.00, 1.00, 11)
                major_ticks = np.linspace(0.00, 1.00, 6)
            elif n_meas == 3:
                ax.set_ylim((0.49, 0.71))
                minor_ticks = np.linspace(0.50, 0.70, 11)
                major_ticks = np.linspace(0.50, 0.70, 6)
            elif n_meas == 4:
                ax.set_ylim((-0.05, 0.70))
                minor_ticks = np.linspace(0.00, 0.70, 11)
                major_ticks = np.linspace(0.00, 0.70, 6)
            elif n_meas == 5:
                ax.set_ylim((-0.05, 0.30))
                minor_ticks = np.linspace(0.00, 0.30, 11)
                major_ticks = np.linspace(0.00, 0.30, 6)
            elif n_meas == 6:
                ax.set_ylim((-0.05, 0.45))
                minor_ticks = np.linspace(0.00, 0.45, 11)
                major_ticks = np.linspace(0.00, 0.45, 6)
            elif n_meas == 7:
                ax.set_ylim((-0.05, 0.50))
                minor_ticks = np.linspace(0.00, 0.50, 11)
                major_ticks = np.linspace(0.00, 0.50, 6)
            elif n_meas == 8:
                ax.set_ylim((-0.05, 0.95))
                minor_ticks = np.linspace(0.00, 0.95, 11)
                major_ticks = np.linspace(0.00, 0.95, 6)
            elif n_meas == 9:
                ax.set_ylim((-0.05, 1.00))
                minor_ticks = np.linspace(0.00, 1.00, 11)
                major_ticks = np.linspace(0.00, 1.00, 6)
            elif n_meas == 10:
                ax.set_ylim((-0.05, 0.90))
                minor_ticks = np.linspace(0.00, 0.90, 11)
                major_ticks = np.linspace(0.00, 0.90, 6)
            elif n_meas == 11:
                ax.set_ylim((-0.05, 0.90))
                minor_ticks = np.linspace(0.00, 0.90, 11)
                major_ticks = np.linspace(0.00, 0.90, 6)
            ax.yaxis.set_major_locator(ticker.FixedLocator(major_ticks))
            ax.yaxis.set_minor_locator(ticker.FixedLocator(minor_ticks))
            ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
            ax.yaxis.grid(True, which="minor", linestyle="--", linewidth=0.5)

    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    plt.tight_layout()
    plt.savefig(str(output_path / "boxplots/bp_by_res.svg"), bbox_inches="tight")
    plt.clf()
    plt.close(fig)

    # Plots of medians:
    fig, axes = plt.subplots(3, 4, figsize=(10, 6))
    axes = axes.flatten()
    for ax, meas in zip(axes, ["TPR", "TNR", "PPV", "bACC", "GM", "F1c", "FMc", "JI", "AHR", "FGC", "F1r", "FMr"]):
        ax.plot([5, 10, 25, 50], [x for x, _, _ in Q2s_by_res[meas]], color=colors[0], linestyle="dashdot")
        ax.plot([5, 10, 25, 50], [x for _, x, _ in Q2s_by_res[meas]], color=colors[1], linestyle="dashdot")
        ax.plot([5, 10, 25, 50], [x for _, _, x in Q2s_by_res[meas]], color=colors[2], linestyle="dashdot")
        ax.plot([5, 10, 25, 50], [x for x, _, _ in Q2s_by_res[meas]], "o", color=colors[0])
        ax.plot([5, 10, 25, 50], [x for _, x, _ in Q2s_by_res[meas]], "o", color=colors[1])
        ax.plot([5, 10, 25, 50], [x for _, _, x in Q2s_by_res[meas]], "o", color=colors[2])
        ax.set_title(meas, fontsize=12)
        ax.grid(color="green", linestyle="--", linewidth=0.5, axis="y")
        ax.xaxis.set_major_locator(ticker.FixedLocator(np.array([5, 10, 25, 50])))
        ax.xaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    plt.suptitle("Change in resolution", fontsize=16)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    plt.tight_layout()
    plt.savefig(str(output_path / "medians/medians_by_res.svg"), bbox_inches="tight")
    plt.clf()
    plt.close(fig)

    # 2) CHANGE IN CONTACT DENSITY

    # Storing the medians for plots!
    Q2s_by_cd = dict()
    for meas in ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]:
        Q2s_by_cd[meas] = []

    # Looping in contact density:
    fig, axes = plt.subplots(12, 4, figsize=(6, 11))
    for n_cd, cd in enumerate(contact_densities):

        # Extract rows referring to current contact density:
        sliced_results = results.loc[results["Contact Density"] == cd]

        # Medians for current contact density:
        Q2s_this_cd = dict()

        for n_meas, clas_meas_name in enumerate(
            ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]
        ):

            # Medians for current resolution and current measure:
            Q2s_this_cd[clas_meas_name] = []
            for method in ["stripepy", "chromosight", "stripecaller"]:
                Q2s_this_cd[clas_meas_name].append(
                    np.median(sliced_results.loc[sliced_results["Method"] == method][clas_meas_name].values)
                )
            Q2s_by_cd[clas_meas_name].append(Q2s_this_cd[clas_meas_name])

            # Box Plots
            ax = axes[n_meas, n_cd]
            sns.boxplot(
                y=clas_meas_name,
                x="Method",
                data=sliced_results,
                ax=ax,
                palette=colors[0:3],
                hue="Method",
                orient="v",
                width=0.5,
                fliersize=3,
            )
            ax.set(xlabel=None, ylabel=None)
            ax.tick_params(labelbottom=False, bottom=False)
            ax.yaxis.set_tick_params(labelsize=7)
            ax.axes.get_xaxis().set_visible(False)
            if n_cd == 0:
                ax.set_ylabel(clas_meas_name, fontsize=12)
            if n_meas == 0:
                ax.set_title(r"$\delta$" f" = {int(cd)}", fontsize=12)
                ax.set_ylim((-0.02, 0.50))
                minor_ticks = np.linspace(0.00, 0.50, 11)
                major_ticks = np.linspace(0.00, 0.50, 6)
            elif n_meas == 1:
                ax.set_ylim((0.90, 1.005))
                minor_ticks = np.linspace(0.90, 1.00, 11)
                major_ticks = np.linspace(0.90, 1.00, 6)
            elif n_meas == 2:
                ax.set_ylim((-0.02, 1.03))
                minor_ticks = np.linspace(0.00, 1.00, 11)
                major_ticks = np.linspace(0.00, 1.00, 6)
            elif n_meas == 3:
                ax.set_ylim((0.49, 0.72))
                minor_ticks = np.linspace(0.50, 0.70, 11)
                major_ticks = np.linspace(0.50, 0.70, 6)
            elif n_meas == 4:
                ax.set_ylim((-0.02, 0.70))
                minor_ticks = np.linspace(0.00, 0.70, 11)
                major_ticks = np.linspace(0.00, 0.70, 6)
            elif n_meas == 5:
                ax.set_ylim((-0.02, 0.325))
                minor_ticks = np.linspace(0.00, 0.30, 11)
                major_ticks = np.linspace(0.00, 0.30, 6)
            elif n_meas == 6:
                ax.set_ylim((-0.03, 0.51))
                minor_ticks = np.linspace(0.00, 0.5, 11)
                major_ticks = np.linspace(0.00, 0.5, 6)
            elif n_meas == 7:
                ax.set_ylim((-0.05, 0.5))
                minor_ticks = np.linspace(0.00, 0.5, 11)
                major_ticks = np.linspace(0.00, 0.5, 6)
            elif n_meas == 8:
                ax.set_ylim((-0.05, 1.00))
                minor_ticks = np.linspace(0.00, 1.00, 11)
                major_ticks = np.linspace(0.00, 1.00, 6)
            elif n_meas == 9:
                ax.set_ylim((-0.05, 1.03))
                minor_ticks = np.linspace(0.00, 1.00, 11)
                major_ticks = np.linspace(0.00, 1.00, 6)
            elif n_meas == 10:
                ax.set_ylim((-0.05, 0.90))
                minor_ticks = np.linspace(0.00, 0.90, 11)
                major_ticks = np.linspace(0.00, 0.90, 6)
            elif n_meas == 11:
                ax.set_ylim((-0.05, 0.90))
                minor_ticks = np.linspace(0.00, 0.90, 11)
                major_ticks = np.linspace(0.00, 0.90, 6)
            ax.yaxis.set_major_locator(ticker.FixedLocator(major_ticks))
            ax.yaxis.set_minor_locator(ticker.FixedLocator(minor_ticks))
            ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
            ax.yaxis.grid(True, which="minor", linestyle="--", linewidth=0.5)

    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    plt.tight_layout()
    plt.savefig(str(output_path / "boxplots/bp_by_cd.svg"), bbox_inches="tight")
    plt.clf()
    plt.close(fig)

    # Plots of medians:
    fig, axes = plt.subplots(3, 4, figsize=(10, 6))
    axes = axes.flatten()
    for ax, meas in zip(axes, ["TPR", "TNR", "PPV", "bACC", "GM", "F1c", "FMc", "JI", "AHR", "FGC", "F1r", "FMr"]):
        ax.plot([1, 5, 10, 15], [x for x, _, _ in Q2s_by_cd[meas]], color=colors[0], linestyle="dashdot")
        ax.plot([1, 5, 10, 15], [x for _, x, _ in Q2s_by_cd[meas]], color=colors[1], linestyle="dashdot")
        ax.plot([1, 5, 10, 15], [x for _, _, x in Q2s_by_cd[meas]], color=colors[2], linestyle="dashdot")
        ax.plot([1, 5, 10, 15], [x for x, _, _ in Q2s_by_cd[meas]], "o", color=colors[0])
        ax.plot([1, 5, 10, 15], [x for _, x, _ in Q2s_by_cd[meas]], "o", color=colors[1])
        ax.plot([1, 5, 10, 15], [x for _, _, x in Q2s_by_cd[meas]], "o", color=colors[2])
        ax.set_title(meas, fontsize=12)
        ax.grid(color="green", linestyle="--", linewidth=0.5, axis="y")
        ax.xaxis.set_major_locator(ticker.FixedLocator(np.array([1, 5, 10, 15])))
        ax.xaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    plt.suptitle("Change in contact density", fontsize=16)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    plt.tight_layout()
    plt.savefig(str(output_path / "medians/medians_by_cd.svg"), bbox_inches="tight")
    plt.clf()
    plt.close(fig)

    # 3) CHANGE IN NOISE LEVELS

    # Storing the medians for plots!
    Q2s_by_ns = dict()
    for meas in ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]:
        Q2s_by_ns[meas] = []

    # Looping in noise levels:
    fig, axes = plt.subplots(12, 4, figsize=(6, 11))

    for n_ns, noise in enumerate(noises):

        # Extract rows referring to current noise:
        sliced_results = results.loc[results["Noise"] == noise]

        # Medians for current contact density:
        Q2s_this_ns = dict()

        for n_meas, clas_meas_name in enumerate(
            ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]
        ):

            # Medians for current resolution and current measure:
            Q2s_this_ns[clas_meas_name] = []
            for method in ["stripepy", "chromosight", "stripecaller"]:
                Q2s_this_ns[clas_meas_name].append(
                    np.median(sliced_results.loc[sliced_results["Method"] == method][clas_meas_name].values)
                )
            Q2s_by_ns[clas_meas_name].append(Q2s_this_ns[clas_meas_name])

            # Box Plots
            ax = axes[n_meas, n_ns]
            sns.boxplot(
                y=clas_meas_name,
                x="Method",
                data=sliced_results,
                ax=ax,
                palette=colors[0:3],
                hue="Method",
                orient="v",
                width=0.5,
                fliersize=3,
            )
            ax.set(xlabel=None, ylabel=None)
            ax.tick_params(labelbottom=False, bottom=False)
            ax.yaxis.set_tick_params(labelsize=7)
            ax.axes.get_xaxis().set_visible(False)
            if n_ns == 0:
                ax.set_ylabel(clas_meas_name, fontsize=12)
            if n_meas == 0:
                ax.set_title(r"$\sigma$" f" = {int(noise / 1000)}k", fontsize=12)
                ax.set_ylim((-0.02, 0.500))
                minor_ticks = np.linspace(0.00, 0.500, 11)
                major_ticks = np.linspace(0.00, 0.500, 6)
            elif n_meas == 1:
                ax.set_ylim((0.90, 1.005))
                minor_ticks = np.linspace(0.90, 1.00, 11)
                major_ticks = np.linspace(0.90, 1.00, 6)
            elif n_meas == 2:
                ax.set_ylim((-0.02, 1.02))
                minor_ticks = np.linspace(0.00, 1.00, 11)
                major_ticks = np.linspace(0.00, 1.00, 6)
            elif n_meas == 3:
                ax.set_ylim((0.49, 0.72))
                minor_ticks = np.linspace(0.50, 0.70, 11)
                major_ticks = np.linspace(0.50, 0.70, 6)
            elif n_meas == 4:
                ax.set_ylim((-0.02, 0.70))
                minor_ticks = np.linspace(0.00, 0.70, 11)
                major_ticks = np.linspace(0.00, 0.70, 6)
            elif n_meas == 5:
                ax.set_ylim((-0.02, 0.325))
                minor_ticks = np.linspace(0.00, 0.3, 11)
                major_ticks = np.linspace(0.00, 0.3, 6)
            elif n_meas == 6:
                ax.set_ylim((-0.03, 0.5))
                minor_ticks = np.linspace(0.00, 0.5, 11)
                major_ticks = np.linspace(0.00, 0.5, 6)
            elif n_meas == 7:
                ax.set_ylim((-0.05, 0.5))
                minor_ticks = np.linspace(0.00, 0.5, 11)
                major_ticks = np.linspace(0.00, 0.5, 6)
            elif n_meas == 8:
                ax.set_ylim((-0.05, 0.975))
                minor_ticks = np.linspace(0.00, 0.975, 11)
                major_ticks = np.linspace(0.00, 0.975, 6)
            elif n_meas == 9:
                ax.set_ylim((-0.05, 1.02))
                minor_ticks = np.linspace(0.00, 1.0, 11)
                major_ticks = np.linspace(0.00, 1.0, 6)
            elif n_meas == 10:
                ax.set_ylim((-0.05, 0.90))
                minor_ticks = np.linspace(0.00, 0.90, 11)
                major_ticks = np.linspace(0.00, 0.90, 6)
            elif n_meas == 11:
                ax.set_ylim((-0.05, 0.90))
                minor_ticks = np.linspace(0.00, 0.90, 11)
                major_ticks = np.linspace(0.00, 0.90, 6)
            ax.yaxis.set_major_locator(ticker.FixedLocator(major_ticks))
            ax.yaxis.set_minor_locator(ticker.FixedLocator(minor_ticks))
            ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
            ax.yaxis.grid(True, which="minor", linestyle="--", linewidth=0.5)

    plt.tight_layout()
    plt.savefig(str(output_path / "boxplots/bp_by_ns.svg"), bbox_inches="tight")
    plt.clf()
    plt.close(fig)

    # Plots of medians:
    fig, axes = plt.subplots(3, 4, figsize=(10, 6))
    axes = axes.flatten()
    for ax, meas in zip(axes, ["TPR", "TNR", "PPV", "bACC", "GM", "F1c", "FMc", "JI", "AHR", "FGC", "F1r", "FMr"]):
        ax.plot([0, 5, 10, 15], [x for x, _, _ in Q2s_by_ns[meas]], color=colors[0], linestyle="dashdot")
        ax.plot([0, 5, 10, 15], [x for _, x, _ in Q2s_by_ns[meas]], color=colors[1], linestyle="dashdot")
        ax.plot([0, 5, 10, 15], [x for _, _, x in Q2s_by_ns[meas]], color=colors[2], linestyle="dashdot")
        ax.plot([0, 5, 10, 15], [x for x, _, _ in Q2s_by_ns[meas]], "o", color=colors[0])
        ax.plot([0, 5, 10, 15], [x for _, x, _ in Q2s_by_ns[meas]], "o", color=colors[1])
        ax.plot([0, 5, 10, 15], [x for _, _, x in Q2s_by_ns[meas]], "o", color=colors[2])
        ax.set_title(meas, fontsize=12)
        ax.grid(color="green", linestyle="--", linewidth=0.5, axis="y")
        ax.xaxis.set_major_locator(ticker.FixedLocator(np.array([0, 5, 10, 15])))
        ax.xaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
    plt.suptitle("Change in noise level", fontsize=16)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    plt.tight_layout()
    plt.savefig(str(output_path / "medians/bp_by_ns.svg"), bbox_inches="tight")
    plt.clf()
    plt.close(fig)
    print("Done.")


def global_boxplots(results: pd.DataFrame, output_path: pathlib.Path):
    """
    This function produces the boxplots of Figure 3 in the manuscript. Also, it prints in screens median values that
    are reported in the manuscript.

    Parameters
    ----------
    results: pd.DataFrame
        dataframe containing measures and other quantities for the callers under analysis (Stripebench benchmark)
    output_path: pathlib.Path
        path to output folder
    """

    # Remove stripenn:
    excluded_methods = ["stripenn"]
    results = results[~results["Method"].isin(excluded_methods)]

    for m in ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]:
        print(f"---Median values for {m} ---")
        print(f"StripePy: {np.median(results.loc[results['Method'] == 'stripepy'][m]):.4f}")
        print(f"Chromosight: {np.median(results.loc[results['Method'] == 'chromosight'][m]):.4f}")
        print(f"StripeCaller: {np.median(results.loc[results['Method'] == 'stripecaller'][m]):.4f}")
        print(f"---Interquantile range for {m} ---")
        print(f"StripePy: {np.percentile(results.loc[results['Method'] == 'stripepy'][m], 75) - np.percentile(results.loc[results['Method'] == 'stripepy'][m], 25):.4f}")
        print(f"Chromosight: {np.percentile(results.loc[results['Method'] == 'chromosight'][m], 75) - np.percentile(results.loc[results['Method'] == 'chromosight'][m], 25):.4f}")
        print(f"StripeCaller: {np.percentile(results.loc[results['Method'] == 'stripecaller'][m], 75) - np.percentile(results.loc[results['Method'] == 'stripecaller'][m], 25):.4f}")


    # GLOBAL BOXPLOTS
    fig, axes = plt.subplots(3, 4, figsize=(6.5, 4.5))
    axes = axes.flatten()
    for ax, clas_meas_name in zip(
        axes, ["TPR", "TNR", "PPV", "bACC", "GM", "F1c", "FMc", "JI", "AHR", "FGC", "F1r", "FMr"]
    ):
        sns.boxplot(
            y=clas_meas_name,
            x="Method",
            data=results,
            ax=ax,
            palette=colors[0:3],
            hue="Method",
            orient="v",
            width=0.5,
            fliersize=3,
        )
        ax.set(xlabel=None, ylabel=None)
        ax.tick_params(labelbottom=False, bottom=False)
        ax.set_yticks([0.05 * i for i in range(21)], minor=True)
        ax.set_yticks([0.2 * i for i in range(6)], minor=False)
        ax.yaxis.grid(True, which="major", linestyle="--", linewidth=0.5)
        ax.yaxis.grid(True, which="minor", linestyle="--", linewidth=0.5)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.axes.get_xaxis().set_visible(False)
        ax.set_title(clas_meas_name, fontsize=10)
        ax.set_ylim((-0.05, 1.05))
    plt.tight_layout()
    plt.savefig(str(output_path / "boxplots/bp.svg"), bbox_inches="tight")
    plt.clf()
    plt.close(fig)
    print("Done.")


def heatmaps(
    results: pd.DataFrame,
    resolutions: List[int],
    contact_densities: List[int],
    noises: List[int],
    output_path: pathlib.Path,
):
    """
    This function the heatmaps in Figure 3 in the manuscript.

    Parameters
    ----------
    results: pd.DataFrame
        dataframe containing measures and other quantities for the callers under analysis (Stripebench benchmark)
    resolutions: List[int]
        list of resolutions in the Stripebench benchmark
    contact_densities: List[int]
        list of contact densities in the Stripebench benchmark
    noises: List[int]
        list of noise levels in the Stripebench benchmark
    output_path: pathlib.Path
        path to output folder
    """

    # Retrieve results of each method:
    M1 = results.loc[results["Method"] == "stripepy"]
    M2 = results.loc[results["Method"] == "chromosight"]
    M3 = results.loc[results["Method"] == "stripecaller"]

    # RECOGNITION

    for key in ["is_anchor_found", "classification_vector"]:

        CFR = np.zeros((2, 2, 2))

        # Generate all combinations of resolutions, contact densities, and noise levels
        combinations = itertools.product(resolutions, contact_densities, noises)

        # Loop over combinations:
        for resolution, contact_density, noise in combinations:

            # Update each entry of CFR:
            if key == "is_anchor_found":

                # Retrieve all entries in M1,...,M5 with key "is_anchor_found"
                V = []
                for M in [M1, M2, M3]:
                    v = M.loc[
                        (M["Resolution"] == resolution)
                        & (M["Contact Density"] == contact_density)
                        & (M["Noise"] == noise)
                    ][f"{key}"].values[0]
                    V.append(v)
                V = np.array(V)

                # Analyse the four possible scenarios:
                for j in range(2):
                    CFR[j] += np.array(
                        [
                            [
                                np.sum((V[0] == 0) & (V[j + 1] == 0)) / len(V[0]) * 100,
                                np.sum((V[0] == 1) & (V[j + 1] == 0)) / len(V[0]) * 100,
                            ],
                            [
                                np.sum((V[0] == 0) & (V[j + 1] == 1)) / len(V[0]) * 100,
                                np.sum((V[0] == 1) & (V[j + 1] == 1)) / len(V[0]) * 100,
                            ],
                        ]
                    )

            else:

                # Retrieve all entries in M2, M3, M4 with key "classification_vector"
                V = []
                for M in [M1, M2, M3]:
                    v = M.loc[
                        (M["Resolution"] == resolution)
                        & (M["Contact Density"] == contact_density)
                        & (M["Noise"] == noise)
                    ][f"{key}"].values[0]
                    V.append(v)
                V = np.array(V)

                # Retrieve ground truth classification:
                GT = np.array(
                    M1.loc[
                        (M1["Resolution"] == resolution)
                        & (M1["Contact Density"] == contact_density)
                        & (M1["Noise"] == noise)
                    ][f"GT_{key}"].values[0]
                )
                is_anchor = np.where(GT == 1)[0]

                # Analyse the four possible scenarios:
                for j in range(2):
                    CFR[j] += np.array(
                        [
                            [
                                np.sum((V[0][is_anchor] == 0) & (V[j + 1][is_anchor] == 0)) / len(is_anchor) * 100,
                                np.sum((V[0][is_anchor] == 1) & (V[j + 1][is_anchor] == 0)) / len(is_anchor) * 100,
                            ],
                            [
                                np.sum((V[0][is_anchor] == 0) & (V[j + 1][is_anchor] == 1)) / len(is_anchor) * 100,
                                np.sum((V[0][is_anchor] == 1) & (V[j + 1][is_anchor] == 1)) / len(is_anchor) * 100,
                            ],
                        ]
                    )

        fig, axes = plt.subplots(1, 2)
        axes = axes.flatten()
        categories = ["Not Identified", "Identified"]

        for num_CFR, (CFR, ax) in enumerate(zip([CFR[0], CFR[1]], axes)):
            sns.heatmap(
                CFR.transpose() / (len(resolutions) * len(contact_densities) * len(noises)),
                annot=True,
                ax=ax,
                cmap="Blues",
                xticklabels=categories,
                yticklabels=categories,
                cbar=False,
                vmin=0,
                vmax=100,
                fmt=".1f",
                annot_kws={"size": 18},
            )
            for t in ax.texts:
                t.set_text(t.get_text() + " %")

            # x- and y-labels:
            if num_CFR == 0:
                ax.set_ylabel("StripePy", fontsize=18)
                ax.set_xlabel("Chromosight", fontsize=18)
            if num_CFR == 1:
                ax.set_xlabel("Stripecaller", fontsize=18)
            ax.tick_params(left=False, bottom=False)

            # Reduce the size of xtick labels for each subplot
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
            ax.axis("scaled")
            ax.invert_xaxis()
            ax.xaxis.set_label_position("top")
            plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
            plt.setp(ax.yaxis.get_majorticklabels(), rotation=90)
        plt.tight_layout()
        plt.savefig(str(output_path / f"heatmaps/heatmap_{key}.svg"), bbox_inches="tight")
        plt.clf()
        plt.close(fig)
    print("Done.")


def plot_RoI_and_predictions(
    Iproc_RoI,
    RoI,
    resolution,
    sites=None,
    ticks_step_size=100000,
    M1=None,
    M2=None,
    M3=None,
    M4=None,
    plot_axis=True,
    output_path=None,
):

    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.matshow(
        Iproc_RoI,
        vmax=np.amax(Iproc_RoI),
        extent=(RoI["genomic"][0], RoI["genomic"][1], RoI["genomic"][3], RoI["genomic"][2]),
        cmap=fruit_punch,
    )

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)

    ax1_btm = divider.new_vertical(size="3%", pad=0.2, pack_start=True)
    fig.add_axes(ax1_btm)
    for site in [site * resolution for site in sites[0]]:
        ax1_btm.plot([site, site], [-0.5, 0.5], color=(0.0, 0.0, 1.0), linestyle="solid", linewidth=4)

    ax1_top = divider.new_vertical(size="3%", pad=0.2, pack_start=False)
    fig.add_axes(ax1_top)
    for site in [site * resolution for site in sites[1]]:
        ax1_top.plot([site, site], [-0.5, 0.5], color=(0.0, 0.0, 1.0), linestyle="solid", linewidth=4)

    ax1_btm.set_xticks([])
    ax1_btm.set_yticks([])
    ax1_top.set_xticks([])
    ax1_top.set_yticks([])

    # stripepy:
    if M1 is not None:

        # Lower, upper or whole matrix?
        HIoIs, VIoIs, pred_clas_vec, seeds = M1
        HIoIs = HIoIs["LT"].tolist() + HIoIs["UT"].tolist()
        VIoIs = VIoIs["LT"].tolist() + VIoIs["UT"].tolist()

        # Consider candidate per candidate:
        for HIoI, VIoI in zip(HIoIs, VIoIs):

            # Correction: re-obtain genomic coordinates
            HIoI = [int(x * resolution) for x in HIoI]
            VIoI = [int(x * resolution) for x in VIoI]

            # Only candidates inside the Region of Interest will be plotted!
            if HIoI[1] < RoI["genomic"][0] or HIoI[0] > RoI["genomic"][1]:
                continue

            upp_idx_up_opt = min(RoI["genomic"][1], max(RoI["genomic"][0], VIoI[0]))
            upp_idx_dw_opt = max(RoI["genomic"][0], min(VIoI[1], RoI["genomic"][1] - resolution))
            upp_idx_sx_opt = min(max(RoI["genomic"][0], HIoI[0]), RoI["genomic"][1])
            upp_idx_dx_opt = max(RoI["genomic"][0], min(HIoI[1], RoI["genomic"][1] - resolution))

            if upp_idx_up_opt != upp_idx_dw_opt and upp_idx_sx_opt != upp_idx_dx_opt:
                ax.add_patch(
                    Rectangle(
                        (upp_idx_sx_opt, upp_idx_dw_opt),
                        upp_idx_dx_opt - upp_idx_sx_opt,
                        upp_idx_up_opt - upp_idx_dw_opt,
                        facecolor=colors[0],
                        alpha=1.0,
                        fill=True,
                        edgecolor="black",
                    )
                )

    # Stripenn -- filtered:
    if M4 is not None:

        # Lower, upper or whole matrix?
        HIoIs, VIoIs, _, _ = M4
        HIoIs = HIoIs["LT"].tolist() + HIoIs["UT"].tolist()
        VIoIs = VIoIs["LT"].tolist() + VIoIs["UT"].tolist()

        # Consider candidate per candidate:
        for HIoI, VIoI in zip(HIoIs, VIoIs):

            # Correction: re-obtain genomic coordinates
            HIoI = [int(x * resolution) for x in HIoI]
            VIoI = [int(x * resolution) for x in VIoI]

            # Only candidates inside the Region of Interest will be plotted!
            if HIoI[1] < RoI["genomic"][0] or HIoI[0] > RoI["genomic"][1]:
                continue

            upp_idx_up_opt = min(RoI["genomic"][1], max(RoI["genomic"][0], VIoI[0]))
            upp_idx_dw_opt = max(RoI["genomic"][0], min(VIoI[1], RoI["genomic"][1] - resolution))
            upp_idx_sx_opt = min(max(RoI["genomic"][0], HIoI[0]), RoI["genomic"][1])
            upp_idx_dx_opt = max(RoI["genomic"][0], min(HIoI[1], RoI["genomic"][1] - resolution))
            if upp_idx_up_opt != upp_idx_dw_opt and upp_idx_sx_opt != upp_idx_dx_opt:
                ax.add_patch(
                    Rectangle(
                        (upp_idx_sx_opt, upp_idx_dw_opt),
                        upp_idx_dx_opt - upp_idx_sx_opt,
                        upp_idx_up_opt - upp_idx_dw_opt,
                        facecolor=colors[3],
                        alpha=1.0,
                        fill=True,
                        edgecolor="black",
                    )
                )

    # Stripecaller:
    if M3 is not None:

        # Lower, upper or whole matrix?
        HIoIs, VIoIs, _, _ = M3
        HIoIs = HIoIs["LT"].tolist() + HIoIs["UT"].tolist()
        VIoIs = VIoIs["LT"].tolist() + VIoIs["UT"].tolist()

        # Consider candidate per candidate:
        for HIoI, VIoI in zip(HIoIs, VIoIs):

            # Correction: re-obtain genomic coordinates
            HIoI = [int(x * resolution) for x in HIoI]
            VIoI = [int(x * resolution) for x in VIoI]

            # Increase by one bin for plot-purposes:
            HIoI[1] = HIoI[1] + resolution

            # Only candidates inside the Region of Interest will be plotted!
            if HIoI[1] < RoI["genomic"][0] or HIoI[0] > RoI["genomic"][1]:
                continue

            upp_idx_up_opt = min(RoI["genomic"][1], max(RoI["genomic"][0], VIoI[0]))
            upp_idx_dw_opt = max(RoI["genomic"][0], min(VIoI[1], RoI["genomic"][1] - resolution))
            upp_idx_sx_opt = min(max(RoI["genomic"][0], HIoI[0]), RoI["genomic"][1])
            upp_idx_dx_opt = max(RoI["genomic"][0], min(HIoI[1], RoI["genomic"][1] - resolution))
            if (
                upp_idx_up_opt != upp_idx_dw_opt
                and upp_idx_sx_opt != upp_idx_dx_opt
                and upp_idx_up_opt < RoI["genomic"][1]
            ):
                if upp_idx_up_opt != upp_idx_dw_opt and upp_idx_sx_opt != upp_idx_dx_opt:
                    ax.add_patch(
                        Rectangle(
                            (upp_idx_sx_opt - resolution, upp_idx_dw_opt + resolution),
                            upp_idx_dx_opt - upp_idx_sx_opt,
                            upp_idx_up_opt - upp_idx_dw_opt,
                            facecolor=colors[2],
                            edgecolor="black",
                            fill=True,
                        )
                    )
    # Chromosight:
    if M2 is not None:

        # Lower, upper or whole matrix?
        HIoIs, VIoIs, _, _ = M2
        HIoIs = HIoIs["LT"].tolist() + HIoIs["UT"].tolist()
        VIoIs = VIoIs["LT"].tolist() + VIoIs["UT"].tolist()

        # Consider candidate per candidate:
        for HIoI, VIoI in zip(HIoIs, VIoIs):

            # Correction: re-obtain genomic coordinates
            HIoI = [int(x * resolution) for x in HIoI]
            VIoI = [int(x * resolution) for x in VIoI]

            # Increase by one bin for plot-purposes:
            HIoI[1] = HIoI[1] + resolution
            VIoI[1] = VIoI[1] + resolution

            # Only candidates inside the Region of Interest will be plotted!
            if (
                HIoI[1] < RoI["genomic"][0]
                or HIoI[0] > RoI["genomic"][1]
                or VIoI[1] < RoI["genomic"][0]
                or VIoI[0] > RoI["genomic"][1]
            ):
                continue

            upp_idx_up_opt = min(RoI["genomic"][1], max(RoI["genomic"][0], VIoI[0]))
            upp_idx_dw_opt = max(RoI["genomic"][0], min(VIoI[1], RoI["genomic"][1] - resolution))
            upp_idx_sx_opt = min(max(RoI["genomic"][0], HIoI[0]), RoI["genomic"][1])
            upp_idx_dx_opt = max(RoI["genomic"][0], min(HIoI[1], RoI["genomic"][1] - resolution))
            if upp_idx_up_opt != upp_idx_dw_opt and upp_idx_sx_opt != upp_idx_dx_opt:
                ax.scatter(
                    [(upp_idx_sx_opt + upp_idx_sx_opt) / 2],
                    [(upp_idx_dw_opt + upp_idx_up_opt) / 2],
                    color=colors[1],
                    marker="o",
                    edgecolors="black",
                    s=75,
                )

    fig.set_dpi(256)
    ax.set_xticks([])
    ax.set_yticks([])

    if plot_axis:
        s_adjusted = int(RoI["genomic"][0] / 100000) * 100000 + (RoI["genomic"][0] % 100000 >= 0) * 100000
        e_adjusted = int(RoI["genomic"][1] / 100000) * 100000 + (RoI["genomic"][1] % 100000 >= 0) * 100000
        ax1_btm.set_xticks(range(s_adjusted, e_adjusted, 2 * ticks_step_size))
        format_ticks(ax1_btm)
        # ax1_btm.tick_params(labelsize=20, length=15, width=4)
        ax1_btm.tick_params(labelsize=25, size=13, length=15)

    ax1_top.set_xlim([RoI["genomic"][0], RoI["genomic"][1] + resolution])
    ax1_btm.set_xlim([RoI["genomic"][0], RoI["genomic"][1] + resolution])
    fig.tight_layout()
    if output_path is not None:
        plt.savefig(output_path, bbox_inches="tight")
    plt.clf()
    plt.close(fig)


def real_data_csv_measures(
    metrics: Dict[str, Dict[str, float]],
    n_found_anchors: Dict[str, int],
    n_predicted_stripes: Dict[str, int],
    wall_clock_times: List[int],
    output_path: pathlib.Path,
):
    """
    This function produces Supplementary Table 6 in the manuscript (one contact map per run).

    Parameters
    ----------
    metrics: Dict[str, Dict[str, float]]
        dictionary of dictionaries, it stores classification and  recognition measures for each stripe caller.
    n_found_anchors: Dict[str, int]
        dictionary that contains, for each stripe caller, the number of anchor sites found inside a predicted stripe
    n_predicted_stripes: Dict[str, int]
        dictionary that contains, for each stripe caller, the number of predicted stripes containing anchor sites
    wall_clock_times: List[int]
        list of the elapsed real times
    output_path: pathlib.Path
        path to output folder
    """

    # List of lists, that will contain each row of the dataframe:
    rows = []

    # Line starts:
    for n_meas, meas in enumerate(["TPR", "TNR", "PPV", "bACC", "GM", "F1c", "FMc", "JI", "AHR", "FGC", "F1r", "FMr"]):

        # New line:
        row = [meas]

        for M in ["M1", "M2", "M3", "M4"]:
            row.append(f"{metrics[M][meas] * 100:>6.2f}")

        rows.append(row)

    # Number of anchors found (nAF) and number of stripes predicted (nSP):
    row = ["nAF"]
    for M in ["M1", "M2", "M3", "M4"]:
        row.append(f"{n_found_anchors[M]}")
    rows.append(row)
    row = ["nSP"]
    for M in ["M1", "M2", "M3", "M4"]:
        row.append(f"{n_predicted_stripes[M]}")
    rows.append(row)

    # Wall-clock times:
    rows.append(["time"] + wall_clock_times)

    # Store results in a pandas dataframe:
    df = pd.DataFrame(rows, columns=["measure", "M1", "M2", "M3", "M4"])

    df.to_csv(output_path, index=False)
    print("Done.")


def csv_normalization_tables(metrics, n_found_anchors, n_predicted_stripes, output_path):

    # Contact maps:
    cmap_names = list(metrics.keys())

    # List of lists, that will contain each row of the dataframe:
    rows = []

    # Line starts:
    for n_meas, meas in enumerate(["TPR", "TNR", "PPV", "bACC", "GM", "F1c", "FMc", "JI", "AHR", "FGC", "F1r", "FMr"]):

        # New line:
        row = [meas]

        for cmap_name in cmap_names:
            for M in ["NONE", "GW_ICE", "GW_SCALE"]:
                row.append(f"{metrics[cmap_name][M][meas] * 100:>6.2f}")

        rows.append(row)

    # Number of anchors found (nAF) and number of stripes predicted (nSP):
    row = ["nAF"]
    for cmap_name in cmap_names:
        for M in ["NONE", "GW_ICE", "GW_SCALE"]:
            row.append(f"{n_found_anchors[cmap_name][M]}")
    rows.append(row)
    row = ["nSP"]
    for cmap_name in cmap_names:
        for M in ["NONE", "GW_ICE", "GW_SCALE"]:
            row.append(f"{n_predicted_stripes[cmap_name][M]}")
    rows.append(row)

    # Store results in a pandas dataframe:
    df = pd.DataFrame(
        rows,
        columns=[
            "measure",
            "NONE",
            "ICE",
            "SCALE",
            "NONE",
            "ICE",
            "SCALE",
            "NONE",
            "ICE",
            "SCALE",
            "NONE",
            "ICE",
            "SCALE",
        ],
    )

    df.to_csv(output_path, index=False)
    print("Done.")
