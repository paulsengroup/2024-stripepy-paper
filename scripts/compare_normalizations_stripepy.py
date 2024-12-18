#!/usr/bin/env python3

# Copyright (C) 2024 Andrea Raffo <andrea.raffo@ibv.uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix

from utils import IO, evaluate

# Resolution
resolution = 10000

# Output path:
output_path = pathlib.Path("./output/real data/tables/")


def my_parse_args():

    class CustomFormatter(argparse.RawTextHelpFormatter):
        def _fill_text(self, text, width, indent):
            return "".join([indent + line + "\n" for line in text.splitlines()])

    def _existing_path(arg):
        path = pathlib.Path(arg).parent
        if path.exists() and path.is_dir():
            return arg

        raise FileNotFoundError(f"Path not reachable: {path}")

    cli = argparse.ArgumentParser(
        description="Routine to evaluate StripePy when different normalizations are adopted. This "
        "script requires to download the following files and folders:\n"
        "• the real contact maps\n"
        "• the ChIP-seq peaks\n"
        "• StripePy's output folder\n"
        "The location of the contact maps and ChIP-seq peaks are documented in the StripePy manuscript.\n"
        "The output folders of the stripe callers are here collected: https://zenodo.org/records/14449731.\n",
        formatter_class=CustomFormatter,
    )

    cli.add_argument(
        "--path-to-real-data",
        type=_existing_path,
        required=True,
        help="Path to the folder containing the used contact maps.",
    )

    cli.add_argument(
        "--path-to-chip-seq-peaks",
        type=_existing_path,
        required=True,
        help="Path to the folder containing the ChIP-seq peaks in .bed.gz format.",
    )

    cli.add_argument(
        "--path-to-stripepy",
        type=_existing_path,
        required=True,
        help="Path to the folder containing the runs produced by stripepy (e.g., /tmp/stripepy/real data/).",
    )

    # Parse the input parameters:
    args = vars(cli.parse_args())

    # Gather input parameters in dictionaries:
    configs_input = {key: args[key] for key in ["path_to_real_data", "path_to_chip_seq_peaks", "path_to_stripepy"]}

    # Print the used parameters (chosen or default-ones):
    print("\nArguments:")
    print(f"contact-map: {configs_input['path_to_real_data']}")
    print(f"chip-seq: {configs_input['path_to_chip_seq_peaks']}")
    print(f"stripepy: {configs_input['path_to_stripepy']}")

    return configs_input


if __name__ == "__main__":

    configs_input = my_parse_args()

    # Define the configuration details in a JSON format
    configs = {
        "4DNFI9GMP2J8": {
            "contact-map": f"{configs_input['path_to_real_data']}/4DNFI9GMP2J8.fixed.mcool",
            "chip-seq-peaks": f"{configs_input['path_to_chip_seq_peaks']}/ENCFF692RPA.bed.gz",
            "path-to-stripepy": f"{configs_input['path_to_stripepy']}/4DNFI9GMP2J8.fixed",
        },
        "4DNFI6HDY7WZ": {
            "contact-map": f"{configs_input['path_to_real_data']}/4DNFI6HDY7WZ.fixed.mcool",
            "chip-seq-peaks": f"{configs_input['path_to_chip_seq_peaks']}/ENCFF692RPA.bed.gz",
            "path-to-stripepy": f"{configs_input['path_to_stripepy']}/4DNFI6HDY7WZ.fixed",
        },
        "ENCFF993FGR": {
            "contact-map": f"{configs_input['path_to_real_data']}/ENCFF993FGR.mcool",
            "chip-seq-peaks": f"{configs_input['path_to_chip_seq_peaks']}/ENCFF796WRU.bed.gz",
            "path-to-stripepy": f"{configs_input['path_to_stripepy']}/ENCFF993FGR",
        },
        "ENCFF216QQM": {
            "contact-map": f"{configs_input['path_to_real_data']}/ENCFF216QQM.fixed.mcool",
            "chip-seq-peaks": f"{configs_input['path_to_chip_seq_peaks']}/ENCFF796WRU.bed.gz",
            "path-to-stripepy": f"{configs_input['path_to_stripepy']}/ENCFF216QQM.fixed",
        },
    }

    # Prefixes:
    normalizations = ["NONE", "GW_ICE", "GW_SCALE"]

    # Metrics:
    metrics = {cm_name: {} for cm_name in configs.keys()}

    # Numbers of anchors/stripes:
    n_found_anchors = {cm_name: {} for cm_name in configs.keys()}  # number of found anchor sites
    n_predicted_stripes = {cm_name: {} for cm_name in configs.keys()}  # number of predicted stripes

    for cm_name in configs.keys():

        print(f"---Contact map: {cm_name}---")

        # Contact map loading:
        c, _, _, _ = IO.cmap_loading(str(configs[cm_name]["contact-map"]), resolution)
        chrs = list(c.chromosomes().keys())
        chr_ids = list(range(len(chrs)))
        c_pairs = list(zip(chr_ids, chrs))

        # Ingredients for recognition measures:
        is_anchor_found = {n: [] for n in normalizations}  # 0-1 classification vector
        is_candidate_good = {n: [] for n in normalizations}  # 0-1 classification vector

        # Numbers of anchors/stripes:
        n_found_anchors[cm_name] = {n: 0 for n in normalizations}  # number of found anchor sites
        n_predicted_stripes[cm_name] = {n: 0 for n in normalizations}  # number of predicted stripes

        # Classification and recognition measures:
        confusion_matrices = {n: {key: 0 for key in ["TP", "TN", "FP", "FN"]} for n in normalizations}
        metrics_keys = ["TPR", "TNR", "PPV", "bACC", "GM", "JI", "F1c", "FMc", "AHR", "FGC", "F1r", "FMr"]
        metrics[cm_name] = {n: {key: 0 for key in metrics_keys} for n in normalizations}

        # ChIP-seq peaks:
        header_names = ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
        peaks = pd.read_csv(
            configs[cm_name]["chip-seq-peaks"], compression="gzip", sep="\t", header=None, names=header_names
        )
        peaks["anchor"] = (peaks["start"] + peaks["end"]) / 2
        peaks = peaks.astype({"anchor": int})

        # Loop over chromosomes:
        for chr_idx, chr in c_pairs[:-3]:

            print(f"Chr: {chr}")

            # Contact map for this chromosome:
            I_NONE = c.fetch(chr, normalization="NONE").to_coo().tocsr()

            # Peaks for this chromosome:
            these_peaks = peaks[peaks["chrom"] == chr]

            # Ground-truth classification vector:
            GT_anchors = np.unique(np.round(these_peaks.anchor.values / resolution)).astype(int).tolist()
            GT_anchors = np.unique(GT_anchors)
            GT_clas_vec = np.isin(np.arange(I_NONE.shape[0]), GT_anchors).astype(int)

            # Method M1 - NONE:
            path = pathlib.Path(f"{configs[cm_name]['path-to-stripepy']}/{resolution}/")
            HIoIs, _, pred_clas_vec, seeds = IO.retrieve_stripepy(path, chr, I_NONE.shape[0], 3.0)

            # Classification measures:
            pred_clas_vec = pred_clas_vec["LT"] + pred_clas_vec["UT"]
            pred_clas_vec = np.where(pred_clas_vec > 1, 1, pred_clas_vec)
            tn, fp, fn, tp = confusion_matrix(GT_clas_vec, pred_clas_vec, labels=[0, 1]).ravel()
            confusion_matrices["NONE"]["TP"] += tp
            confusion_matrices["NONE"]["TN"] += tn
            confusion_matrices["NONE"]["FP"] += fp
            confusion_matrices["NONE"]["FN"] += fn

            # Recognition measures:
            x = evaluate._is_anchor_in_stripes(
                GT_anchors=GT_anchors, pred_HIoIs=HIoIs["LT"].tolist() + HIoIs["UT"].tolist()
            )
            is_anchor_found["NONE"] += x[0].tolist()
            is_candidate_good["NONE"] += x[1].tolist()

            # Numbers:
            n_found_anchors[cm_name]["NONE"] += np.sum(x[0])
            n_predicted_stripes[cm_name]["NONE"] += HIoIs["LT"].shape[0] + HIoIs["UT"].shape[0]

            # Method M1 - GW_ICE:
            path = pathlib.Path(f"{configs[cm_name]['path-to-stripepy']}/{resolution}-GW_ICE/")
            HIoIs, _, pred_clas_vec, seeds = IO.retrieve_stripepy(path, chr, I_NONE.shape[0], 3.0)

            # Classification measures:
            pred_clas_vec = pred_clas_vec["LT"] + pred_clas_vec["UT"]
            pred_clas_vec = np.where(pred_clas_vec > 1, 1, pred_clas_vec)
            tn, fp, fn, tp = confusion_matrix(GT_clas_vec, pred_clas_vec, labels=[0, 1]).ravel()
            confusion_matrices["GW_ICE"]["TP"] += tp
            confusion_matrices["GW_ICE"]["TN"] += tn
            confusion_matrices["GW_ICE"]["FP"] += fp
            confusion_matrices["GW_ICE"]["FN"] += fn

            # Recognition measures:
            x = evaluate._is_anchor_in_stripes(
                GT_anchors=GT_anchors, pred_HIoIs=HIoIs["LT"].tolist() + HIoIs["UT"].tolist()
            )
            is_anchor_found["GW_ICE"] += x[0].tolist()
            is_candidate_good["GW_ICE"] += x[1].tolist()

            # Numbers:
            n_found_anchors[cm_name]["GW_ICE"] += np.sum(x[0])
            n_predicted_stripes[cm_name]["GW_ICE"] += HIoIs["LT"].shape[0] + HIoIs["UT"].shape[0]

            # Method M1 - GW_SCALE:
            path = pathlib.Path(f"{configs[cm_name]['path-to-stripepy']}/{resolution}-GW_SCALE/")
            HIoIs, _, pred_clas_vec, seeds = IO.retrieve_stripepy(path, chr, I_NONE.shape[0], 3.0)

            # Classification measures:
            pred_clas_vec = pred_clas_vec["LT"] + pred_clas_vec["UT"]
            pred_clas_vec = np.where(pred_clas_vec > 1, 1, pred_clas_vec)
            tn, fp, fn, tp = confusion_matrix(GT_clas_vec, pred_clas_vec, labels=[0, 1]).ravel()
            confusion_matrices["GW_SCALE"]["TP"] += tp
            confusion_matrices["GW_SCALE"]["TN"] += tn
            confusion_matrices["GW_SCALE"]["FP"] += fp
            confusion_matrices["GW_SCALE"]["FN"] += fn

            # Recognition measures:
            x = evaluate._is_anchor_in_stripes(
                GT_anchors=GT_anchors, pred_HIoIs=HIoIs["LT"].tolist() + HIoIs["UT"].tolist()
            )
            is_anchor_found["GW_SCALE"] += x[0].tolist()
            is_candidate_good["GW_SCALE"] += x[1].tolist()

            # Numbers:
            n_found_anchors[cm_name]["GW_SCALE"] += np.sum(x[0])
            n_predicted_stripes[cm_name]["GW_SCALE"] += HIoIs["LT"].shape[0] + HIoIs["UT"].shape[0]

        # Compute the averages of chromosomes:
        evaluate.compute_classification_measures(confusion_matrices, metrics[cm_name])
        evaluate.compute_recognition_measures(is_anchor_found, is_candidate_good, metrics[cm_name])

    IO.csv_normalization_tables(metrics, n_found_anchors, n_predicted_stripes, output_path / "table7.csv")
