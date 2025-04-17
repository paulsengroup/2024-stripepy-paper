#!/usr/bin/env python3

# Copyright (C) 2024 Andrea Raffo <andrea.raffo@ibv.uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import itertools
import pathlib

import hictkpy
import pandas as pd
from utils import IO, evaluate

## PARAMETERS

# Parameters for the StripeBench benchmark:
file_name_base = "grch38_h1_rad21"  # matrices are named "f{file_name_base}_{contact_density}_{noise}.mcool"
GT_name = "GRCh38_H1_RAD21_occupancy.bed"  # ground truth annotations
prob_threshold = 0.70  # to extract anchor points from the ground truth
noises = [0, 5000, 10000, 15000]  # levels of noise
resolutions = [5000, 10000, 25000, 50000]  # resolutions of Hi-C matrices
contact_densities = [1, 5, 10, 15]  # contact densities (i.e., sequencing depth)

# Relative change threshold for StripePy:
rel_change_threshold = 5.0

# Output folder
output_path = pathlib.Path("./output/StripeBench/")


def my_parse_args():

    def _dir_checked(arg: str) -> pathlib.Path:
        input_dir = pathlib.Path(arg)
        if input_dir.exists() and input_dir.is_dir():
            return pathlib.Path(arg)

        raise FileNotFoundError(f'Input folder "{arg}" is not reachable: folder does not exist')

    class CustomFormatter(argparse.RawTextHelpFormatter):
        def _fill_text(self, text, width, indent):
            return "".join([indent + line + "\n" for line in text.splitlines()])

    cli = argparse.ArgumentParser(
        description="Routine to evaluate StripePy, Chromosight, StripeCaller and Stripenn over StripeBench. This "
        "script requires to download the following folders:\n"
        "• StripeBench\n"
        "• StripePy's output folder\n"
        "• Chromosight's output folder\n"
        "• StripeCaller's output folder\n"
        "• Stripenn's output folder\n"
        "StripeBench is available here: https://zenodo.org/records/14448329.\n"
        "The output folders of the stripe callers are here collected: https://zenodo.org/records/14449731.\n",
        formatter_class=CustomFormatter,
    )

    cli.add_argument(
        "--stripebench",
        type=_dir_checked,
        required=True,
        help="Path to the folder containing the StripeBench benchmark (e.g., '/tmp/StripeBench').",
    )

    cli.add_argument(
        "--stripepy",
        type=_dir_checked,
        required=True,
        help="Path to the folder containing the output files generated by StripePy when executed on StripeBench (e.g., '/tmp/stripepy').",
    )

    cli.add_argument(
        "--chromosight",
        type=_dir_checked,
        required=True,
        help="Path to the folder containing the output files generated by Chromosight when executed on StripeBench (e.g., '/tmp/chromosight').",
    )

    cli.add_argument(
        "--stripecaller",
        type=_dir_checked,
        required=True,
        help="Path to the folder containing the output files generated by StripeCaller when executed on StripeBench (e.g., '/tmp/stripecaller').",
    )

    cli.add_argument(
        "--stripenn",
        type=_dir_checked,
        required=True,
        help="Path to the folder containing the output files generated by Stripenn when executed on StripeBench (e.g., '/tmp/stripenn').",
    )

    return cli.parse_args()


if __name__ == "__main__":

    # Parse the input parameters:
    args = vars(my_parse_args())

    # Gather input parameters in dictionaries:
    configs_paths = {key: args[key] for key in ["stripebench", "stripepy", "chromosight", "stripecaller", "stripenn"]}

    # Print the used parameters (chosen or default-ones):
    print("\nArguments:")
    print(f"Path to StripeBench: {configs_paths['stripebench']}")
    print(f"Path to StripePy: {configs_paths['stripepy']}")
    print(f"Path to Chromosight: {configs_paths['chromosight']}")
    print(f"Path to StripeCaller: {configs_paths['stripecaller']}")
    print(f"Path to Stripenn: {configs_paths['stripenn']}")

    # Dataframe initialization:
    result_rows = []
    column_names = [
        "Method",
        "Resolution",
        "Noise",
        "Contact Density",
        "is_anchor_found",
        "classification_vector",
        "GT_classification_vector",
        "TPR",
        "TNR",
        "PPV",
        "bACC",
        "GM",
        "JI",
        "AHR",
        "FGC",
        "F1c",
        "FMc",
        "F1r",
        "FMr",
    ]

    print(f"\nCOMPUTE ALL MEASURES FOR ALL COMBINATION OF SINGLE FACTORS...")

    # Generate all combinations of resolutions, contact densities, and noise levels
    combinations = itertools.product(resolutions, contact_densities, noises)

    for resolution, contact_density, noise in combinations:

        # Current Hi-C matrix:
        file_name = f"{file_name_base}_{contact_density}_{noise}"

        # Paths to relevant files:
        path2mcool = configs_paths["stripebench"] / "data" / f"{file_name}"
        path2GT = configs_paths["stripebench"]
        path2stripepy = configs_paths["stripepy"] / "StripeBench" / f"{file_name}" / f"{resolution}"
        path2chromosight = configs_paths["chromosight"] / "StripeBench" / f"{file_name}" / f"{resolution}"
        path2stripecaller = configs_paths["stripecaller"] / "StripeBench" / f"{file_name}" / f"{resolution}"
        path2stripenn = configs_paths["stripenn"] / "StripeBench" / f"{file_name}" / f"{resolution}"

        # Update on resolution:
        print(f"Resolution: {resolution}, contact density: {contact_density}, noise level: {noise}...")

        # 1) DATA LOADING
        c = hictkpy.File(str(path2mcool / f"{file_name}.mcool::resolutions" / f"{resolution}"), resolution)
        c_names = list(c.chromosomes().keys())
        c_ids = list(range(len(c_names)))
        c_pairs = list(zip(c_ids, c_names))

        # We use the following encoding:
        # M1 = StripePy, M2 = Chromosight, M3 = StripeCaller, M4 = Stripenn
        M1 = evaluate.initialize_evaluation_dictionary()
        M2 = evaluate.initialize_evaluation_dictionary()
        M3 = evaluate.initialize_evaluation_dictionary()
        M4 = evaluate.initialize_evaluation_dictionary()

        # Ground-truth classification vector:
        GT_clas_vec = []

        # 2) RECOGNITION MEASURES
        tot_n_gt_anchors = 0  # Number of ground truth anchor sites computed over all chromosomes

        # Looping over chromosomes:

        for cur_chr in c_names:
            # Number of bins in the current chromosome:
            bins = c.bins().to_df()
            n_bins = len(bins[bins["chrom"] == cur_chr])

            # Ground truth:
            GT_chr = IO.retrieve_ground_truth(path2GT, GT_name, cur_chr, resolution, prob_threshold, n_bins)
            tot_n_gt_anchors += len(GT_chr["anchors"]["LT"].tolist()) + len(GT_chr["anchors"]["UT"].tolist())

            # Retrieve data for the methods:
            M1_pred = IO.retrieve_stripepy(path2stripepy, cur_chr, n_bins, rel_change_threshold)
            M2_pred = IO.retrieve_chromosight(path2chromosight, cur_chr, n_bins, resolution)
            M3_pred = IO.retrieve_stripecaller(path2stripecaller, cur_chr, n_bins, resolution)
            M4_pred = IO.retrieve_stripenn(path2stripenn, cur_chr, n_bins, resolution, filter_=True)

            # Compare prediction to ground truth values:
            # Order of the input is:  HIoIs, VIoIs, clas_vec, seeds
            M1_chr = evaluate._compare_predictions_to_StripeBench_GT(GT_chr, HIoIs=M1_pred[0], clas_vec=M1_pred[2])
            M2_chr = evaluate._compare_predictions_to_StripeBench_GT(GT_chr, HIoIs=M2_pred[0], clas_vec=M2_pred[2])
            M3_chr = evaluate._compare_predictions_to_StripeBench_GT(GT_chr, HIoIs=M3_pred[0], clas_vec=M3_pred[2])
            M4_chr = evaluate._compare_predictions_to_StripeBench_GT(GT_chr, HIoIs=M4_pred[0], clas_vec=M4_pred[2])

            # For global evaluation over all chromosomes:
            for key in ["is_candidate_good", "is_anchor_found", "TN", "FP", "FN", "TP"]:
                M1[key] += M1_chr[key]
                M2[key] += M2_chr[key]
                M3[key] += M3_chr[key]
                M4[key] += M4_chr[key]

            # Predicted classification vectors:
            M1["classification_vector"] += M1_pred[2]["LT"].tolist() + M1_pred[2]["UT"].tolist()
            M2["classification_vector"] += M2_pred[2]["LT"].tolist() + M2_pred[2]["UT"].tolist()
            M3["classification_vector"] += M3_pred[2]["LT"].tolist() + M3_pred[2]["UT"].tolist()
            M4["classification_vector"] += M4_pred[2]["LT"].tolist() + M4_pred[2]["UT"].tolist()

            # GT classification vectors (useful for heatmaps):
            M1["GT_classification_vector"] += GT_chr["clas_vec"]["LT"].tolist() + GT_chr["clas_vec"]["UT"].tolist()
            M2["GT_classification_vector"] += GT_chr["clas_vec"]["LT"].tolist() + GT_chr["clas_vec"]["UT"].tolist()
            M3["GT_classification_vector"] += GT_chr["clas_vec"]["LT"].tolist() + GT_chr["clas_vec"]["UT"].tolist()
            M4["GT_classification_vector"] += GT_chr["clas_vec"]["LT"].tolist() + GT_chr["clas_vec"]["UT"].tolist()

        print(f"->) This map has {tot_n_gt_anchors} ground truth anchors")

        # Loop over all methods:
        for M, M_name in zip([M1, M2, M3, M4], ["stripepy", "chromosight", "stripecaller", "stripenn"]):
            TPR, TNR, PPV, bACC, GM, JI, AHR, FGC, F1c, FMc, F1r, FMr = evaluate._compute_measures_StripeBench(M)

            result_rows.append(
                [
                    M_name,
                    resolution,
                    noise,
                    contact_density,
                    M["is_anchor_found"],
                    M["classification_vector"],
                    M["GT_classification_vector"],
                    TPR,
                    TNR,
                    PPV,
                    bACC,
                    GM,
                    JI,
                    AHR,
                    FGC,
                    F1c,
                    FMc,
                    F1r,
                    FMr,
                ]
            )

    # Store results in a pandas dataframe:
    results = pd.DataFrame(result_rows, columns=column_names)

    # TABLES AND PLOTS
    print(f"\nCOUNT TABLE...")
    IO.StripeBench_csv_counts(results, resolutions, contact_densities, noises, output_path)

    print(f"\nTABLES CONTAINING CLASSIFICATION AND RECOGNITION MEASURES...")
    IO.StripeBench_csv_measures(results, resolutions, contact_densities, noises, output_path)

    # MARGINAL BOX PLOTS AND PLOTS OF THE MEDIANS
    print(f"\nMARGINAL PLOTS...")
    IO.marginal_plots(results, resolutions, contact_densities, noises, output_path)

    # GLOBAL BOXPLOTS
    print(f"\nGLOBAL BOXPLOTS...")
    IO.global_boxplots(results, output_path)

    # HEATMAPS
    print(f"\nHEATMAPS...")
    IO.heatmaps(results, resolutions, contact_densities, noises, output_path)

    print("\n\n-------------------")
    print("----SUPPLEMENTARY TABLES AT DIFFERENT RESOLUTIONS----")
    IO.StripeBench_LaTex_tables(results, resolutions, contact_densities, noises)
    print("-------------------")

    print("\n\n-------------------")
    print("----COUNT TABLE AT DIFFERENT RESOLUTIONS----")
    IO.StripeBench_LaTex_table_counts(results, resolutions, contact_densities, noises)
    print("-------------------")
