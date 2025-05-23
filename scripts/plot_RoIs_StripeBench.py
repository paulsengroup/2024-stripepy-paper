#!/usr/bin/env python3

# Copyright (C) 2024 Andrea Raffo <andrea.raffo@ibv.uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import os
import pathlib

import numpy as np
import pandas as pd

from utils import IO, evaluate

# Output folder
output_path = pathlib.Path("./output/StripeBench/RoIs/")


def parse_args():

    class CustomFormatter(argparse.RawTextHelpFormatter):
        def _fill_text(self, text, width, indent):
            return "".join([indent + line + "\n" for line in text.splitlines()])

    def _existing_path(arg):
        path = pathlib.Path(arg)
        if path.exists() and path.is_dir():
            return path

        raise FileNotFoundError(f"Path not reachable: {path}")

    cli = argparse.ArgumentParser(
        description="Routine to plot RoIs from the StripeBench benchmark. This "
        "script requires to download the following folders:\n"
        "• StripeBench\n"
        "• StripePy's output folder\n"
        "• Chromosight's output folder\n"
        "• StripeCaller's output folder\n"
        "StripeBench is available here: https://zenodo.org/records/14448329.\n"
        "The output folders of the stripe callers are here collected: https://zenodo.org/records/14449731.\n",
        formatter_class=CustomFormatter,
    )

    cli.add_argument(
        "--stripebench",
        type=_existing_path,
        required=True,
        help="Path to the folder containing the StripeBench benchmark (e.g., '/tmp/StripeBench').",
    )

    cli.add_argument(
        "--stripepy",
        type=_existing_path,
        required=True,
        help="Path to the folder containing the output files generated by StripePy when executed on StripeBench "
        "(e.g., /tmp/stripepy/StripeBench).",
    )

    cli.add_argument(
        "--chromosight",
        type=_existing_path,
        required=True,
        help="Path to the folder containing the output files generated by Chromosight when executed on StripeBench "
        "(e.g., /tmp/chromosight/StripeBench).",
    )

    cli.add_argument(
        "--stripecaller",
        type=_existing_path,
        required=True,
        help="Path to the folder containing the output files generated by StripeCaller when executed on StripeBench "
        "(e.g., /tmp/stripecaller/StripeBench).",
    )

    # Parse the input parameters:
    args = vars(cli.parse_args())

    # Gather input parameters in dictionaries:
    configs_paths = {key: args[key] for key in ["stripebench", "stripepy", "chromosight", "stripecaller"]}

    # Print the used parameters (chosen or default-ones):
    print("\nArguments:")
    print(f"Path to StripeBench: {configs_paths['stripebench']}")
    print(f"Path to StripePy: {configs_paths['stripepy']}")
    print(f"Path to Chromosight: {configs_paths['chromosight']}")
    print(f"Path to StripeCaller: {configs_paths['stripecaller']}")

    return configs_paths


if __name__ == "__main__":

    probability_score = 0.7

    # Import parameters:
    configs_paths = parse_args()

    for which_roi in ["roi_1", "roi_2", "roi_3"]:

        # RoI 1 (contact density = 1, noise = 5000, resolution = 5000)
        if which_roi == "roi_1":
            chromosome = "chr2"
            contact_density = 1
            noise_level = 5000
            resolution = 5000

        # RoI 2 (contact density = 5, noise = 10000, resolution = 10000)
        elif which_roi == "roi_2":
            chromosome = "chr5"
            contact_density = 5
            noise_level = 10000
            resolution = 10000

        # RoI 3 (contact density = 15, noise = 15000, resolution = 25000)
        elif which_roi == "roi_3":
            chromosome = "chr9"
            contact_density = 15
            noise_level = 15000
            resolution = 25000
        elif which_roi == "roi_4":
            chromosome = "chr3"
            contact_density = 1
            noise_level = 5000
            resolution = 5000

        # Data loading:
        path2mcool = (
            configs_paths["stripebench"]
            / "data"
            / f"grch38_h1_rad21_{contact_density}_{noise_level}"
            / f"grch38_h1_rad21_{contact_density}_{noise_level}.mcool"
        )
        c, chr_starts, chr_ends, bp_lengths = IO.cmap_loading(str(path2mcool), resolution)

        if which_roi == "roi_1":
            RoI_length = int(1000000 / resolution)
            RoI_m = int((chr_ends[1] - chr_starts[1]) / 2)
            RoI_s = RoI_m + int(5 * RoI_length)
            RoI_e = RoI_m + int(6 * RoI_length)
            ticks_step_size = 100000
        elif which_roi == "roi_2":
            RoI_length = int(1000000 / resolution)
            RoI_m = int((chr_ends[1] - chr_starts[1]) / 2)
            RoI_s = RoI_m + int(2 * RoI_length)
            RoI_e = RoI_m + int(2 * RoI_length) + int(RoI_length)
            ticks_step_size = 100000
        elif which_roi == "roi_3":
            RoI_length = int(1000000 / resolution)
            RoI_m = int((chr_ends[1] - chr_starts[1]) / 2)
            RoI_s = RoI_m
            RoI_e = RoI_m + int(RoI_length)
            ticks_step_size = 100000
        elif which_roi == "roi_4":
            RoI_length = int(1000000 / resolution)
            RoI_m = int((chr_ends[1] - chr_starts[1]) / 2)
            RoI_s = RoI_m + int(5 / 2 * RoI_length)
            RoI_e = RoI_m + int(RoI_length * 4 / 6) + int(5 / 2 * RoI_length)
            ticks_step_size = 100000
        else:
            RoI_length = int(1000000 / resolution)
            RoI_m = int((chr_ends[1] - chr_starts[1]) / 2)
            RoI_s = RoI_m
            RoI_e = RoI_m + int(RoI_length * 4 / 6)
            ticks_step_size = 100000

        # Retrieve contact map:
        I = c.fetch(chromosome).to_coo().tolil()
        number_of_bins = I.shape[0]

        # Symmetrization:
        I += I.T
        I.setdiag(I.diagonal() / 2)
        I = I.tocsr()

        # Log-transform and standardization:
        Iproc = I.log1p()
        Iproc /= Iproc.max()

        # ChIP-seq peaks:
        path2chip = configs_paths["stripebench"] / "GRCh38_H1_RAD21_occupancy.bed"
        header_names = ["chrom", "start", "end", "?", "score", "sign"]
        df = pd.read_csv(path2chip, sep="\t", header=None, names=header_names)
        df["anchor"] = (df["start"] + df["end"]) / 2
        df = df.astype({"anchor": int})
        df = df[df["chrom"] == chromosome]
        df = df[df["score"] >= probability_score]
        LT_df = df[df["sign"] == "+"]
        UT_df = df[df["sign"] == "-"]

        # Ground-truth classification vector:
        LT_GT_anchors = np.unique(np.round(LT_df.anchor.values / resolution)).astype(int).tolist()
        UT_GT_anchors = np.unique(np.round(UT_df.anchor.values / resolution)).astype(int).tolist()

        # RoI extraction:
        RoI = dict()
        RoI["matrix"] = [RoI_s, RoI_e, RoI_s, RoI_e]
        RoI["genomic"] = [int(roi * resolution) for roi in RoI["matrix"]]
        rows = cols = slice(RoI["matrix"][0], RoI["matrix"][1])
        Iproc_RoI = Iproc[rows, cols].toarray()

        # Ground truth anchors in RoI:
        LT_GT_anchors_in_RoI = [x for x in LT_GT_anchors if RoI["matrix"][0] <= x <= RoI["matrix"][1]]
        UT_GT_anchors_in_RoI = [x for x in UT_GT_anchors if RoI["matrix"][0] <= x <= RoI["matrix"][1]]
        LT_GT_anchors_in_RoI = evaluate.average_within_threshold(LT_GT_anchors_in_RoI, 1)
        UT_GT_anchors_in_RoI = evaluate.average_within_threshold(UT_GT_anchors_in_RoI, 1)
        # IO.HiC(Iproc_RoI, RoI['genomic'], [LT_GT_anchors_in_RoI, UT_GT_anchors_in_RoI], resolution,
        #        plot_in_bp=True, show=True)

        # Get the file name without extension
        file_name = os.path.splitext(os.path.basename(path2mcool))[0]

        # Retrieve output of stripe callers:
        path2M1 = configs_paths["stripepy"] / f"{file_name}" / f"{resolution}" / ""
        path2M2 = configs_paths["chromosight"] / f"{file_name}" / f"{resolution}" / ""
        path2M3 = configs_paths["stripecaller"] / f"{file_name}" / f"{resolution}" / ""
        M1 = IO.retrieve_stripepy(path2M1, chromosome, number_of_bins, 4.0)
        M2 = IO.retrieve_chromosight(path2M2, chromosome, I.shape[0], resolution)
        M3 = IO.retrieve_stripecaller(path2M3, chromosome, I.shape[0], resolution)

        # IO.plot_RoI_and_predictions(Iproc_RoI, RoI, resolution, [LT_GT_anchors_in_RoI, UT_GT_anchors_in_RoI])
        IO.plot_RoI_and_predictions(
            Iproc_RoI,
            RoI,
            resolution,
            [LT_GT_anchors_in_RoI, UT_GT_anchors_in_RoI],
            ticks_step_size,
            M1=M1,
            M2=M2,
            M3=M3,
            output_path=pathlib.Path(f"{output_path}/{which_roi}.svg"),
        )
    print("Done.")
