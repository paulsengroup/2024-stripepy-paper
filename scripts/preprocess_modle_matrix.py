#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import logging
import pathlib

import h5py
import numpy as np


def positive_int(arg) -> float:
    if (n := int(arg)) > 0:
        return n

    raise ValueError("Not a positive int")


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser(
        'Extract the single-resolution matrices from .mcool files and add "fake" balancing weights.'
    )

    cli.add_argument(
        "mcools",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more .mcool files.",
    )
    cli.add_argument(
        "--resolutions",
        required=True,
        type=positive_int,
        nargs="+",
        help="One or more resolution(s) to be extracted.",
    )
    cli.add_argument(
        "--output-dir",
        type=pathlib.Path,
        help="Path to a folder where to write the output files. When not provided, oputput files are stored in the same folder as that of the input files.",
    )
    cli.add_argument(
        "--weight-name",
        type=str,
        default="weight",
        help="Name of the dataset where weight should be stored.",
    )

    return cli


def copy_cooler(src: pathlib.Path, dest: pathlib.Path, resolution: int, weight_name: str):
    logging.info("Copying %s::/resolutions/%d to %s...", src, resolution, dest)

    group = f"/resolutions/{resolution}"
    with h5py.File(src, "r") as fin:
        if group not in fin:
            raise RuntimeError(f"File {src} does not contain data for {resolution} resolution")
        with h5py.File(dest, "w") as fout:
            for grp in fin[group].keys():
                fin.copy(f"{group}/{grp}", fout, name=grp)
            fout.attrs.update(fin[group].attrs)

            logging.info("Adding balancing weights under %s::/bins/%s", dest, weight_name)
            fout.create_dataset(f"bins/{weight_name}", data=np.ones_like(fout["/bins/chrom"], dtype=float))


def main():
    args = vars(make_cli().parse_args())
    weight_name = args["weight_name"]

    for mcool in args["mcools"]:
        output_dir = args["output_dir"]
        if output_dir is None:
            output_dir = mcool.parent

        for resolution in args["resolutions"]:
            output_path = output_dir / f"{mcool.stem}_{resolution}.cool"
            copy_cooler(mcool, output_path, resolution, weight_name)


def setup_logger(level: str):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(level=level, format=fmt)
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger("INFO")
    main()
