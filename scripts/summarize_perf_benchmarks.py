#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import sys

import pandas as pd


def existing_file(arg):
    if (path := pathlib.Path(arg)).exists():
        return path

    raise FileNotFoundError(f'File "{arg}" does not exists')


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsvs",
        nargs="+",
        type=existing_file,
        help="One or more TSV file to be processed.",
    )

    return cli


def format_timedelta(td: pd.Timedelta) -> str:
    days = td.days
    hours, remainder = divmod(td.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    seconds += (td.microseconds / 1.0e6) % 1

    parts = []
    if days > 0:
        parts.append(f"{days}d")
    if hours > 0:
        parts.append(f"{hours}h")
    if minutes > 0:
        parts.append(f"{minutes}m")
    if seconds > 0:
        if seconds.is_integer():
            parts.append(f"{seconds:.0f}s")
        else:
            parts.append(f"{seconds:.4g}s")

    assert len(parts) > 0

    return " ".join(parts)


def format_value(mean, std: float) -> str:
    is_timedelta = isinstance(mean, pd.Timedelta)

    if is_timedelta:
        mean = format_timedelta(mean)
    elif isinstance(mean, float):
        mean = round(mean, 2)
        if mean.is_integer():
            mean = int(mean)

    if is_timedelta:
        return f"{mean} ± {std:.2f}s"

    return f"{mean} ± {std:.2f}"


def summarize(df_time: pd.DataFrame, df_float: pd.DataFrame) -> pd.DataFrame:
    mean = df_time.groupby(["matrix_file", "resolution", "nproc", "tool"]).mean()
    std = df_float.groupby(["matrix_file", "resolution", "nproc", "tool"]).sem()

    df = mean.copy().astype(object)
    for row in mean.index:
        for col in mean.columns:
            x = mean.loc[row, col]
            y = std.loc[row, col]

            df.loc[row, col] = format_value(x, y)

    df.columns = [col.removesuffix("_sec") for col in df.columns]

    baseline = mean.loc[mean.index.get_level_values("tool") == "StripePy", "elapsed_real_time_sec"].iloc[0]
    df["perf_improvement"] = [f"{round(x, 2)}x" for x in mean["elapsed_real_time_sec"] / baseline]

    return df


def main():
    args = vars(make_cli().parse_args())

    df1 = pd.concat([pd.read_table(tsv) for tsv in args["tsvs"]]).drop(columns=["iteration"])
    df2 = df1.copy()

    for col in df1:
        if col.endswith("_sec"):
            df1[col] = pd.to_timedelta(df1[col], unit="s")

    for df in (df1, df2):
        df["max_rss_kb"] = df["max_rss_kb"] / 1.0e3
        df.rename(columns={"max_rss_kb": "max_rss_mb"}, inplace=True)

    summarize(df1, df2).to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    main()
