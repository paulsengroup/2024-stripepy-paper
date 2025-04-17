#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import itertools
import json
import logging
import pathlib
import shutil
import subprocess as sp
import tempfile
import textwrap
from typing import Any, Callable, Dict, List


def existing_file(arg):
    if (path := pathlib.Path(arg)).exists():
        return path

    raise FileNotFoundError(f'File "{arg}" does not exists')


def num_cpus(arg: str) -> int:
    import multiprocessing as mp

    try:
        n = int(arg)
        if 0 < n <= mp.cpu_count():
            return n
    except:  # noqa
        pass

    raise argparse.ArgumentTypeError(
        f"Not a valid number of CPU cores (allowed values are integers between 1 and {mp.cpu_count()})"
    )


def positive_int(arg: str) -> int:
    try:
        n = int(arg)
        if n > 0:
            return n
    except:  # noqa
        pass

    raise argparse.ArgumentTypeError("Not a positive integer")


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "matrix-file",
        type=existing_file,
        help="Path to a file in .[m]cool or .hic format to use for benchmarking.",
    )
    cli.add_argument(
        "resolution",
        type=positive_int,
        help="Resolution in bp to be used for benchmarking.",
    )
    cli.add_argument(
        "--tool",
        type=str,
        choices={"stripepy", "chromosight", "stripecaller", "stripenn"},
        required=True,
        nargs="+",
        help="Name of one or more tools to be benchmarked.",
    )
    cli.add_argument(
        "--nproc",
        type=num_cpus,
        default=8,
        help="Maximum number of parallel processes to use (default: %(default)s).",
    )
    cli.add_argument(
        "--iterations",
        type=positive_int,
        default=5,
        help="Number of iterations to run for each tool.",
    )
    cli.add_argument(
        "--suppress-output",
        action="store_true",
        default=False,
        help="Capture output generated on stderr and stdout by the benchmarked tool(s).",
    )

    return cli


@functools.cache
def find_gnu_time() -> pathlib.Path:
    path = shutil.which("time")
    if path is None:
        path = shutil.which("gtime")

        if path is None:
            raise RuntimeError("Unable to find time or gtime executable in your PATH!")

    return pathlib.Path(path).resolve()


def get_gnu_time_args(dest: pathlib.Path | None) -> List[str]:
    fields = textwrap.dedent(
        """
        {"cpu_time_kernel_sec": %S,
        "cpu_time_user_sec": %U,
        "elapsed_real_time_sec": %e,
        "max_rss_kb": %M}
        """
    ).strip()

    args = [str(find_gnu_time()), "--format", fields]

    if dest is not None:
        args.extend(("--output", str(dest)))

    return args


def find_tools(tools: List[str]) -> Dict[str, pathlib.Path]:
    found_tools = {}
    for i, tool in enumerate(tools):
        if tool == "stripecaller":
            path = shutil.which("call-stripes")
        else:
            path = shutil.which(tool)
        if path is None:
            raise RuntimeError(f"Unable to find {tool} executable in your PATH!")
        found_tools[tool] = pathlib.Path(path).resolve()

    return found_tools


def compute_cpu_pct(result: Dict[str, Any]) -> float:
    return 100 * (result["cpu_time_kernel_sec"] + result["cpu_time_user_sec"]) / result["elapsed_real_time_sec"]


def run_tool(args: List[str], tmpdir: pathlib.Path, suppress_output: bool) -> Dict[str, Any]:
    with tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8", dir=tmpdir) as f:
        metric_file = pathlib.Path(f.name)
        time_args = get_gnu_time_args(metric_file)

        args = time_args + args

        if not suppress_output:
            sp.check_call(args)
        else:
            try:
                sp.check_output(args, stderr=sp.STDOUT)
            except sp.CalledProcessError as e:
                logging.error(e.output)
                raise

        f.seek(0)
        res = json.load(f)

    res["cpu_pct"] = compute_cpu_pct(res)
    return res


def run_stripepy(
    stripepy: pathlib.Path,
    matrix_file: pathlib.Path,
    resolution: int,
    nproc: int,
    iteration: int,
    suppress_output: bool,
) -> Dict[str, Any]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        args = [
            str(stripepy),
            "call",
            str(matrix_file),
            str(resolution),
            "--output",
            str(tmpdir / "out.hdf5"),
            "--nproc",
            str(nproc),
        ]

        result = run_tool(args, tmpdir, suppress_output)
        return result | {
            "matrix_file": str(matrix_file.name),
            "resolution": resolution,
            "nproc": nproc,
            "tool": "StripePy",
            "iteration": iteration,
        }


def run_chromosight(
    chromosight: pathlib.Path,
    matrix_file: pathlib.Path,
    resolution: int,
    nproc: int,
    iteration: int,
    suppress_output: bool,
) -> Dict[str, Any]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        args_upper = [
            str(chromosight),
            "detect",
            f"{matrix_file}::/resolutions/{resolution}",
            str(tmpdir / "out1.txt"),
            "--pattern",
            "stripes_left",
            "--min-dist",
            "20000",
            "--max-dist",
            "200000",
            "--threads",
            str(nproc),
        ]

        args_lower = [
            str(chromosight),
            "detect",
            f"{matrix_file}::/resolutions/{resolution}",
            str(tmpdir / "out2.txt"),
            "--pattern",
            "stripes_right",
            "--min-dist",
            "20000",
            "--max-dist",
            "200000",
            "--threads",
            str(nproc),
        ]

        result1 = run_tool(args_upper, tmpdir, suppress_output)
        result2 = run_tool(args_lower, tmpdir, suppress_output)
        result = {}

        for k in result1:
            if k.endswith("_sec"):
                result[k] = result1[k] + result2[k]

        result["max_rss_kb"] = max(result1["max_rss_kb"], result2["max_rss_kb"])
        result["cpu_pct"] = compute_cpu_pct(result1)

        return result | {
            "matrix_file": str(matrix_file.name),
            "resolution": resolution,
            "nproc": nproc,
            "tool": "Chromosight",
            "iteration": iteration,
        }


def run_stripenn(
    stripenn: pathlib.Path,
    matrix_file: pathlib.Path,
    resolution: int,
    nproc: int,
    iteration: int,
    suppress_output: bool,
) -> Dict[str, Any]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        args = [
            str(stripenn),
            "compute",
            "--cool",
            f"{matrix_file}::/resolutions/{resolution}",
            "--out",
            str(tmpdir / "out"),
            "--maxpixel",
            "0.95,0.96,0.97,0.98,0.99",
            "--numcores",
            str(nproc),
        ]

        result = run_tool(args, tmpdir, suppress_output)
        return result | {
            "matrix_file": str(matrix_file.name),
            "resolution": resolution,
            "nproc": nproc,
            "tool": "Stripenn",
            "iteration": iteration,
        }


def run_stripecaller(
    stripecaller: pathlib.Path,
    matrix_file: pathlib.Path,
    resolution: int,
    nproc: int,
    iteration: int,
    suppress_output: bool,
) -> Dict[str, Any]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)

        args = [
            str(stripecaller),
            "--path",
            f"{matrix_file}::/resolutions/{resolution}",
            "--output",
            str(tmpdir / "out.bedpe"),
            "--logFile",
            str(tmpdir / "out.log"),
            "--nproc",
            str(nproc),
        ]

        result = run_tool(args, tmpdir, suppress_output)
        return result | {
            "matrix_file": str(matrix_file.name),
            "resolution": resolution,
            "nproc": nproc,
            "tool": "StripeCaller",
            "iteration": iteration,
        }


def get_runner_fx(tool: str) -> Callable:
    try:
        return globals()[f"run_{tool}"]
    except KeyError as e:
        raise NotImplementedError(f"Unknown tool {tool}") from e


def setup_logger(level: str = "INFO"):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(level=level, format=fmt)
    logging.getLogger().setLevel(level)


class _JSONEncoder(json.JSONEncoder):
    """
    An encoder class that can serialize pathlib.Path objects as JSON.
    """

    def default(self, o: Any):
        if isinstance(o, pathlib.Path):
            return str(o)
        return super().default(o)


def main():
    setup_logger()
    args = vars(make_cli().parse_args())

    tools = find_tools(args["tool"])
    matrix_file = args["matrix-file"]
    resolution = args["resolution"]
    nproc = args["nproc"]
    suppress_output = args["suppress_output"]
    iterations = args["iterations"]

    logging.info("\nPARAMS:\n%s", json.dumps(args, indent=2, sort_keys=True, cls=_JSONEncoder))

    metrics = []
    for i, tool_name in itertools.product(range(1, iterations + 1), tools):
        tool_bin = tools[tool_name]
        logging.info("[%d/%d] running %s (%s)...", i, iterations, tool_name, tool_bin)
        runner = get_runner_fx(tool_name)

        result = runner(
            tool_bin,
            matrix_file,
            resolution,
            nproc,
            i,
            suppress_output,
        )

        logging.info("[%d/%d] running %s took %ss!", i, iterations, tool_name, result["elapsed_real_time_sec"])
        metrics.append(result)

    assert len(metrics) != 0

    print("\t".join(metrics[0].keys()))
    for m in metrics:
        print("\t".join(str(x) for x in m.values()))


if __name__ == "__main__":
    main()
