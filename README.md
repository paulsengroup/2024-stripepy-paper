<!--
Copyright (C) 2024 Andrea Raffo <andrea.raffo@ibv.uio.no>

SPDX-License-Identifier: MIT
-->

# Synopsis

This repository contains the source code used for:

- extracting individual resolutions from the .mcool files contained in the StripeBench benchmark and adding a vector of 1.0s as weights (required for some stripe callers);
- the data analysis for the StripePy manuscript (preprint available at [doi.org/10.1101/2024.12.20.629789](https://doi.org/10.1101/2024.12.20.629789)).
- the analysis of performance.

Input data consists of:

- The StripeBench benchmark, see [doi.org/10.5281/zenodo.14448329](https://doi.org/10.5281/zenodo.14448329);
- The real contact maps and their ground truth annotations, as detailed in the StripePy manuscript;
- The output of StripePy and other stripe callers over the benchmark and maps mentioned in the previous two points, see [10.5281/zenodo.15308825](https://doi.org/10.5281/zenodo.15308825);
  To reproduce StripePy's output, you can use any version within the 1.x.x series;
- To run the performance benchmarks, the contact map ENCFF993FGR in mcool format should be located inside this repository under `data/ENCFF993FGR.mcool`.

## Requirements

This repository includes an `environment.yml` file, which contains the dependencies for the project environment.

- To create the environment, run the following command: `conda env create -f environment.yml`
- Activate the environment using the command: `conda activate 2024-stripepy-paper`

Once the repository is downloaded, **the user is tasked** with creating a folder named `output` inside it, together with a
set of subfolders which are represented in the following tree:

```txt
output
├── StripeBench
│   ├── RoIs
│   ├── boxplots
│   ├── heatmaps
│   ├── medians
│   └── tables
└── real data
    ├── RoIs
    └── tables
```

On UNIX systems, this can be accomplished with:

```bash
mkdir -p output/StripeBench/{RoIs,boxplots,heatmaps,medians,tables} output/real\ data/{RoIs,tables}
```

Plase note that to run the `scripts/run_perf_benchmarks.sh` script, Docker is required.

## Scripts

The following scripts can be found inside the `scripts/` folder:

- `preprocess_modle_matrix.py` extracts individual resolutions from the .mcool files contained in the StripeBench benchmark and adds a vector of 1.0s as weights
  (required for some stripe callers). Example usage:

  ```bash
  scripts/preprocess_modle_matrix.py StripeBench/data/grch38_h1_rad21_*/*.mcool --resolutions 5000 10000 25000 50000
  ```

- `run_evaluation_StripeBench.py` generates Figure 3B-Q, Supplementary Figures 1-4, and Supplementary Tables 1-6.
- `plot_RoIs_StripeBench.py` generates Figure 3A.
- `run_evaluation_real_data.py` generates Supplementary Table 7.
- `plot_RoIs_real_data.py` generates Supplementary Figure 5.
- `run_perf_benchmarks.sh` generates Supplementary Table 8.
- `compare_normalizations_stripepy.py` generates Supplementary Table 9.
- `plot_RoIs_real_data_normalizations.py` generates Supplementary Figure 6.

To check the input requirements for the various script, activate the environment and then run

```bash
scripts/script_name.py --help
```
