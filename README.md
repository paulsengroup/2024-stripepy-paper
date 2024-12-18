# Synopsis

This repository contains the source code used for the data analysis for the StripePy manuscript (preprint available soon).

Input data consists of:
- The StripeBench benchmark, see https://doi.org/10.5281/zenodo.14448329
- The real contact maps and their ground truth annotations, as detailed in the StripePy manuscript
- The output of StripePy and other stripe callers over the benchmark and maps mentioned in the previous two points, see https://doi.org/10.5281/zenodo.14449731

## Requirements
This repository includes an `environment.yml` file, which contains the dependencies for the project environment.
- To create the environment, run the following command: `conda env create -f environment.yml`
- Activate the environment using the command: `conda activate 2024-stripepy-paper`

Once the repository is downloaded, **the user is tasked** with creating a folder named `output` inside it, together with a
set of subfolders which are represented in the following tree:
```
output
├── StripeBench
│   ├── RoIs
│   ├── boxplots
│   ├── heatmaps
│   ├── medians
│   └── tables
└── real data
    ├── RoIs
    └── tables
```

## Scripts
The following scripts are included:
- `run_evaluation_StripeBench.py` generates Figures 3B-Q and 4, Extended Data Figures 1-3, and Tables 1-5.
- `plot_RoIs_StripeBench.py` generates Figure 1A
- `run_evaluation_real_data.py` generates Table 6
- `plot_RoIs_real_data.py` generates Figure 5
- `compare_normalizations_stripepy.py` generates Table 7.
- `plot_RoIs_real_data_normalizations.py` generates Extended Data Figure 4.

To check the input requirements for the various script, activate the environment and then run
```bash
python script_name.py --help
```
