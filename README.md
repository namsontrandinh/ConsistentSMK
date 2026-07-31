# Knapsack-Constrained Submodular Maximization via Dual Cumulative Sets in Streaming

This repository contains the experimental C++ implementation used in the paper **“Knapsack-Constrained Submodular Maximization via Dual Cumulative Sets in Streaming.”**

The code studies **general submodular maximization under a knapsack constraint (GSMK)** in the streaming setting. It implements the following proposed methods:

- **DCS**: Dual Cumulative Sets, a one-pass streaming algorithm
- **DCMS**: Dual Cumulative Sets Multi-Stream, a multi-pass streaming algorithm

The repository also includes the baseline methods used in the experiments:

- **EDL**
- **TwinGreedy**
- **algo9**
- **algo10**

The implementation supports three benchmark applications:

- **IMK**: Influence Maximization under a Knapsack constraint
- **MCK**: Maximum Weighted Cut under a Knapsack constraint
- **RMK**: Revenue Maximization under a Knapsack constraint

## Repository Structure

The project is organized around a small set of C++ executables built from `src/main.cpp` together with two preprocessing utilities.

### Main executables

After running `make`, the following binaries are generated:

- `maxcut`
- `revenue`
- `ic`
- `preproc`
- `preproc_ic`

### Important files

- `Makefile` — build rules for all executables
- `runIC.sh` — batch commands for IMK experiments on `email.bin` and `fb.bin`
- `runMC.sh` — batch commands for MCK experiments on `ER.bin`
- `runMC_Astro.sh` — batch commands for MCK experiments on `Astro.bin`
- `runRM.sh` — batch commands for RMK experiments on `GrQc.bin`
- `runRM_Hept.sh` — batch commands for RMK experiments on `Hept.bin`

### Input datasets

The repository currently uses preprocessed binary graph files:

- `Astro.bin`
- `ER.bin`
- `fb.bin`
- `email.bin`
- `GrQc.bin`
- `Hept.bin`

## Requirements

To build and run the code, you need:

- GNU `g++` with C++17 support
- GNU `make`
- OpenMP support in your compiler for the `revenue` and `ic` binaries
- A Unix-like environment for running the provided shell scripts

The `Makefile` uses the flags `-std=c++17 -O2 -Wall`, links `pthread`, and enables `-fopenmp` for `revenue` and `ic`.

## Build

Compile all executables with:

```bash
make
```

This builds:

- `maxcut`
- `revenue`
- `ic`
- `preproc`
- `preproc_ic`

To remove compiled files:

```bash
make clean
```

## Input Format and Preprocessing

Experiments are run on binary graph files (`*.bin`).

Two preprocessing tools are provided:

- `preproc` — converts a general edge-list input into a `.bin` graph file
- `preproc_ic` — preprocessing tool for influence maximization graphs; according to the build configuration, it always treats the graph as directed and normalizes incoming weights

If you already have the provided `.bin` files, you can run experiments directly without preprocessing.

## Running Experiments

The executable depends on the application:

- `./ic` for **IMK**
- `./maxcut` for **MCK**
- `./revenue` for **RMK**

The general command format is:

```bash
./<binary> --graph <dataset.bin> --B_factor <budget_ratio> --alg <algorithm_name> [other options] --csv <output.csv>
```

### Common options

- `--graph` — input binary graph file
- `--B_factor` — budget factor
- `--alg` — algorithm name (`dcs`, `dcms`, `edl`, `twin_greedy`, `algo9`, `algo10`)
- `--w` — window parameter used by `dcs` and `dcms`
- `--alpha` — additional parameter used by `dcms`
- `--csv` — output CSV file

---

## 1. Influence Maximization (IMK)

The script `runIC.sh` evaluates:

- `email.bin`
- `fb.bin`

with budget factors from **0.01 to 0.05**.

The parameter settings used in the script are:

- `w = 2` for `dcs`
- `w = 2`, `alpha = 2` for `dcms`

### Example

```bash
./ic --graph email.bin --B_factor 0.03 --alg dcs --w 2 --csv email.csv
./ic --graph email.bin --B_factor 0.03 --alg dcms --w 2 --alpha 2 --csv email.csv
./ic --graph email.bin --B_factor 0.03 --alg edl --csv email.csv
```

### Run with Slurm

```bash
sbatch runIC.sh
```

The provided script requests 32 tasks on the `gpu` partition.

---

## 2. Maximum Cut (MCK)

Two scripts are provided:

- `runMC.sh` for `ER.bin`
- `runMC_Astro.sh` for `Astro.bin`

Both scripts evaluate budget factors from **0.10 to 0.65**, using:

- `w = 2` for `dcs`
- `w = 2`, `alpha = 2` for `dcms`

### Example

```bash
./maxcut --graph ER.bin --B_factor 0.25 --alg dcs --w 2 --csv er.csv
./maxcut --graph ER.bin --B_factor 0.25 --alg dcms --w 2 --alpha 2 --csv er.csv
./maxcut --graph ER.bin --B_factor 0.25 --alg twin_greedy --csv er.csv
```

### Run with Slurm

```bash
sbatch runMC.sh
sbatch runMC_Astro.sh
```

---

## 3. Revenue Maximization (RMK)

Two scripts are provided:

- `runRM.sh` for `GrQc.bin`
- `runRM_Hept.sh` for `Hept.bin`

For `GrQc.bin`, the script evaluates budget factors from **0.10 to 0.65**.

For `Hept.bin`, the script keeps only the larger-budget runs active. In particular:

- the smallest budget settings are commented out
- from `0.35` onward, both `dcs` and `dcms` are included
- in the `0.30` block, only the baseline methods remain active

### Example

```bash
./revenue --graph GrQc.bin --B_factor 0.20 --alg dcs --w 2 --csv GrQc.csv
./revenue --graph GrQc.bin --B_factor 0.20 --alg dcms --w 2 --alpha 2 --csv GrQc.csv
./revenue --graph GrQc.bin --B_factor 0.20 --alg edl --csv GrQc.csv
```

### Run with Slurm

```bash
sbatch runRM.sh
sbatch runRM_Hept.sh
```

## Output

Each run appends its results to the CSV file specified by `--csv`.

The provided scripts generate files such as:

- `email.csv`
- `fb.csv`
- `er.csv`
- `Astro.csv`
- `GrQc.csv`
- `Hept.csv`

## Notes

- This repository is centered on the compiled C++ implementation and binary datasets.
- The provided shell scripts reflect the command structure used in the experiments.
- If you want to use your own datasets, you may need to generate new `.bin` files using `preproc` or `preproc_ic` before running the corresponding executable.

## Citation

If you use this repository, please cite the associated paper:

```bibtex
@article{tran2026dcs,
  title   = {Knapsack-Constrained Submodular Maximization via Dual Cumulative Sets in Streaming},
  author  = {Tran, Uyen T. T. and Ha, Dung T. K. and Pham, Canh V.},
  journal = {to appear},
  year    = {2026}
}
```

You may replace the placeholder citation above with the final bibliographic entry once the paper is publicly available.
