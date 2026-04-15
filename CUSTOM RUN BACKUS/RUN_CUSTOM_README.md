# `run_custom` README

This file explains how to use [run_custom](/home/shaharm/Desktop/github_repos/Yoni/run_custom) and the single-permutation helper [run_single_permutation_slurm.py](/home/shaharm/Desktop/github_repos/Yoni/run_single_permutation_slurm.py).

## Overview

The workflow is split into two layers:

- [run_custom](/home/shaharm/Desktop/github_repos/Yoni/run_custom) builds the requested target once and submits all requested permutations to Slurm.
- [run_custom](/home/shaharm/Desktop/github_repos/Yoni/run_custom) builds the requested target once and submits the requested configuration to Slurm.
- [run_single_permutation_slurm.py](/home/shaharm/Desktop/github_repos/Yoni/run_single_permutation_slurm.py) handles one scalar permutation: one mode, one MPI layout, and one grid size.

The base [src/Datafiles/datafile.json](/home/shaharm/Desktop/github_repos/Yoni/src/Datafiles/datafile.json) is not modified in place. Each submitted job gets its own copied `datafile.json` inside its output directory.

## Main Command

```bash
./run_custom -mode [serial|gpu] -vendor [nvidia|intel] -physics [lagrange|euler] -nx N -ny N -nz N <-npx N> <-npy N> <-npz N> <-max_cycles N> <-check> <-clean> <-profile_ncu> <-profile_nsys> <-acc_notify> <-additional_flags_serial "(...)"> <-additional_flags_gpu "(...)"> <-gpu_request_mode [none|gpus|gpus-per-node|gres]> <-gpu_gres_name NAME> <-partition NAME> <-account NAME> <-output_dir PATH> <-main_dir PATH>
```

Notation:

- `[]` = choose one mandatory value
- `<>` = optional

## Important Shell Syntax

Only the additional compilation-flag arguments use parenthesized lists.

Because `bash` treats parentheses specially, quote the entire flag list expression.

Examples:

```bash
-additional_flags_serial "(-Mnovect -g)"
-additional_flags_gpu "(-gpu=cc80 -Minfo=accel)"
```

## Modes

- `serial`
  Uses the CPU executable.
  Fixed thread count:
  - `1`

- `gpu`
  Uses the GPU executable.
  Fixed thread count:
  - `1`

There is no `cpugpu` mode anymore. Run `serial` and `gpu` separately.

## Parameters

- `-mode [serial|gpu]`
  Mandatory.
  Selects the executable and runtime mode.

- `-vendor [nvidia|intel]`
  Mandatory.
  Selects the compiler family for the build.

- `-physics [lagrange|euler]`
  Mandatory.
  Sets `rezone_advect.rezone_type` in each job-local datafile.

- `-npx N`
- `-npy N`
- `-npz N`
  Optional.
  Scalar MPI decomposition values.
  Defaults:
  - `-npx 1`
  - `-npy 1`
  - `-npz 1`

  For each submitted job:
  - `parallel.npx = npx`
  - `parallel.npy = npy`
  - `parallel.npz = npz`
  - `parallel.np = npx * npy * npz`

- `-nx N`
- `-ny N`
- `-nz N`
  Mandatory.
  Scalar grid dimensions.

  For each submitted job:
  - `number_cells_i = nx`
  - `number_cells_j = ny`
  - `number_cells_k = nz`

- `-max_cycles N`
  Optional.
  Sets `simulation_parameters.max_ncyc = N` in each job-local datafile.

- `-check`
  Optional.
  Enables the correctness/debug-oriented compile flags from [src/CMakeLists.txt](/home/shaharm/Desktop/github_repos/Yoni/src/CMakeLists.txt).

- `-clean`
  Optional.
  Forces a clean rebuild through [src/Scripts/clean.sh](/home/shaharm/Desktop/github_repos/Yoni/src/Scripts/clean.sh).
  If omitted, the wrapper prefers [src/Scripts/make.sh](/home/shaharm/Desktop/github_repos/Yoni/src/Scripts/make.sh).

- `-profile_ncu`
- `-profile_nsys`
- `-acc_notify`
  Optional.
  NVIDIA-only GPU runtime options.
  They are valid only with:
  - `-mode gpu`
  - `-vendor nvidia`

- `-additional_flags_serial "(...)"`
  Optional.
  Extra compilation flags appended only to the serial build.

- `-additional_flags_gpu "(...)"`
  Optional.
  Extra compilation flags appended only to the GPU build.

- `-gpu_request_mode [none|gpus|gpus-per-node|gres]`
  Optional.
  Controls how GPU resources are requested from Slurm.
  Default:
  - `none`

  This means no explicit GPU resource line is emitted and the GPU partition is relied on instead.

- `-gpu_gres_name NAME`
  Optional.
  Used only when:
  - `-gpu_request_mode gres`

- `-partition NAME`
  Optional.
  Slurm partition override for the selected mode.
  Defaults:
  - `serial` -> `cluster`
  - `gpu` -> `gpua100`

- `-account NAME`
  Optional.
  Slurm account override for the selected mode.
  Default:
  - `rcl`

- `-output_dir PATH`
  Optional.
  Output root directory.
  Default:
  - `outputs` under `-main_dir`

- `-main_dir PATH`
  Optional.
  Absolute or relative path to the main project folder.
  Default:
  - `.`

  This is the directory that contains:
  - `src/`
  - `run_custom`
  - `run_single_permutation_slurm.py`

## Build Behavior

The wrapper builds only the requested target for the selected mode:

- `serial` -> `main_cpu`
- `gpu` -> `main_gpu`

It rebuilds only when needed. If vendor, `-check`, or additional compilation flags changed since the last build of that target, it automatically switches to a clean rebuild.

Before submitting a new GPU job, `run_custom` checks the current user's active Slurm jobs. If one exists, it submits the new GPU job with an `afterany` dependency on the latest active job ID. Serial jobs are submitted without this dependency.

The stored build-state file is:

- [.run_custom_build_state.json](/home/shaharm/Desktop/github_repos/Yoni/.run_custom_build_state.json)

## Submission Logic

The wrapper submits one Slurm job per invocation for:

- one scalar MPI layout from `-npx/-npy/-npz`
- one scalar grid from `-nx/-ny/-nz`

The thread count is not part of the permutation:

- `serial` always uses `threads = 1`
- `gpu` always uses `threads = 1`

## Output Layout

Each invocation now creates the job directory itself directly under the chosen output root, for example:

```text
outputs/gpu_nvidia_euler_threads1_mpi1x1x1_grid20x20x20/
```

The directory name includes:

- mode
- vendor
- physics (`lagrange` or `euler`)
- thread count
- `npx`, `npy`, `npz`
- `nx`, `ny`, `nz`

If the same directory name already exists, `run_custom` appends `_2`, `_3`, and so on.

Each directory contains:

- `build.log`
- `manifest.json`
- `datafile.json`
- `submit.sbatch`
- `submission.json`
- `solver.log`
- `slurm-<jobid>.out`
- `slurm-<jobid>.err`

Example subdirectory names:

- `serial_nvidia_euler_threads1_mpi2x1x1_grid40x40x40/`
- `gpu_nvidia_euler_threads1_mpi2x2x1_grid40x40x40/`

## Examples

Serial:

```bash
./run_custom -mode serial -vendor nvidia -physics euler -npx 2 -npy 1 -npz 1 -nx 64 -ny 64 -nz 64
```

GPU:

```bash
./run_custom -mode gpu -vendor nvidia -physics lagrange -npx 1 -npy 1 -npz 1 -nx 20 -ny 20 -nz 20 -profile_ncu
```

Different project root and output root:

```bash
./run_custom -mode serial -vendor nvidia -physics euler -npx 1 -npy 1 -npz 1 -nx 40 -ny 40 -nz 40 -main_dir /absolute/path/to/project -output_dir /absolute/path/to/results
```

## Single-Permutation Helper

[run_single_permutation_slurm.py](/home/shaharm/Desktop/github_repos/Yoni/run_single_permutation_slurm.py) is the internal helper used by `run_custom`.

It accepts one scalar MPI layout and one scalar grid size:

```bash
./run_single_permutation_slurm.py -mode [serial|gpu] -vendor [nvidia|intel] -physics [lagrange|euler] -grid_size "(nx ny nz)" -threads N -process_layout "(npx npy npz)" <-max_cycles N> <-gpu_request_mode [none|gpus|gpus-per-node|gres]> <-gpu_gres_name NAME> <-partition NAME> <-account NAME> <-profile_ncu> <-profile_nsys> <-acc_notify> -job_name NAME -output_dir PATH -base_datafile PATH -executable_path PATH
```

Normally you do not need to call it directly.

Slurm allocation details:

- `-npx * npy * npz` defines the total MPI rank count.
- The generated `submit.sbatch` uses `#SBATCH -n total_ranks`.
- The generated `submit.sbatch` also uses `#SBATCH -N total_ranks`, so each MPI rank maps to one node.

## Help

The full CLI help is available with:

```bash
./run_custom --help
```

and also:

```bash
./run_custom -help
```
