#!/usr/bin/env python3

import argparse
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple


SERIAL_FIXED_THREADS = 1
GPU_FIXED_THREADS = 1


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError(f"expected a positive integer, got {value}")
    return number


def parse_grid(text: str) -> Tuple[int, int, int]:
    match = re.fullmatch(r"\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", text.strip())
    if not match:
        raise argparse.ArgumentTypeError(
            "grid size must look like '(i j k)', for example '(20 20 20)'"
        )
    return tuple(int(group) for group in match.groups())


def parse_triplet(text: str, label: str) -> Tuple[int, int, int]:
    match = re.fullmatch(r"\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", text.strip())
    if not match:
        raise argparse.ArgumentTypeError(
            "{0} must look like '(a b c)', for example '(2 1 1)'".format(label)
        )
    values = tuple(int(group) for group in match.groups())
    if any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("{0} values must be positive integers".format(label))
    return values


def module_lines(vendor: str) -> List[str]:
    if vendor == "nvidia":
        return [
            "source /etc/profile.d/lmod.sh >/dev/null 2>&1 || true",
            "ml purge",
            "ml pgi mlog nvhpc-hpcx",
        ]
    return [
        "source /etc/profile.d/lmod.sh >/dev/null 2>&1 || true",
        "ml purge",
        "ml intel/18.0.1.163 openmpi/4.0.4_intel mpi/impi-intel2018",
    ]


def patch_datafile(
    base_datafile: Path,
    output_datafile: Path,
    grid: Tuple[int, int, int],
    process_layout: Tuple[int, int, int],
    threads: int,
    physics: str,
    max_cycles: Optional[int],
) -> None:
    with base_datafile.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    data.setdefault("layers_materials", {})
    data["layers_materials"]["number_cells_i"] = [grid[0]]
    data["layers_materials"]["number_cells_j"] = [grid[1]]
    data["layers_materials"]["number_cells_k"] = [grid[2]]

    data.setdefault("parallel", {})
    total_processes = process_layout[0] * process_layout[1] * process_layout[2]
    data["parallel"]["np"] = total_processes
    data["parallel"]["npx"] = process_layout[0]
    data["parallel"]["npy"] = process_layout[1]
    data["parallel"]["npz"] = process_layout[2]
    data["parallel"]["threads"] = threads

    data.setdefault("rezone_advect", {})
    data["rezone_advect"]["rezone_type"] = physics

    if max_cycles is not None:
        data.setdefault("simulation_parameters", {})
        data["simulation_parameters"]["max_ncyc"] = max_cycles

    with output_datafile.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=4)
        handle.write("\n")


def build_job_script(args: argparse.Namespace, job_dir: Path, datafile_path: Path) -> str:
    solver_log = job_dir / "solver.log"
    profile_prefix = job_dir / "profile"
    total_processes = args.process_layout[0] * args.process_layout[1] * args.process_layout[2]
    effective_cpus = args.threads
    requested_gpus = total_processes if args.mode == "gpu" else 0

    command_parts = [
        "mpirun",
        "-n",
        str(total_processes),
        shlex.quote(str(args.executable_path)),
        shlex.quote(str(datafile_path)),
    ]
    solver_command = " ".join(command_parts)

    if args.profile_ncu:
        solver_command = (
            "ncu --force-overwrite true "
            f"-o {shlex.quote(str(profile_prefix))} "
            f"{solver_command}"
        )
    elif args.profile_nsys:
        solver_command = (
            "nsys profile -t cuda,openmp,nvtx --stats=true --force-overwrite true "
            f"-o {shlex.quote(str(profile_prefix))} "
            f"{solver_command}"
        )

    sbatch_lines = [
        "#!/bin/bash",
        "#SBATCH -J {0}".format(args.job_name),
        "#SBATCH -N {0}".format(total_processes),
        "#SBATCH -n {0}".format(total_processes),
        "#SBATCH --cpus-per-task={0}".format(effective_cpus),
        "#SBATCH --output={0}".format(job_dir / "slurm-%j.out"),
        "#SBATCH --error={0}".format(job_dir / "slurm-%j.err"),
    ]

    if args.partition:
        sbatch_lines.insert(2, "#SBATCH -p {0}".format(args.partition))
    if args.account:
        insert_index = 3 if args.partition else 2
        sbatch_lines.insert(insert_index, "#SBATCH -A {0}".format(args.account))

    if args.mode == "gpu":
        if args.gpu_request_mode == "gpus":
            sbatch_lines.append("#SBATCH --gpus={0}".format(requested_gpus))
        elif args.gpu_request_mode == "gpus-per-node":
            sbatch_lines.append("#SBATCH --gpus-per-node={0}".format(requested_gpus))
        elif args.gpu_request_mode == "gres":
            sbatch_lines.append(
                "#SBATCH --gres={0}:{1}".format(args.gpu_gres_name, requested_gpus)
            )

    sbatch_lines.extend(
        [
            "",
            "set -euo pipefail",
            *module_lines(args.vendor),
            "",
            "cd {0}".format(shlex.quote(str(job_dir))),
            "export OMP_NUM_THREADS={0}".format(args.threads),
            "export SCALSALE_DATAFILE={0}".format(shlex.quote(str(datafile_path))),
        ]
    )

    if args.mode == "gpu":
        sbatch_lines.append("export OMP_TARGET_OFFLOAD=MANDATORY")
        if args.acc_notify:
            sbatch_lines.append("export NVCOMPILER_ACC_NOTIFY=3")

    sbatch_lines.extend(
        [
            "",
            "{",
            "  echo '=== Run Metadata ==='",
            "  echo 'Mode: {0}'".format(args.mode),
            "  echo 'Vendor: {0}'".format(args.vendor),
            "  echo 'Physics: {0}'".format(args.physics),
            "  echo 'Grid: {0}x{1}x{2}'".format(*args.grid),
            "  echo 'Partition: {0}'".format(args.partition if args.partition else "cluster-default"),
            "  echo 'Account: {0}'".format(args.account if args.account else "cluster-default"),
            "  echo 'Processes: {0}'".format(total_processes),
            "  echo 'Nodes: {0}'".format(total_processes),
            "  echo 'MPI Layout: {0}x{1}x{2}'".format(*args.process_layout),
            "  echo 'Threads: {0}'".format(args.threads),
            "  echo 'effective_cpus_per_task: {0}'".format(effective_cpus),
            "  echo 'num_gpus: {0}'".format(requested_gpus),
            "  echo 'gpu_request_mode: {0}'".format(args.gpu_request_mode if args.mode == "gpu" else "n/a"),
            "  echo 'max_cycles: {0}'".format(
                args.max_cycles if args.max_cycles is not None else "keep"
            ),
            "  echo 'datafile: {0}'".format(shlex.quote(str(datafile_path))),
            "  echo 'command: {0}'".format(solver_command),
            "  echo '=== Solver Output ==='",
            "  {0}".format(solver_command),
            "}} > {0} 2>&1".format(shlex.quote(str(solver_log))),
            "",
        ]
    )

    return "\n".join(sbatch_lines)


def parse_args() -> argparse.Namespace:
    description = """Submit one ScalSALE permutation to Slurm.

This helper is intentionally scalar: one mode, one fixed thread count, one MPI layout,
and one grid size. The main ./run_custom script calls it once per permutation.

Default GPU request mode is 'none', which emits no Slurm GPU resource line and relies on
the selected GPU partition instead.
GPU jobs also force threads=1 regardless of the passed -threads value.
"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=description,
    )
    parser.add_argument("-mode", choices=["serial", "gpu"], required=True)
    parser.add_argument("-vendor", choices=["nvidia", "intel"], required=True)
    parser.add_argument("-physics", choices=["lagrange", "euler"], required=True)
    parser.add_argument("-grid_size", type=parse_grid, required=True)
    parser.add_argument("-threads", type=positive_int, required=True)
    parser.add_argument(
        "-process_layout",
        type=lambda text: parse_triplet(text, "process layout"),
        required=True,
        help="MPI layout written as '(npx npy npz)'",
    )
    parser.add_argument("-max_cycles", type=positive_int)
    parser.add_argument("-partition", help="Slurm partition name; if omitted, sbatch uses the cluster default")
    parser.add_argument("-account", help="Slurm account name; if omitted, sbatch uses the cluster default")
    parser.add_argument(
        "-gpu_request_mode",
        choices=["none", "gpus", "gpus-per-node", "gres"],
        default="none",
        help="Slurm GPU request syntax; default omits GPU resource flags and relies on the selected partition",
    )
    parser.add_argument(
        "-gpu_gres_name",
        default="gpu",
        help="GRES resource name used only with -gpu_request_mode gres",
    )
    parser.add_argument("-profile_ncu", action="store_true")
    parser.add_argument("-profile_nsys", action="store_true")
    parser.add_argument("-acc_notify", action="store_true")
    parser.add_argument(
        "-dependency",
        help="optional sbatch dependency string, for example 'afterany:12345'",
    )
    parser.add_argument("-job_name", required=True)
    parser.add_argument("-output_dir", type=Path, required=True)
    parser.add_argument("-base_datafile", type=Path, required=True)
    parser.add_argument("-executable_path", type=Path, required=True)
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.profile_ncu and args.profile_nsys:
        raise SystemExit("use only one of -profile_ncu or -profile_nsys")

    if args.mode == "serial" and (args.profile_ncu or args.profile_nsys or args.acc_notify):
        raise SystemExit(
            "serial jobs do not support -profile_ncu, -profile_nsys, or -acc_notify"
        )

    if args.vendor != "nvidia" and (args.profile_ncu or args.profile_nsys or args.acc_notify):
        raise SystemExit(
            "-profile_ncu, -profile_nsys, and -acc_notify are supported only with -vendor nvidia"
        )


def extract_job_id(stdout_text: str) -> Optional[str]:
    match = re.search(r"Submitted batch job (\d+)", stdout_text)
    if match:
        return match.group(1)
    return None


def main() -> int:
    args = parse_args()
    args.grid = args.grid_size
    if args.mode == "gpu":
        args.threads = GPU_FIXED_THREADS
    else:
        args.threads = SERIAL_FIXED_THREADS
    validate_args(args)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    datafile_path = args.output_dir / "datafile.json"
    patch_datafile(
        base_datafile=args.base_datafile,
        output_datafile=datafile_path,
        grid=args.grid,
        process_layout=args.process_layout,
        threads=args.threads,
        physics=args.physics,
        max_cycles=args.max_cycles,
    )

    sbatch_path = args.output_dir / "submit.sbatch"
    sbatch_path.write_text(
        build_job_script(args, args.output_dir, datafile_path),
        encoding="utf-8",
    )
    sbatch_path.chmod(0o755)

    sbatch_command = ["sbatch"]
    if args.dependency:
        sbatch_command.extend(["--dependency", args.dependency])
    sbatch_command.append(str(sbatch_path))

    submission = subprocess.run(
        sbatch_command,
        cwd=args.output_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=False,
    )

    job_id = extract_job_id(submission.stdout)
    submission_record = {
        "command": sys.argv[1:],
        "returncode": submission.returncode,
        "stdout": submission.stdout,
        "stderr": submission.stderr,
        "job_id": job_id,
        "job_name": args.job_name,
        "mode": args.mode,
        "vendor": args.vendor,
        "physics": args.physics,
        "grid": list(args.grid),
        "process_layout": list(args.process_layout),
        "processes": args.process_layout[0] * args.process_layout[1] * args.process_layout[2],
        "npx": args.process_layout[0],
        "npy": args.process_layout[1],
        "npz": args.process_layout[2],
        "threads": args.threads,
        "num_cpus": args.threads,
        "num_gpus": args.process_layout[0] * args.process_layout[1] * args.process_layout[2]
        if args.mode == "gpu"
        else 0,
        "max_cycles": args.max_cycles,
        "dependency": args.dependency,
        "sbatch_command": sbatch_command,
    }
    (args.output_dir / "submission.json").write_text(
        json.dumps(submission_record, indent=2) + "\n",
        encoding="utf-8",
    )

    if submission.returncode != 0:
        sys.stderr.write(submission.stderr)
        return submission.returncode

    sys.stdout.write(submission.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
