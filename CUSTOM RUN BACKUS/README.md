# Offloading ScalSALE to GPUs with OpenMP

This project is a work in progress focused on offloading the hydrodynamic computation phase of the ScalSALE code to GPUs using OpenMP 5.0 target offloading.

The code is structured to support three execution modes:
- GPU execution with OpenMP offload
- Multithreaded CPU execution using OpenMP
- Serial CPU execution

## Hydrodynamic Step

The table below lists the main functions involved in the hydrodynamic step, their offload status (updated with each Git push), and relevant implementation notes.

| Function Name                               | Offloaded | Notes                       |
|---------------------------------------------|-----------|-----------------------------|
| `this%Calculate_thermodynamics`             | No        | none                        |
| `time%Calculate_dt`                         | No        | none                        |
| `this%Calculate_acceleration_3d`            | Yes       | none                        |
| `this%Calculate_artificial_viscosity_3d`    | Yes       | none                        |
| `this%Calculate_velocity_3d`                | Yes       | none                        |
| `this%rezone%Calculate_rezone_3d`           | Yes       | none                        |
| `this%Calculate_energy_3d`                  | No        | none                        |
| `this%Calculate_mesh_3d`                    | No        | none                        |

### Additional Notes
blah blah

Note that a function marked as "offloaded" does not necessarily mean the entire function is executed on the GPU. Rather, it indicates that the function has been reviewed and offloaded where possible. Partial offloading may occur depending on data dependencies and parallelism opportunities.

Further details can be found in the "Notes" column of the table above or in dedicated implementation comments.

## How to Run

For the new strong-scaling wrapper, see [RUN_CUSTOM_README.md](/home/shaharm/Desktop/github_repos/Yoni/RUN_CUSTOM_README.md).

-Make sure on NegevHPC you load:

```
ml pgi nvhpc-openmpi3
```

-Compile:
```
./make.sh
```
**Use the provided `run.sh` script to execute the simulation in different modes:**

- Run on GPU (with OpenMP offloading):
  ```bash
  ./run.sh GPU
  ```

- Run on CPU with 32 OpenMP threads:
  ```bash
  ./run.sh CPU
  ```

- Run serially (no OpenMP threading):
  ```bash
  ./run.sh SERIAL
  ```

Note: The current implementation supports only a single MPI rank, which offloads computations to a single GPU. Support for multiple MPI ranks and multi-GPU execution will be added in the future.


## GPU Profiling:

export NVCOMPILER_ACC_NOTIFY=1
export NVCOMPILER_ACC_TIME=1

nsys profile -t cuda,openmp,nvtx --stats=true -o profile_output ./run.sh GPU

  ## V&V
Pay attention to the V&V results, which will be marked as either PASSED (green) or FAILED (red).
The verification compares checksums of the following arrays: ```x, y, z, velocity_x, velocity_y, velocity_z, acceleration_x, acceleration_y, acceleration_z, density, pressure```.

A checksum is considered valid if it differs from the reference by no more than a single bit in the least significant position. This tolerance accounts for minor floating-point variations. Any larger difference results in a failure. 

The V&V is based on a benchmark case with a 360×360×360 Lagrangian mesh.
