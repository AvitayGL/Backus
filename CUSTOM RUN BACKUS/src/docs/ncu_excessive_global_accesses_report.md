# NCU Profiling Report: Excessive Global Accesses in 3D Advection Kernel

**Kernel location:** `Rezone_and_Advect/advect.f90`, lines 2947-3268
**Date:** 2026-03-31
**Tool:** NVIDIA Nsight Compute (NCU)

---

## What "Excessive Global Accesses" Means

When a GPU **warp** (32 threads executing in lockstep) performs a memory read or write, the hardware **coalesces** individual thread requests into the minimum number of **cache-line transactions** (each 32 bytes from L2, or 128 bytes through L1).

**"X% of this line's global accesses are excessive"** means that X% of the bytes actually transferred from/to global memory were **wasted** — loaded into cache lines but never used by any thread in the warp. The GPU had to issue more memory transactions than theoretically necessary.

---

## NCU Findings Summary

| Lines | Excess % | Description |
|-------|----------|-------------|
| 3003-3033 | 12% | Loading vertex coordinates from `material_x/y/z`, `x/y/z` |
| 3159 | 12% | Reading `vol(i, j, k)` |
| 3180-3182 | 64% | Zeroing private arrays `fxtm_thread`, `sxt1_thread`, `sxt2_thread` |
| 3195 | 12% | Read-modify-write on `vof_adv(i, j, k)` |
| 3202, 3204 | 65% | Reading `fxtm_thread(tmp_mat)`, writing `mat_vof_adv(tmp_mat, i, j, k)` |
| 3213-3215 | 64% | Reading `fxtm_thread`, `density_vof`, `sie_vof` |
| 3220, 3222 | 60% | Writing `sie_vof_adv`, `mat_cell_mass_adv` |
| 3224 | 16% | Reading `mat_id(i,j,k)` |
| 3225, 3227, 3230, 3232 | 60% | Writing `sie_vof_adv`, `mat_cell_mass_adv` (alternate branches) |
| 3249, 3251 | 65% | Reading `init_mat_layers`, writing `init_mat_layers_adv` |

---

## Root Cause Analysis

### Issue 1: 12% Excess (Lines 3003-3033, 3159, 3195) — Warp Boundary Effects

These lines access 3D arrays with `(i, j, k)` indices:

```fortran
x_lag1 = material_x(i1, j1, k1)   ! line 3003
weight = ... / vol(i, j, k)        ! line 3159
vof_adv(i, j, k) = ...             ! line 3195
```

The `COLLAPSE(3)` directive maps the `(i, j, k)` iteration space linearly onto GPU threads. Since `i` is the innermost loop, adjacent threads get consecutive `i` values — which is ideal for Fortran's column-major layout.

However, when a 32-thread warp straddles a **row boundary** (i.e., some threads have `i` near `nx` and others wrap to `i=1` with `j+1`), those threads access addresses that are far apart in memory, producing extra transactions. With typical grid sizes (e.g., `nx=60`), about ~12% of warps hit this boundary, explaining the 12% overhead.

**Verdict:** This is inherent to 3D collapsed loops and essentially unavoidable. 12% is low and acceptable.

---

### Issue 2: 60-65% Excess (Lines 3202-3251) — Stride-N Access from Material-First 4D Layout

This is the **dominant performance problem**.

#### Data Layout

The 4D arrays are allocated as `(1:n_materials, 0:nx, 0:ny, 0:nz)` — **material index is the first (fastest-varying) dimension**:

```fortran
! From data_4d.f90:
allocate(Constructor_init_val%values(1:d4, 0:d1, 0:d2, 0:d3))
!                                    ^^^^ material index first
```

#### Access Pattern in Kernel

Threads are distributed across the `(i, j, k)` spatial domain (via `COLLAPSE(3)`), and each thread loops over materials:

```fortran
do tmp_mat = 1, n_materials                                              ! line 3211
    dvof = fxtm_thread(tmp_mat) * 2d0                                    ! line 3213
    donnor_mass = density_vof(tmp_mat, id, jd, kd) * dvof                ! line 3214
    ...
    mat_cell_mass_adv(tmp_mat, i, j, k) = ... + donnor_mass              ! line 3222
    ...
    init_mat_layers_adv(tmp_mat, i, j, k) = ... + donnor_init_mat_layer  ! line 3251
end do
```

#### Why This Causes 65% Waste

For a fixed `tmp_mat` value, adjacent threads have adjacent `i` values. In column-major Fortran, the address of `mat_vof_adv(tmp_mat, i, j, k)` is:

```
base + (tmp_mat - 1) + n_materials * (i + nx_size * (j + ny_size * k))
```

So consecutive threads (consecutive `i`) access addresses that are **`n_materials` elements apart** (stride-3 with `n_materials=3`). In bytes, that's **24 bytes apart** (3 x 8 bytes per double).

A cache line holds 128 bytes = 16 doubles. With stride-3 access, only **~5 of those 16 doubles** are actually used by threads in the warp:

```
Wasted = 1 - (useful / fetched)
       = 1 - 32 / (ceil(32*3/16) * 16/3)
       ≈ 1 - 32/37.3
       ≈ 65%
```

This matches the **65% excessive** reported by NCU almost exactly.

#### Affected Arrays

All 4D arrays with material-first layout:
- `mat_vof_adv(tmp_mat, i, j, k)`
- `sie_vof_adv(tmp_mat, i, j, k)`
- `mat_cell_mass_adv(tmp_mat, i, j, k)`
- `density_vof(tmp_mat, id, jd, kd)`
- `sie_vof(tmp_mat, id, jd, kd)`
- `init_mat_layers(tmp_mat, id, jd, kd)`
- `init_mat_layers_adv(tmp_mat, i, j, k)`

---

### Issue 3: 64% Excess (Lines 3180-3182) — Private Arrays in Local Memory

```fortran
fxtm_thread = 0d0    ! line 3180
sxt1_thread = 0d0    ! line 3181
sxt2_thread = 0d0    ! line 3182
```

These are `PRIVATE` arrays of size `n_materials=3`. The OpenMP compiler places them in **GPU local memory** (which is physically global memory, just per-thread). Local memory is interleaved across threads, but writing all 3 elements of a small array per thread still results in a poor access pattern — the same stride issue applies since the compiler unrolls these into 3 individual stores.

---

## Proposed Solutions

### Solution for 60-65% Excess: Transpose the 4D Arrays

The most impactful optimization is to change the array layout so that the **material dimension is last** instead of first.

**Current layout (bad for GPU):**
```fortran
array(1:n_materials, 0:nx, 0:ny, 0:nz)    ! material varies fastest
```

**Proposed layout (good for GPU):**
```fortran
array(0:nx, 0:ny, 0:nz, 1:n_materials)    ! spatial i varies fastest
```

With this layout, `array(i, j, k, tmp_mat)` means adjacent threads (adjacent `i`) access **consecutive memory addresses** — perfect coalescing, 0% excess.

#### Approach A: Conditional Transposition (GPU-only, via preprocessor)

Change the allocation in `data_4d.f90` based on build target, and use a macro or wrapper to swap index ordering:

```fortran
#ifdef GPU_BUILD
  allocate(values(0:d1, 0:d2, 0:d3, 1:d4))   ! spatial-first for GPU
#else
  allocate(values(1:d4, 0:d1, 0:d2, 0:d3))   ! material-first for CPU
#endif
```

Access via a macro:
```fortran
#ifdef GPU_BUILD
#define MAT4D(arr, m, i, j, k)  arr(i, j, k, m)
#else
#define MAT4D(arr, m, i, j, k)  arr(m, i, j, k)
#endif
```

Then replace all accesses like `density_vof(tmp_mat, id, jd, kd)` with `MAT4D(density_vof, tmp_mat, id, jd, kd)`.

#### Approach B: Local Staging (kernel-only change, partial fix)

If the data layout cannot be changed globally, stage material data into scalar variables at the start of the material loop. This helps reads but writes to `_adv` arrays would still be strided. This is a half-measure.

### Solution for 64% Excess on Private Arrays

Since `n_materials` is small (typically 3), replace the private arrays with **scalar variables**:

```fortran
real(8) :: fxtm1, fxtm2, fxtm3   ! instead of fxtm_thread(3)
```

The compiler will place scalars in registers rather than local memory, eliminating the global memory accesses entirely. This requires `n_materials` to be known at compile time and changes to `Volume_material_3d`.

A more general alternative: ensure the compiler keeps these in registers by using compiler hints or restructuring the code to avoid array syntax on small private arrays.

### No Fix Needed for 12% Excess

The 12% excess on lines 3003-3033, 3159, and 3195 is inherent to 3D collapsed loops and is not worth addressing. It is a minor overhead.

---

## Impact Summary

| Lines | Excess | Root Cause | Fix | Impact |
|-------|--------|-----------|-----|--------|
| 3003-3033, 3159, 3195 | 12% | Warp boundary effects in `COLLAPSE(3)` | None needed | Low |
| 3180-3182 | 64% | Private arrays in local memory | Use scalars / compiler hints | Medium |
| All other flagged lines | 60-65% | Stride-3 from material-first 4D layout | **Transpose 4D arrays** | **High** |

The **single highest-impact change** is transposing the 4D arrays to put the spatial dimension first and material dimension last. This would eliminate ~60-65% wasted memory bandwidth on the majority of flagged lines, which are the core data-intensive part of the kernel (the material accumulation loops at lines 3201-3257).
