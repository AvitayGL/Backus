#!/bin/bash

VTUNE=/opt/sw/compilers/oneAPI/2024/vtune/2024.2/bin64/vtune

"$VTUNE" collect hotspots \
    -r vtune_cpu_threads${OMP_NUM_THREADS} \
    -- ../exec/main

