#!/bin/bash

########################################
# Mode selection
# Usage:
#   ./run_all.sh        -> ALL modes
#   ./run_all.sh MPI    -> only MPI
#   ./run_all.sh CPU    -> only CPU
#   ./run_all.sh SERIAL -> only SERIAL
#   ./run_all.sh GPU    -> only GPU
########################################

mode="${1:-ALL}"
mode=$(echo "$mode" | tr '[:lower:]' '[:upper:]')

echo "run_all.sh mode: $mode"

########################################
# Config
########################################

# Problem sizes to sweep (must match: datafile_<size>_<ranks>.json)
sizes=(400)

# CPU thread counts to sweep (same list as MPI ranks)
cpu_threads=(1 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64)

# MPI ranks to sweep
mpi_procs=(1 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64)

json_dir="../Datafiles"
json_file="${json_dir}/datafile.json"
backup_file="${json_dir}/datafile_backup.json"

########################################
# Backup original JSON (if any)
########################################

if [[ -f "$json_file" ]]; then
    cp "$json_file" "$backup_file"
    echo "Backed up original $json_file to $backup_file"
else
    backup_file=""
    echo "No original $json_file to back up"
fi

for size in "${sizes[@]}"; do
    echo "==============================="
    echo "Running for size $size"
    echo "==============================="

    ########################################
    # For modes other than pure MPI:
    # use the 1-rank JSON as the baseline config
    ########################################
    if [[ "$mode" != "MPI" ]]; then
        base_json="${json_dir}/datafile_${size}_1.json"
        if [[ ! -f "$base_json" ]]; then
            echo "ERROR: Base JSON not found: $base_json"
            [[ -n "$backup_file" ]] && cp "$backup_file" "$json_file"
            exit 1
        fi
        cp "$base_json" "$json_file"
        echo "Using base JSON: $base_json -> $json_file"
    fi

    ########################################
    # SERIAL
    ########################################
    if [[ "$mode" == "ALL" || "$mode" == "SERIAL" ]]; then
        echo "Running: ./run.sh SERIAL"
        serial_output=$(./run.sh SERIAL 2>&1)
        serial_status=$?

        # Only V&V / mismatches / Total cycle time
        echo "$serial_output" | grep -Ei 'Mismatch|V&V|Total cycle time' || true

        if [[ $serial_status -ne 0 ]]; then
            echo "SERIAL FAILED with exit code $serial_status"
            [[ -n "$backup_file" ]] && cp "$backup_file" "$json_file"
            exit 1
        fi

        echo "SERIAL finished"

        # Move and rename the vnv* file
        mkdir -p "baseline_vnv_size_$size"
        vnv_file=$(ls vnv* 2>/dev/null | head -n 1)

        if [[ -n "$vnv_file" ]]; then
            mv "$vnv_file" "baseline_vnv_size_$size/$vnv_file"
            # baseline_vnv is a FILE (baseline reference)
            cp "baseline_vnv_size_$size/$vnv_file" "baseline_vnv"
            echo "Moved $vnv_file to baseline_vnv_size_$size/"
            echo "Copied V&V file to baseline_vnv (baseline file)"
        else
            echo "No vnv* output found after SERIAL"
            [[ -n "$backup_file" ]] && cp "$backup_file" "$json_file"
            exit 1
        fi
    fi

    ########################################
    # CPU: vary number of threads
    # (print only V&V + Total cycle time lines)
    ########################################
    if [[ "$mode" == "ALL" || "$mode" == "CPU" ]]; then
        echo "==============================="
        echo "CPU runs (threads = ${cpu_threads[*]}) for size $size"
        echo "==============================="

        base_json="${json_dir}/datafile_${size}_1.json"
        cp "$base_json" "$json_file"

        for nt in "${cpu_threads[@]}"; do
            echo "Running: ./run.sh CPU $nt"
            cpu_output=$(./run.sh CPU "$nt" 2>&1)
            cpu_status=$?

            # Only V&V / mismatches / Total cycle time
            echo "$cpu_output" | grep -Ei 'Mismatch|V&V|Total cycle time' || true

            if [[ $cpu_status -ne 0 ]]; then
                echo "CPU (threads=$nt) FAILED with exit code $cpu_status"
            else
                echo "CPU (threads=$nt) finished"
            fi
            echo ""
        done
    fi

    ########################################
    # GPU (single run)
    # (print only V&V + Total cycle time lines)
    ########################################
    if [[ "$mode" == "ALL" || "$mode" == "GPU" ]]; then
        base_json="${json_dir}/datafile_${size}_1.json"
        cp "$base_json" "$json_file"

        echo "Running: ./run.sh GPU"
        gpu_output=$(./run.sh GPU 2>&1)
        gpu_status=$?

        echo "$gpu_output" | grep -Ei 'Mismatch|V&V|Total cycle time' || true

        if [[ $gpu_status -ne 0 ]]; then
            echo "GPU run FAILED with exit code $gpu_status"
        else
            echo "GPU finished"
        fi
    fi

    ########################################
    # MPI: vary number of ranks
    # (V&V + min/avg/max Total cycle time only)
    ########################################
    if [[ "$mode" == "ALL" || "$mode" == "MPI" ]]; then
        echo "==============================="
        echo "MPI runs (procs = ${mpi_procs[*]}) for size $size"
        echo "==============================="

        for procs in "${mpi_procs[@]}"; do
            mpi_json="${json_dir}/datafile_${size}_${procs}.json"

            if [[ ! -f "$mpi_json" ]]; then
                echo "WARNING: MPI JSON not found for procs=$procs: $mpi_json"
                echo "Skipping MPI $procs"
                echo ""
                continue
            fi

            cp "$mpi_json" "$json_file"

            echo "Running: ./run.sh MPI $procs"
            mpi_output=$(./run.sh MPI "$procs" 2>&1)
            mpi_status=$?

            # Show only V&V / mismatch info
            echo "$mpi_output" | grep -Ei 'Mismatch|V&V' || true

            # Compute min/avg/max over all Total cycle time values in this run
            stats=$(echo "$mpi_output" | awk '
/Total cycle time:/ {
    v = $NF + 0.0;
    print v;
}' | awk '
BEGIN {min=1e99; max=-1e99; sum=0.0; n=0}
{
    val = $1 + 0.0;
    if (val < min) min = val;
    if (val > max) max = val;
    sum += val;
    n++;
}
END {
    if (n > 0) {
        printf "%.9f %.9f %.9f %d\n", min, sum/n, max, n;
    }
}
')
            if [[ -n "$stats" ]]; then
                read min avg max n <<< "$stats"
                echo "MPI $procs: Total cycle time statistics (N=$n)"
                echo "  min = $min"
                echo "  avg = $avg"
                echo "  max = $max"
            else
                echo "MPI $procs: No Total cycle time lines found"
            fi

            # Do NOT print “finished (FAILED, exit code …)”
            # (V&V lines above already indicate pass/fail)
            echo ""
        done
    fi

    echo ""
done

########################################
# Restore the original JSON if backup exists
########################################
if [[ -n "$backup_file" && -f "$backup_file" ]]; then
    cp "$backup_file" "$json_file"
    echo "JSON file restored from backup"
else
    echo "No original JSON backup to restore"
fi

