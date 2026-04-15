#!/bin/bash
set -euo pipefail

json_dir="../Datafiles"
json_file="${json_dir}/datafile.json"
backup_file="${json_dir}/datafile_backup.json"

CPU_THREADS=(1 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64)
MPI_PROCS=(1 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64)

# Lines to extract from runs
GREP_PAT='Mismatch|V&V|Total Time:|Total cycle time'

# --- helpers ---
restore_json() { [[ -f "$backup_file" ]] && mv -f "$backup_file" "$json_file"; }
latest_vnv() {
  shopt -s nullglob
  local files=(vnv*)
  shopt -u nullglob
  [[ ${#files[@]} -eq 0 ]] && return 1
  ls -1t "${files[@]}" | head -n 1
}
run_and_report() {
  # $1 = mode, $2.. = extra args
  local mode="$1"; shift || true
  echo "→ Running: ./run.sh ${mode} ${*:-}"
  ./run.sh "$mode" "$@" 2>&1 | tee -a "logs_${mode}.log" | grep -Ei "$GREP_PAT" || true
  echo "✓ ${mode} ${*:-} finished"
}

# --- backup current json and ensure we restore it at the end ---
cp -f "$json_file" "$backup_file"
trap restore_json EXIT

echo "==============================="
echo "Baseline (SERIAL) with datafile_400_1.json"
echo "==============================="

# Use datafile_400_1.json for baseline + GPU + CPU runs
cp -f "${json_dir}/datafile_400_1.json" "$json_file"

# SERIAL (baseline vnv + time)
run_and_report SERIAL

# Capture baseline vnv artifact once for all following runs
mkdir -p baseline_vnv_400
if vnv_file="$(latest_vnv)"; then
  cp -f "$vnv_file" "baseline_vnv_400/$vnv_file"
  cp -f "baseline_vnv_400/$vnv_file" "baseline_vnv"
  echo "✓ Saved baseline V&V: $vnv_file → baseline_vnv_400/ (and copied to ./baseline_vnv)"
else
  echo "✗ No vnv* output found after SERIAL"; exit 1
fi

# Also print a one-liner baseline time summary
baseline_time=$(grep -E 'Total Time:' -m1 "logs_SERIAL.log" | awk -F'Total Time:' '{gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2}' || true)
[[ -n "${baseline_time:-}" ]] && echo "BASELINE Total Time: ${baseline_time}"

echo
echo "==============================="
echo "GPU run (same datafile_400_1.json)"
echo "==============================="
run_and_report GPU
echo

echo "==============================="
echo "CPU runs (threads = ${CPU_THREADS[*]}) with datafile_400_1.json"
echo "==============================="
for t in "${CPU_THREADS[@]}"; do
  run_and_report CPU "$t"
done
echo

echo "==============================="
echo "MPI runs (procs = ${MPI_PROCS[*]})"
echo "  Using datafile_400_<t>.json per run"
echo "==============================="
for t in "${MPI_PROCS[@]}"; do
  # swap JSON to the per-proc file before each MPI run
  src="${json_dir}/datafile_400_${t}.json"
  if [[ ! -f "$src" ]]; then
    echo "✗ Missing ${src}"; exit 1
  fi
  cp -f "$src" "$json_file"
  run_and_report MPI "$t"
done

echo
echo "✓ Done. JSON restored."

