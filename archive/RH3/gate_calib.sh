#!/bin/zsh
# Calibration-only pass per night, stock settings: keck_kcwi.py picks
# keck_kcrm_RH3.fits automatically for dispname RH3.
ENV=/opt/miniconda3/envs/pypeit/bin
BASE="$(cd "$(dirname "$0")" && pwd)"
for n in "$@"; do
  d="$BASE/$n/pypeit_run/keck_kcrm_B"
  echo "=== $n keck_kcrm_B  $(date +%H:%M:%S) ==="
  (cd "$d" && $ENV/run_pypeit keck_kcrm_B.pypeit -c > calib.log 2>&1; echo "rc=$?")
  wc=$(ls "$d"/Calibrations/WaveCalib_*.fits 2>/dev/null | head -1)
  if [ -n "$wc" ]; then
    $ENV/python "$BASE/../pypeit_test/pick_seed_slits.py" "$wc"
  else
    echo "  NO WaveCalib written"
  fi
done
echo "GATE DONE $(date +%H:%M:%S)"
