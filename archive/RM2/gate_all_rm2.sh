#!/bin/zsh
# Gate the shipped keck_kcrm_RM2.fits on every RM2 science config, three streams
# wide.  Calibrations only (-c), so nothing here touches science frames.
R=/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization
G="$R/pypeit_test/gate_rm2_template.py"
PY=/opt/miniconda3/envs/pypeit/bin/python
OUT="$R/RM2 pypeit run/gatelogs"; mkdir -p "$OUT"
gate () {  # night setup
  $PY "$G" --night $1 --setup $2 > "$OUT/$1_$2.txt" 2>&1
  echo "$1/$2 rc=$? :: $(grep -E '^(PASS|FAIL)' "$OUT/$1_$2.txt" | head -1)"
}
streamA () { gate 2024-01-04 keck_kcrm_B; gate 2024-04-30 keck_kcrm_B; }
streamB () { gate 2024-06-11 keck_kcrm_B; gate 2025-01-01 keck_kcrm_B; gate 2025-01-02 keck_kcrm_B; }
streamC () { gate 2024-06-10 keck_kcrm_A; gate 2024-12-24 keck_kcrm_A; gate 2024-05-09 keck_kcrm_A; }
streamA & A=$!
streamB & Bp=$!
streamC & C=$!
wait $A $Bp $C
echo "GATE ALL RM2 COMPLETE"
