#!/bin/zsh
# Round 2: gate the configs the shipped template failed against
# keck_kcrm_RM2_8900.fits (seeded from 2024-01-04 B alone), two streams wide.
# 2025-01-02 is gated separately once round 1 has released its directory.
R=/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization
G="$R/pypeit_test/gate_rm2_template.py"
PY=/opt/miniconda3/envs/pypeit/bin/python
OUT="$R/RM2 pypeit run/gatelogs"; mkdir -p "$OUT"
gate () {
  $PY "$G" --night $1 --setup $2 > "$OUT/$1_$2_own.txt" 2>&1
  echo "$1/$2 rc=$? :: $(grep -E '^(PASS|FAIL)' "$OUT/$1_$2_own.txt" | head -1)"
}
( gate 2024-01-04 keck_kcrm_B; gate 2024-06-11 keck_kcrm_B ) & P1=$!
( gate 2024-04-30 keck_kcrm_B; gate 2025-01-01 keck_kcrm_B ) & P2=$!
wait $P1 $P2
echo "GATE MEDIUM RM2 COMPLETE"
