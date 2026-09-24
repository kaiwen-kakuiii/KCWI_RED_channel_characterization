#!/bin/zsh
# Gate all three cenwave-8600 nights against keck_kcrm_RH3_8600.fits, with the
# ThArRH3 catalogue (the arc frames are ThAr).  2024-10-29 seeded the template,
# so its gate is self-consistency; the other two are independent.
cd "$(dirname "$0")/.."
for n in 2023-11-12 2024-11-06 2024-10-29; do
  echo "########## $n  $(date +%H:%M:%S)"
  /opt/miniconda3/envs/pypeit/bin/python pypeit_test/gate_rh3_template.py \
      --night "$n" --cenwave 8600 2>&1 | grep -vE "INFO\]"
done
echo "GATES DONE $(date +%H:%M:%S)"
