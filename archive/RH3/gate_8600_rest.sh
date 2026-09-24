#!/bin/zsh
# Gate the remaining cenwave-8600 nights against keck_kcrm_RH3_8600.fits.
# 2024-10-29 is the night the template was seeded from, so its gate is a
# self-consistency check; 2024-11-06 is independent.  gate_rh3_template.py says
# which is which in its own output.
cd "$(dirname "$0")/.."
for n in 2024-11-06 2024-10-29; do
  echo "########## $n  $(date +%H:%M:%S)"
  /opt/miniconda3/envs/pypeit/bin/python pypeit_test/gate_rh3_template.py \
      --night "$n" --cenwave 8600 2>&1 | grep -vE "INFO\]"
done
echo "GATES DONE $(date +%H:%M:%S)"
