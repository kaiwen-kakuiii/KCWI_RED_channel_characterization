#!/bin/zsh
# pypeit_setup -c all for every RM2 night.  One night per directory, so a
# config letter never has to be stable across nights (RM1.md section 1).
R=/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization
BASE="$R/RM2 pypeit run"
ENV=/opt/miniconda3/envs/pypeit/bin
for night in $(ls "$R/fits/by_night/RM2"); do
  mkdir -p "$BASE/$night/pypeit_run"
  echo "=== $night"
  $ENV/pypeit_setup -s keck_kcrm -r "$R/fits/by_night/RM2/$night" -c all \
      -d "$BASE/$night/pypeit_run" > "$BASE/$night/setup.log" 2>&1
  echo "  rc=$? setups: $(ls -d "$BASE/$night/pypeit_run"/keck_kcrm_* 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"
done
echo "SETUP RM2 COMPLETE"
