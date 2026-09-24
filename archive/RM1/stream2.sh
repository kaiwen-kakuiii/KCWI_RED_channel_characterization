#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM1 pypeit run"
for n in 2024-03-15 2023-12-10 2023-12-09; do
  /opt/miniconda3/envs/pypeit/bin/python run_rm1.py --only $n
done
echo "STREAM2 COMPLETE"
