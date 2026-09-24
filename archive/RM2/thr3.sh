#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM2 pypeit run"
for n in 2023-09-23 2024-06-10 2024-12-24 2024-05-09; do
  /opt/miniconda3/envs/pypeit/bin/python run_rm2_throughput.py --only $n
done
echo "THR3 COMPLETE"
