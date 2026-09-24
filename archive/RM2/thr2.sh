#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM2 pypeit run"
for n in 2024-06-11 2025-01-01 2025-01-02; do
  /opt/miniconda3/envs/pypeit/bin/python run_rm2_throughput.py --only $n
done
echo "THR2 COMPLETE"
