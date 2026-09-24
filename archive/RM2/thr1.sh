#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM2 pypeit run"
for n in 2024-01-04 2024-04-30; do
  /opt/miniconda3/envs/pypeit/bin/python run_rm2_throughput.py --only $n
done
echo "THR1 COMPLETE"
