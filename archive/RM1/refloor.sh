#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM1 pypeit run"
for n in 2024-03-15 2024-05-02; do
  /opt/miniconda3/envs/pypeit/bin/python run_rm1_throughput.py --only $n
done
echo "REFLOOR COMPLETE"
