#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM1 pypeit run"
for n in 2024-05-02 2023-12-11; do
  /opt/miniconda3/envs/pypeit/bin/python run_rm1.py --only $n
done
echo "STREAM3 COMPLETE"
