#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM2 pypeit run"
/opt/miniconda3/envs/pypeit/bin/python run_rm2.py --only 2024-06-11
/opt/miniconda3/envs/pypeit/bin/python run_rm2.py --only 2025-01-01
echo "SCI2 COMPLETE"
