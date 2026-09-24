#!/bin/zsh
cd "/Users/kaiwenzhang/PycharmProjects/KCWI/KCWI_RED_channel_characterization/RM1 pypeit run"
/opt/miniconda3/envs/pypeit/bin/python order_study_rm1.py --orders 3 5 7 9 15
echo "ORDERSTUDY COMPLETE"
