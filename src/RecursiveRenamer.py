import os
from pathlib import Path

"""Based on set conditions renames all files in a directory."""

base_dir = "/Users/dennis/Dev/STARE_GAZE/Zenodo/INVOKE/INVOKE_HockerCicero_Prom400_AllTSSAll_10Fold/"
all_files = [str(x) for x in list(Path(base_dir).rglob("*.*")) if 'DS_Store' not in str(x)]

for file in all_files:
    os.rename(file, '/'.join(file.split('/')[:-1]) + '/' + file.split('/')[-1].replace('HockerCicero_Prom400_AllTSS', 'Hocker_Cicero'))
