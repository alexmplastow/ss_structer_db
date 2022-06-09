#!/bin/bash

cd /N/slate/aplastow/PycharmProjects_slate/s_structure_db

source *venv/bin/activate

python ss_structure_db.py -i  hg38* -c -s 15 -ff 30 -cof calibration_output.csv -o calibration_structure.csv -rs

python folding_time.py -it *0221.fasta -ico calibration_output.csv
