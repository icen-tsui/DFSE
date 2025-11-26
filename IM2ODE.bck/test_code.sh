#!/bin/bash

make
cp de.x ../../test
cd ../../test
sbatch de.sh
tail -f run_de.out

