#!/bin/bash

for i in {1..49}
do
sbatch run_model_Ne278_1070loci.sh
sleep 1

done
