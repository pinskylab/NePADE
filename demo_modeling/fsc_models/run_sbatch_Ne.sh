#!/bin/bash

for i in {1..49}
do
sbatch run_model_Ne279_1068loci.sh
sleep 1

done
