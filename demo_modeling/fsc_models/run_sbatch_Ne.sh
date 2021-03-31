#!/bin/bash

for i in {1..50}
do
sbatch run_model6_singlethread.sh
sleep 1

done
