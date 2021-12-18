#!/bin/bash
mu=3.18e-06
output_every_Xgen=6143
numgen_inN=12.10171
start_output=0.12286
cost=0.000814
for seed in 100
do
echo "
$seed
$mu
$cost
$output_every_Xgen
$numgen_inN
$start_output
" | ./Scripts/Viralevolution_100000 >./Output/Simulation/SimData/Data_T_0.318_cost_0.000814.txt
done
