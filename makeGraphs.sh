#!/bin/bash

#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y
module load conda_R

Rscript "makeGraphs.R"

