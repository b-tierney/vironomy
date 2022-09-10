#!/bin/bash

#SBATCH --job-name=vironomy_treebuild
#SBATCH --output=logs/vironomy_treebuild_%A_%a.log
#SBATCH --time=1-11:59:00 ### CHANGE THIS TO WHATEVER YOU WANT
#SBATCH --mem=32G ### CHANGE THIS TO WHATEVER YOU WANT
#SBATCH --cpus-per-task=12 ### CHANGE THIS TO WHATEVER YOU WANT
#SBATCH --partition=panda_physbio ### CHANGE THIS TO WHATEVER YOU WANT

samplefile=$1
treename=$(cat  $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")
treealgorithm=$(cat  $samplefile | cut -f2 |sed -n "${SLURM_ARRAY_TASK_ID}p")
alignmentloc=$(cat  $samplefile | cut -f3 |sed -n "${SLURM_ARRAY_TASK_ID}p")
treeoutloc=$(cat  $samplefile | cut -f4 |sed -n "${SLURM_ARRAY_TASK_ID}p")

threads=12 ### CHANGE THREADS HERE TOO

mkdir -p $treeoutloc

if [[ "$treealgorithm" == "fasttree" ]]; then
	export OMP_NUM_THREADS=$threads
	fasttree $alignmentloc > "$treeoutloc".fasttree.tree

if [[ "$treealgorithm" == "iqtree" ]]; then
	iqtree -s $alignmentloc --prefix $treeoutloc -m MFP --seqtype AA -T $threads









