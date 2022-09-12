#!/bin/bash

#SBATCH --job-name=vironomy_treebuild
#SBATCH --output=logs/vironomy_treebuild_%A_%a.log
#SBATCH --time=1-11:59:00 
#SBATCH --mem=128G 
#SBATCH --cpus-per-task=12
#SBATCH --partition=panda_physbio 

samplefile=$1
treename=$(cat  $samplefile | cut -f1 |sed -n "${SLURM_ARRAY_TASK_ID}p")
treealgorithm=$(cat  $samplefile | cut -f2 |sed -n "${SLURM_ARRAY_TASK_ID}p")
alignmentloc=$(cat  $samplefile | cut -f3 |sed -n "${SLURM_ARRAY_TASK_ID}p")
treeoutloc=$(cat  $samplefile | cut -f4 |sed -n "${SLURM_ARRAY_TASK_ID}p")

threads=24 ### CHANGE THREADS HERE TOO

mkdir -p logs
mkdir -p $treeoutloc

if [[ "$treealgorithm" == "fasttree" ]]; then
	export OMP_NUM_THREADS=$threads
	FastTreeMP $alignmentloc > "$treeoutloc".fasttree.tree
fi

if [[ "$treealgorithm" == "iqtree" ]]; then
	iqtree -s $alignmentloc --prefix $treeoutloc -m MFP --seqtype AA -T $threads
fi
