#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt
#SBATCH --mem=26000
#SBATCH --time=00:50:00

let START=$SLURM_ARRAY_TASK_ID*$N_PER_JOB
let END=$START+${N_PER_JOB}-1

# if end>total tasks then end

for SUB_JOB_NUM in `seq $START $END`
do
	cd $SUB_JOB_NUM
	${ANTIGEN_ROOT}/run.py > "out.bashscreen"
	excessDiversityFile="out.tmrcaLimit"
	extinctFile="out.extinct"
	excessDiversity=${SWEEP_DIR}/results/$SUB_JOB_NUM/$excessDiversityFile
	extinct=${SWEEP_DIR}/results/$SUB_JOB_NUM/$extinctFile
	cd ${SWEEP_DIR}
	if [ ! -f "$excessDiversity" -a ! -f "$extinct" ]
	then
	MathKernel -noprompt -run < antigen-analysis-noFig.m ${SWEEP_DIR}/results/$SUB_JOB_NUM --display `hostname`:1 > "out.math"
	fi
	cd ${SWEEP_DIR}/results/
done

# Use this instead to restart runs that go extinct:
# ${ANTIGEN_ROOT}/run_py --with-restarts
