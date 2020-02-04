#!/bin/bash

#SBATCH -n 1
#SBATCH -p sched_mit_chisholm
#SBATCH --time=48:00:00

source activate snakemake;

snakemake --use-conda --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.cores} --mem={cluster.mem} --time={cluster.time} --exclude=node420,node421,node422" --jobname {rulename}.{jobid} --jobs 250 --keep-going
