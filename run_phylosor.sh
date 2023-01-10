#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --partition=shared
#SBATCH --output=/gpfs/home/natem/logs/phylosor.txt

module load python/3.8.3 &&
snakemake -k -j 20 --use-conda --cluster-config config/cluster.json --cluster "sbatch --time={cluster.walltime} --mem={cluster.mem} -c {cluster.n} --partition={cluster.queue} --output={cluster.logfile}"
