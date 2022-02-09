#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --partition=shared
#SBATCH --output=/gpfs/home/natem/logs/tialoc.txt

source /opt/miniconda3/bin/activate &&
snakemake -k -j 50 --use-conda