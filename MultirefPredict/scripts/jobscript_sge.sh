#!/usr/bin/env bash
#$ -S /bin/bash
#$ -N diagnostics
#$ -R y
#$ -cwd
#$ -l h_rt=64:00:00
#$ -l h_rss=8G
#$ -q cpus
#$ -l cpus=1
#$ -pe smp 1
# -fin run_diagnostics_ebased.py
#$ -cwd

module load anaconda
source activate qcengine

python run_diagnostics_ebased.py > $SGE_O_WORKDIR/diag.out

mv * $SGE_O_WORKDIR
