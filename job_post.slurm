#!/bin/bash
#SBATCH --job-name=GEE_POST    # job's name
# --output = name of the output file  --error= name of error file (%j = jobID )
#SBATCH --output=POST_out  # nom du fichier de sortie
#SBATCH --error=POST_err   # nom du fichier d'erreur (fusion avec la sortie)
#SBATCH -N 1                        # number of nodes (or --nodes=1)
#SBATCH -n 2                       # number of cores
#SBATCH --mem=40G                  # --mem or--mem-per-cpu but not both
#SBATCH --time=23:59:00             # Wall Time 48h
#SBATCH --account=eolab       # MANDATORY : account  ( launch myaccounts to list your accounts)
#SBATCH --partition=cpu2022

module load python
module load geos
source /your/path/toward/venv/shoreliner
type -t module >/dev/null && module use /work/EOLAB/tools/FES/FES2022/modulefiles
module load ecola

python postprocess.py
