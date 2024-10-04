#!/bin/bash
#SBATCH --job-name=GEE    # job's name
# --output = name of the output file  --error= name of error file (%j = jobID )
#SBATCH --output=GEE_out  # nom du fichier de sortie
#SBATCH --error=GEE_err   # nom du fichier d'erreur (fusion avec la sortie)
#SBATCH -N 1                        # number of nodes (or --nodes=1)
#SBATCH -n 1                       # number of cores
#SBATCH --mem=5G                  # --mem or--mem-per-cpu but not both
#SBATCH --time=71:59:59             # Wall Time 48h
#SBATCH --account=eolab       # MANDATORY : account  ( launch myaccounts to list your accounts)
#SBATCH --partition=cpu2022

python download.py 
python main_process.py 

