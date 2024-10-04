#!/bin/bash
#SBATCH --job-name=launcher    # job's name
# --output = name of the output file  --error= name of error file (%j = jobID )
#SBATCH --output=launcher_out  # nom du fichier de sortie
#SBATCH --error=launcher_err   # nom du fichier d'erreur (fusion avec la sortie)
#SBATCH -N 1                        # number of nodes (or --nodes=1)
#SBATCH -n 1                       # number of cores
#SBATCH --mem=1G                  # --mem or--mem-per-cpu but not both
#SBATCH --time=00:09:00             # Wall Time 48h
#SBATCH --account=eolab       # MANDATORY : account  ( launch myaccounts to list your accounts)
#SBATCH --partition=cpu2022

export http_proxy="http://user:password@proxy-surf.loc.cnes.fr:8050"
export https_proxy="http://user:password@proxy-surf.loc.cnes.fr:8050"
export ftp_proxy="${http_proxy}"
export no_proxy=cnes.fr,sis.cnes.fr,gitlab.cnes.fr

module load python
module load geos
source /your/path/toward/venv/shoreliner

/your/path/toward/shoreliner/folder
python launcher.py

