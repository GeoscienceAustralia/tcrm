#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -N tcrm_ssp126_wind
#PBS -m ae
#PBS -M liang.hu@ga.gov.au
#PBS -lwalltime=48:00:00
#PBS -lmem=190GB,ncpus=48,jobfs=400GB
#PBS -joe
#PBS -lstorage=gdata/w85+scratch/w85+gdata/xp65
#PBS -l wd


module purge
module load pbs
module load dot
module load openmpi

module use /g/data/xp65/public/modules
module load conda/analysis3
source /scratch/w85/lh9573/EXTRA_PYTHON_LIBS/tcrm/bin/activate

mpirun -np 48 python3 tc_track_interpolation.py -c example/AU_ssp126.ini
mpirun -np 48 python3 tc_wind_generation.py -c example/AU_ssp126.ini