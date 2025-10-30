#!/bin/bash
#PBS -P w85
#PBS -q normal
#PBS -N tcrm_WP_ssp126_18_reform
#PBS -m ae
#PBS -M liang.hu@ga.gov.au
#PBS -lwalltime=48:00:00
#PBS -lmem=190GB,ncpus=48,jobfs=400GB
#PBS -joe
#PBS -lstorage=gdata/w85+scratch/w85+gdata/hh5
#PBS -l wd


module purge
module load pbs
module load dot
module load openmpi

module use /g/data/hh5/public/modules
module load conda/analysis3
source /scratch/w85/lh9573/EXTRA_PYTHON_LIBS/tcrm/bin/activate

# mpirun -np 48 python3 tcrm.py -c example/WP_ssp126_18.ini
mpirun -np 24 python3 tc_track_interpolation.py -c example/WP_ssp126_10000_reform.ini