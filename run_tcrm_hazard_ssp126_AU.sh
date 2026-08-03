#!/bin/bash
#PBS -P w85
#PBS -q hugemem
#PBS -N tcrm_hazard_126
#PBS -m ae
#PBS -M liang.hu@ga.gov.au
#PBS -lwalltime=48:00:00
#PBS -lmem=1470GB,ncpus=48,jobfs=1400GB
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

mpirun -np 48 python3 tcrm.py -c example/AU_ssp126.ini