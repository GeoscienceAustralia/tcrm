#!/bin/bash
#PBS -Pw85
#PBS -qnormal
#PBS -N calc_era5_mslp
#PBS -m ae
#PBS -M craig.arthur@ga.gov.au
#PBS -lwalltime=48:00:00
#PBS -lmem=64GB,ncpus=16,jobfs=8000MB
#PBS -W umask=0002
#PBS -joe
#PBS -lstorage=scratch/w85+gdata/ub4
#PBS -v STARTYEAR

module purge
module load pbs
module load dot

module load netcdf/4.6.3
module load cdo/1.9.8
module load nco/4.9.2
module load openmpi

# Suppresses an error related to HDF5 libraries:
export HDF5_DISABLE_VERSION_CHECK=2
ECHO=/bin/echo

cd ${PBS_JOBFS}

BASEPATH=/g/data/ub4/era5/netcdf/surface/msl
SCRATCH=/scratch/w85/

#STARTYEAR=1979
COUNTER=0
LOGFILE=$HOME/tcrm/calc_era5_mslp.$STARTYEAR.log
echo "Starting calculation of daily long term mean MSL" > $LOGFILE
for i in {0..40..1}; do
    for m in {1..12..1}; do
        YEAR=$(($STARTYEAR + $i))
        ENDDATE=`date -d "$YEAR/$m/1 +1 month -1 day" "+%Y%m%d"`
        STARTDATE=`date -d "$YEAR/$m/1" "+%Y%m%d"`
        DATESTR=$STARTDATE\_$ENDDATE
        INPUTFILE=$BASEPATH/$YEAR/msl_era5_global_$DATESTR.nc
        MONTHFMT=`date -d "$YEAR/$m/1" "+%Y%m"`
        OUTPUTFILE=${PBS_JOBFS}/msl_era5_global.$MONTHFMT.nc
        echo $INPUTFILE >> $LOGFILE
        echo $OUTPUTFILE >> $LOGFILE
        cdo daymean $INPUTFILE $OUTPUTFILE
        if [[ $? -ne 0 ]]; then
            echo "Looks like the command failed when processing $INPUTFILE" >> $LOGFILE
        else
            ((COUNTER+=1))
            echo "Processed $INPUTFILE ($COUNTER)" >> $LOGFILE
        fi
    done
    # concatenate individual monthly files into annual files
    cdo -b F64 mergetime ${PBS_JOBFS}/msl_era5_global.$YEAR[0-9][0-9].nc $SCRATCH/msl/msl_era5_global.$YEAR.nc
    # Remove individual monthly files:
    rm ${PBS_JOBFS}/msl_era5_global.$YEAR[0-9][0-9].nc
done

$ECHO "Finished processing $COUNTER files for daily MSLP data" >> $LOGFILE

#$ECHO "Now merging into a single daily long-term mean dataset" >> $LOGFILE

#cd $SCRATCH/msl

#cdo -P ${PBS_NCPUS} ydaymean -cat 'msl_era5_global.*.nc' msl_era5_global_dailyltm.nc

#if [[ $? -ne 0 ]]; then
#    $ECHO "cdo ydaymean command failed" >> $LOGFILE
#else
#    $ECHO "Finished processing $COUNTER files to create daily long-term mean MSL file" >> $LOGFILE
#fi
