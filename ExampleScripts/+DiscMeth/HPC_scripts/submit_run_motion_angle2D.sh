#!/bin/bash -l
#
# allocate 1 node (4 CPUs) for 24 hours
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=work
#SBATCH --constraint=icx
#SBATCH --output=slurm-%x-%j.out
# first non-empty non-comment line ends SBATCH  options

# load module (check Matlab versions with "module avail")
module load matlab/R2024b
module load gcc/11.2.0

echo ${SUBJECT}
echo ${MOTION}
echo ${TRIAL}
echo ${DISCMETH}
echo ${NNODES}

# build matlab code which should be called
WORKDIR="'~/DiscMeth2D_Workdirectory/'"
PATHIPOPT="'~/mexIPOPT/'"
FUNCTION="DiscMeth.run_motion_angle_2D"
MATLABCODE="addpath(genpath(${PATHIPOPT})); addpath(genpath(${WORKDIR})); ${FUNCTION}(${WORKDIR}, '${SUBJECT}', '${MOTION}', '${TRIAL}', '${DISCMETH}', '${NNODES}'); quit;"
echo ${MATLABCODE}

# submit job
matlab -nodisplay -nosplash -r "${MATLABCODE}"
