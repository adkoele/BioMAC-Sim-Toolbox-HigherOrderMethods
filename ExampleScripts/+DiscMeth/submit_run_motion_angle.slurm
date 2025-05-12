#!/bin/bash -l
#
# allocate 1 node (4 CPUs) for 24 hours
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:0:0
#SBATCH --cpus-per-task=4
#
# first non-empty non-comment line ends SBATCH  options

# load module (check Matlab versions with "module avail")
module load matlab/R2023a

echo ${SUBJECT}
echo ${MOTION}
echo ${TRIAL}
echo ${DISCMETH}
echo ${NNODES}

# build matlab code which should be called
WORKDIR="'~/DiscMeth_Workdirectory/'"
PATHIPOPT="'~/mexIPOPT/'"
FUNCTION="DiscMeth.run_motion_angle"
MATLABCODE="addpath(genpath(${PATHIPOPT})); addpath(genpath(${WORKDIR})); ${FUNCTION}(${WORKDIR}, '${SUBJECT}', '${MOTION}', '${TRIAL}'); quit;"
echo ${MATLABCODE}

# submit job
matlab -nodisplay -nosplash -r "${MATLABCODE}"

