#!/bin/bash
###Purpose: BWA and Samtools SCons wrapper run with GNU Parallel + Slurm
#For use on NERSC
#Original code written by William Arndt @ NERSC
#Modified by Ian Rambo

#SBATCH -t 00:30:00
#SBATCH -q debug
#SBATCH --nodes=2
#SBATCH --job-name=scons_map_scale_test
#SBATCH --mail-user=imrambo@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -C haswell
#SBATCH --output=slurm-%A_%a.out

module load python
module load parallel

source activate scons_map

srun -k -c 64 -n 16 -W 0 scons_map_parallel_payload.sh $1
