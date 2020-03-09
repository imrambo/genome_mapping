#!/bin/bash
#For use on NERSC
#Original code written by William Arndt @ NERSC
#Modified by Ian Rambo

module load parallel

if [[ -z "${SLURM_NODEID}" ]]; then
echo "need \$SLURM_NODEID set"
exit
fi
if [[ -z "${SLURM_NNODES}" ]]; then
echo "need \$SLURM_NNODES set"
exit
fi

WORKDIR=/global/cscratch1/sd/imrambo/coassembly_MANERR_02-12-2020/genome_mapping
OUTDIR=/global/cscratch1/sd/imrambo/coassembly_MANERR_02-12-2020/manerr_coassembly_scons_map_debug
ASSEMBLY=/global/cscratch1/sd/imrambo/coassembly_MANERR_02-12-2020/data/final_assembly-min1000.fa

#Create the output directory if it does not exist
test -d $OUTDIR || mkdir -p $OUTDIR

cd $WORKDIR

SAM_THREAD_MEM=$(( ${SLURM_MEM_PER_NODE}/${SLURM_CPUS_ON_NODE} ))

PARALLEL_LOG=${HOME}/joblogs/scons_map_TEST_${SLURM_NODEID}.joblog
SCONS_LOG=${HOME}/joblogs/scons_map_TEST_${SLURM_NODEID}.log

sort $1 | awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" 'NR % NNODE == NODEID' | \
    parallel --jobs $SLURM_NNODES --joblog $PARALLEL_LOG \
    scons --dry-run --assembly=${ASSEMBLY} \
    --read_percent_id=97 --logfile=${SCONS_LOG} \
    --tmpdir=${SCRATCH}/tmp --nheader=20000 --samsort_mem=${SAM_THREAD_MEM} \
    --samsort_thread=${SLURM_CPUS_ON_NODE} --fastq_dir={//} --sampleids={/} \
    --outdir=${OUTDIR} --rm_local_build=0 --noIntraDepthVariance=1 --markdup=1 \
    --align_thread=${SLURM_CPUS_ON_NODE}
