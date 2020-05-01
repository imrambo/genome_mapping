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

#Directory containing SCons wrapper; script needs to run from here
workdir=/global/cscratch1/sd/imrambo/coassembly_MANERR_02-12-2020/genome_mapping
#Output directory for SCons wrapper
outdir=/global/cscratch1/sd/imrambo/coassembly_MANERR_02-12-2020/manerr_coassembly_scons_map_debug_run00
#FASTA assembly to map to
assembly=/global/cscratch1/sd/imrambo/coassembly_MANERR_02-12-2020/data/final_assembly-min1000_toy1000.fna

#Create the output directory if it does not exist
test -d $outdir || mkdir -p $outdir

#Memory per samtools sort thread
sam_thread_mem=$(( ${SLURM_MEM_PER_NODE}/${SLURM_CPUS_ON_NODE} ))

#Logfiles
parallel_log=${HOME}/joblogs/scons_map_TEST_${SLURM_NODEID}.joblog
scons_log=${HOME}/joblogs/scons_map_TEST_${SLURM_NODEID}.log

cd $workdir && \
    sort $1 | awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" 'NR % NNODE == NODEID' | \
    parallel --jobs $SLURM_NNODES --joblog $parallel_log \
    scons --dry-run --assembly=${assembly} \
    --read_percent_id=97 --logfile=${scons_log} \
    --tmpdir=${SCRATCH}/tmp --nheader=20000 --samsort_mem=${sam_thread_mem} \
    --samsort_thread=${SLURM_CPUS_ON_NODE} --fastq_dir={//} --sampleids={/} \
    --outdir=${outdir} --rm_local_build=0 --noIntraDepthVariance=1 --markdup=1 \
    --align_thread=${SLURM_CPUS_ON_NODE}
