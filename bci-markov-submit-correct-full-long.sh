#!/bin/bash

#$ -N bci.markov

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.out
#$ -e /work/$USER/$JOB_NAME-$JOB_ID.err

#$ -l h_rt=400:00:00
#$ -l h_vmem=2G
##$ -l avx # avx is cpu feature that is only available in the newer generation of the hardware

#$ -pe smp 2

export MC_CORES=${NSLOTS:-1}

module load R

Rscript ~/bci/bci.markov.correct.full.r "$@" /work/purschke/bci-output/$JOB_NAME-$JOB_ID-$(basename $1 .Rdata)_$(basename $2 .Rdata)_${3}_${5}.Rdata
