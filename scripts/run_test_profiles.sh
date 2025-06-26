#!/usr/bin/env bash

OUTDIR='test_profiles_out'

PROFILES=$(grep "includeConfig 'conf/test" nextflow.config | grep '_full' -v | awk '{print $1}')
while IFS= read -r PROF; do
    echo "=================================================================================================="
    echo "Running profile ${PROF} and saving output to ${OUTDIR}/${PROF}"
    echo "=================================================================================================="
    nextflow run main.nf -profile ${PROF},conda -resume --outdir ${OUTDIR}/${PROF}
done <<<"${PROFILES}"
