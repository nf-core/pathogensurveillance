#! /bin/env bash

# This script runs the pipeline for every profile and copies the main report input to into the direcory of test data (see variables below for specific paths used)

# Configuration variables
PROJECT_ROOT="pathogensurveillance"             # The name of the root directory for the project
PROFILES=(xanthomonas)                          # The names of the profiles to be used
RUNTIME_PROFILE="docker"                        # The name of the profile that supplies a way to run the pipeline (conda, docker, etc)
OUTPUT_DIR="test/output"                        # Where the output is to be saved. The output for each profile will be put in a directory named by the profile
REPORT_TEST_DIR="assets/main_report/_test_data" # Where this script will save a directory for each profile, with a subdirectory for each report group

# Check if we are in the right directory
if [[ "$(basename $PWD)" != "$PROJECT_ROOT" ]]; then
    printf "ERROR: You do not seem to be in the root of the pathogensurveillance directory. This script must be run from that location. If the name of the project directory has changed for some reason, then this script will have to be modified."
    exit 1
fi

# Run pathogensurveillance for each profile
for PROFILE in ${PROFILES[@]}; do
    printf "\\n\\n============ Runnig profile $PROFILE ============\\n\\n"
    # Run nextfow for this profiles
    PROFILE_OUT_DIR="$OUTPUT_DIR/$PROFILE"
    nextflow run main.nf -profile "$PROFILE,$RUNTIME_PROFILE" -resume --outdir "$PROFILE_OUT_DIR" --bakta_db assets/bakta_db/db-light
    REPORT_INPUT_DEST="${REPORT_TEST_DIR}/${PROFILE}"
    printf "\\n\\n============ Saving report input for $PROFILE to $REPORT_INPUT_DEST ============\\n\\n"
    mkdir -p $REPORT_INPUT_DEST
    rsync -av --delete --copy-links ${PROFILE_OUT_DIR}/prepare_report_input/* $REPORT_INPUT_DEST
done
