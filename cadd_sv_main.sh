#! /bin/bash
# Run this script to do CADD-SV analysis
# STEP1. Make sure you have run customVCFanalysis_PacBio.nf
# STEP2. Make sure you have created BED files using tsv_to_BED_for_cadd_sv.R
# STEP3. Make sure you have copied the resulting files and config.yaml to /data/humangen_mouse/CADD-SV/
# STEP4. Close this file and execute it using ./cadd_sv_main.sh

sbatch cadd_sv_core.sh