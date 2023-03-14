#!/bin/bash
#$ -cwd -V
#$ -o joblog.$JOB_ID
#$ -e joblog.$JOB_ID
#$ -l h_rt=96:00:00,h_data=4G,highp
#$ -pe shared 16
#$ -M mgandal@gmail.com
#$ -m bea


echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

. /u/local/Modules/default/init/modules.sh

module load anaconda3
module load samtools
source activate /u/project/gandalm/gandalm/mike_conda_base/


REF_GENOME="/u/project/gandalm/jops/isoseq/ref/GRCh37.primary_assembly.genome.fa"
REF_GENOME_NAME="GRCh37"
REF_ANNOT="/u/project/gandalm/jops/isoseq/ref/gencode.v33lift37.annotation.gtf"
REF_ANNOT_NAME="gencode.v33lift37"
BAM_File="/u/project/gandalm/jops/ucdavis_cpvz_isoseq/merged_bam/CP_VZ.hifi.bam"

isoquant.py \
--reference $REF_GENOME \
--genedb $REF_ANNOT \
--complete_genedb \
--bam $BAM_File \
--data_type pacbio_ccs \
--stranded forward \
--fl_data \
--sqanti_output \
--check_canonical \
--count_exons \
--read_group tag:RG \
--threads $(nproc --all) \
-o $SCRATCH/isoquant_hifi

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "