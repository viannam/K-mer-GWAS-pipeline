#!/bin/bash
#SBATCH --job-name=kmers
#SBATCH --mail-user=viannam@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=120GB
#SBATCH --qos=mresende-b
#SBATCH --account=mresende
#SBATCH --time=96:00:00
#SBATCH --output=kmers_%A.out
#SBATCH --error=kmers_%A.err
##SBATCH --array=1-693  # Uncomment for per-sample parallel processing

# ------------------------------------------------------------------------------
# K-mers-based GWAS Pipeline
# Based on: https://github.com/voichek/kmersGWAS
# Author: Mariana Vianna (viannam@ufl.edu)
# ------------------------------------------------------------------------------

# === Load required modules ===
module load R/4.1
module load gemma/0.98.1
module load gcc/9.3.0
module load bamtools

# === Variables ===
WORKDIR="/blue/mresende/share/viannam/Kmers"
KMERS_BIN="${WORKDIR}/bin"
GENOTYPE_LIST="${WORKDIR}/Sample-genotypes.txt"
KMERS_LIST_PATH="${WORKDIR}/kmers_list_paths.txt"
KMERS_TABLE="kmers_table"
KMER_LENGTH=31
THREADS=10

# === Step 1: Generate sample-specific input files (if using array jobs) ===
# Un/Comment the following lines to enable per-sample execution using SLURM array
SAMPLE_ID=$(cut -f1 ${GENOTYPE_LIST} | sed -n ${SLURM_ARRAY_TASK_ID}p)
SAMPLE_PATH="${WORKDIR}/02.kmers/Genotypes/${SAMPLE_ID}"
mkdir -p ${SAMPLE_PATH}
 ls /orange/mresende/viannam/02.fasta/${SAMPLE_ID}/*${SAMPLE_ID}*R1*.fastq* \
    /orange/mresende/viannam/02.fasta/${SAMPLE_ID}/*${SAMPLE_ID}*R2*.fastq* \
    > ${SAMPLE_PATH}/input_file_${SAMPLE_ID}.txt

# === Step 2: Run KMC for canonical and non-canonical k-mer counting ===
KMC v3 command example (commented out, customize as needed)
${WORKDIR}/external_programs/kmc_v3 -t${THREADS} -k${KMER_LENGTH} -ci2 @input_file_${SAMPLE_ID}.txt output_kmc_canon_${SAMPLE_ID} ./ 
${WORKDIR}/external_programs/kmc_v3 -t${THREADS} -k${KMER_LENGTH} -ci0 -b @input_file_${SAMPLE_ID}.txt output_kmc_all_${SAMPLE_ID} ./ 

# === Step 3: Add strand information ===
# ${KMERS_BIN}/kmers_add_strand_information -c output_kmc_canon_${SAMPLE_ID} -n output_kmc_all_${SAMPLE_ID} -k ${KMER_LENGTH} -o kmers_with_strand_${SAMPLE_ID}

# === Step 4: List k-mers found in multiple samples ===
cd ${WORKDIR}
${KMERS_BIN}/list_kmers_found_in_multiple_samples \
    -l ${KMERS_LIST_PATH} \
    -k ${KMER_LENGTH} \
    --mac 5 \
    -p 0.2 \
    -o kmers_to_use

# === Step 5: Build k-mer presence/absence table ===
${KMERS_BIN}/build_kmers_table -l ${KMERS_LIST_PATH} -k ${KMER_LENGTH} -a kmers_to_use -o ${KMERS_TABLE}

# === Step 6: Calculate kinship matrix ===
${KMERS_BIN}/emma_kinship_kmers -t ${KMERS_TABLE} -k ${KMER_LENGTH} --maf 0.05 > ${KMERS_TABLE}.kinship

# === Step 7: Convert k-mer table to PLINK format ===
${KMERS_BIN}/kmers_table_to_bed -t ${KMERS_TABLE} -k ${KMER_LENGTH} -p germ.pheno --maf 0.01 --mac 1 -b 10000000 -o Output/germ/PLINK_kmers_germ

# === Step 8: Optional - extract subset of k-mers ===
${KMERS_BIN}/filter_kmers -t ${KMERS_TABLE} -k kmers_list.txt -o output.txt

# === Step 9: Run k-mers GWAS with permutation-based threshold ===
python2.7 kmers_gwas.py --pheno germ.pheno \
     --kmers_table ${KMERS_TABLE} \
     -l ${KMER_LENGTH} \
     -p ${THREADS} \
     -k 100000 \
     --dont_remove_intermediates \
     --maf 0.01 \
     --pattern_counter \
     --outdir GWAS_results/germ/
