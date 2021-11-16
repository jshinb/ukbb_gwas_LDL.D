#!/bin/bash 
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=1:00:00
#SBATCH --job-name=LDL.D
#SBATCH --output=gwas_LDL.D_%j.txt
#SBATCH --mail-type=FAIL

# MAC filter removed
module load NiaEnv/2019b intel/2019u5 plink/.experimental-2.00alpha

resdir=$1
famfile=$2
phenofile=$3
covfile=$4
chr=$5

echo "SNPxSex analysis starts"
cd ${resdir}
$echo pwd
$echo ls -l $famfile
$echo ls -l $phenofile
$echo ls -l $covfile
echo $chr

mybfile=ukb43688_imp_chr${chr}_v3_filtered
plink2 --bgen /project/t/tpaus/tpaus/UKBB/datasets/gene_data/ukb_imp_chr${chr}_v3.bgen \
--sample /project/t/tpaus/tpaus/UKBB/datasets/gene_data/gene_data_backup/backup_gene_data_sample/ukb43688_imp_chr${chr}_v3_s487314.sample \
--maf 0.01 \
--geno 0.05 \
--hwe 1e-10 \
--remove /home/t/tpaus/mbernard/ukbiobank/subjects_to_remove_list.txt \
--keep $famfile \
--pheno $phenofile \
--covar $covfile \
--glm interaction --parameters 1-7 --tests 1, 3 --out SNPxSex_${chr}
