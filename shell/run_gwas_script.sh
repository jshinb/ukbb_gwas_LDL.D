#####
#https://rgcgithub.github.io/regenie/recommendations/
#http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas
sfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_LDL.D_gwas/gwas_script.sh

#----------------------- VARIABLE ORDER --------------------------------------#
#resdir=$1
#famfile=$2
#phenofile=$3
#covfile=$4
#chr=$5
#-----------------------------------------------------------------------------#

#INPUT (4):
dirname='UKB_avg.LDL.D_20211102'
famfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_LDL.D_gwas/data/FAM.txt
phenofile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_LDL.D_gwas/data/PHENO.txt
covfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_LDL.D_gwas/data/COV.txt

resdir=${SCRATCH}/PLINK/$dirname
mkdir $resdir
cd $resdir

for chr in {2..22}
do
sbatch $sfile $resdir $famfile $phenofile $covfile $chr
done

# runsfile=/gpfs/fs1/home/t/tpaus/jshinb/ukbb_LDL.D_gwas/run_gwas_script.sh;sh $runsfile
