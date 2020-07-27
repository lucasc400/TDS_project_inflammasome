#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N prot_stability

cd /rds/general/user/sc5119/home/TDS_project/Scripts/Lucas/sPLS
module load anaconda3/personal

Rscript pls_stability_protein.R


