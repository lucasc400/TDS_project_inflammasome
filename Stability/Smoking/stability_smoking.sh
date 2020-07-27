#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N stability_smoking

cd /rds/general/user/sc5119/home/TDS_project/Scripts/Lucas/Stability/Smoking
module load anaconda3/personal

Rscript Stability_smoking.R


