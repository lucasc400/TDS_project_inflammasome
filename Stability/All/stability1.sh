#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -N stability1

cd /rds/general/user/sc5119/home/TDS_project/Scripts/Lucas/Stability/All
module load anaconda3/personal

Rscript Stability1.R


