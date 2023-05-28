#!/bin/sh
#SBATCH -J "demo check"
#SBATCH -N 4
#SBATCH -n 32
#SBATCH -A hpc-prf-lgrsfpga
#SBATCH -p normal
#SBATCH --mem-per-cpu 2G
#SBATCH -t 6:00:00
#SBATCH -o %x-%j.log
#SBATCH -e %x-%j.log



cd /scratch/hpc-prf-lgrsfpga/yashwanth/codethesis/switch_Floorplanner

make GUROBI_USE=1

FILES = 'MBLA_Data/*'

for file in $FILES
do
        ./final_fpgafloorplanner "$file" 
        
done


exit 0