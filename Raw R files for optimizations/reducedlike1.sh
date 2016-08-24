#PBS -l pmem=350mb 
#PBS -l nodes=1:ppn=1 
#PBS -W group_list=zferguson
#PBS -m abe 
#PBS -l walltime=64:00:00 
 
cd /scratch/lfs/rzenil/bichrom/reducedlikeED 
module load R 
R CMD BATCH explorelikereduced1.R