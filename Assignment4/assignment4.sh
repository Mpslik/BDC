#!/bin/bash
#SBATCH --job-name=BDC_Assignment4
#SBATCH --nodes=1
#SBATCH --partition=assemblix
#SBATCH --nodelist=assemblix2019
#SBATCH --ntasks=5
#SBATCH --time=00:10:00
#SBATCH --mem=4G


# Example arguments to assignment4.py:
# -n 4 : 4 chunks per file (like assignment1).
# --files : your FastQ files
# --output : if you want CSV output instead of printing to stdout

mpirun -n 5 python3 assignment4.py \
    -n 4 \
    --files /commons/Themas/Thema12/HPC/rnaseq.fastq

echo "Done with assignment4 job!"