#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --partition=largemem
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
##SBATCH --mail-user=%u@oist.jp
##SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=%j.out
##SBATCH --error=job_%j.err

. $HOME/.bashrc 
. ~/sasha_env/bin/activate

snakemake -j 999 -p --cluster-config cluster.json --cluster "sbatch  -p {cluster.partition} --ntasks {cluster.n} --mem-per-cpu {cluster.mem-per-cpu} --time {cluster.time} --nodes {cluster.nodes} --cpus-per-task {cluster.cpus-per-task} "
