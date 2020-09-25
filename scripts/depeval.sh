#!/bin/bash
#SBATCH --account=def-mmehride
#SBATCH --job-name="polyfusion_code"
#SBATCH --output="pf.%j.%N.out"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --constraint=broadwell
#SBATCH --mem=125G
#SBATCH --export=ALL
#SBATCH -t 10:00:00

binFile=$1
matrixPath=$2
opt=$3
par=$4

for f in ${matrixPath}/*.mtx
do
	$binFile $f $opt $par
done
