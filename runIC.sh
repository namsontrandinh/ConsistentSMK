#!/bin/bash
#SBATCH --job-name=IC-U
#SBATCH --partition=gpu
#SBATCH --ntasks=32
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err


./ic --graph email.bin --B_factor 0.01 --alg edl --csv email.csv
./ic --graph email.bin --B_factor 0.02 --alg edl --csv email.csv
./ic --graph email.bin --B_factor 0.03 --alg edl --csv email.csv
./ic --graph email.bin --B_factor 0.04 --alg edl --csv email.csv
./ic --graph email.bin --B_factor 0.05 --alg edl --csv email.csv


./ic --graph fb.bin --B_factor 0.01 --alg edl --csv fb.csv
./ic --graph fb.bin --B_factor 0.02 --alg edl --csv fb.csv
./ic --graph fb.bin --B_factor 0.03 --alg edl --csv fb.csv
./ic --graph fb.bin --B_factor 0.04 --alg edl --csv fb.csv
./ic --graph fb.bin --B_factor 0.05 --alg edl --csv fb.csv



