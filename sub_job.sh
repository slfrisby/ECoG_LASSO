#!/bin/bash

dirp=/group/mlr-lab/Saskia/ECoG_LASSO
if [ ! -d $dirp/work/logs/ ]; then
mkdir -p $dirp/work/logs/
fi

# test
#for q in 01; do
# run all participants
for q in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 17 20 21 22; do
# run original elife participants
#for q in 01 02 03 04 05 07 09 10; do
# run participants with EOG
#for q in 01 02 03 06 22; do

echo "$q"

# runs hands-off section of preprocessing 
# sbatch -o $dirp/work/logs/"$q"preprocess.out -c 16 --job-name=ECoG_"$q" --export=p=${q} $dirp/scripts/sub_matlabjob.sh

# runs wavelet transform
#sbatch -o $dirp/work/logs/"$q"wavelet.out -c 16 --job-name=ECoG_"$q" --export=p=${q} $dirp/scripts/sub_matlabjob.sh

# fits models (standard queue or lopri)
#sbatch -o $dirp/work/logs/"$q"modelfit.out -c 24 -t 4-0:00 --job-name=ECoG_"$q" --export=q=${q} $dirp/scripts/sub_rjob.sh 
#sbatch -o $dirp/work/logs/"$q"modelfit.out -c 24 -q lopri -t 4-0:00 --job-name=ECoG_"$q" --export=q=${q} $dirp/scripts/sub_rjob.sh 

# extracts waveletted data at multiple preprocessing stages
# sbatch -o $dirp/work/logs/"$q"DDP.out -c 16 --job-name=ECoG_"$q" --export=p=${q} $dirp/scripts/sub_matlabjob.sh

# fits models at multiple preprocessing stages
#sbatch -o $dirp/work/logs/"$q"DDPmodels.out -c 24 -t 7-0:00 --job-name=ECoG_"$q" --export=q=${q} $dirp/scripts/sub_rjob.sh 
#sbatch -o $dirp/work/logs/"$q"DDPmodels.out -c 24 -q lopri -t 7-0:00 --job-name=ECoG_"$q" --export=q=${q} $dirp/scripts/sub_rjob.sh 

# extracts preprocessed electrooculogram (EOG)
#sbatch -o $dirp/work/logs/"$q"eyes.out -c 16 --job-name=00_"$q" --export=p=${q} $dirp/scripts/sub_matlabjob.sh

# fits models to preprocessed electrooculogram
#sbatch -o $dirp/work/logs/"$q"eyemodels.out -c 24 -t 4-0:00 --job-name=00_"$q" --export=q=${q} $dirp/scripts/sub_rjob.sh 
#sbatch -o $dirp/work/logs/"$q"eyemodels.out -c 24 -q lopri -t 4-0:00 --job-name=00_"$q" --export=q=${q} $dirp/scripts/sub_rjob.sh 

# makes coefficient surface meshes for plotting
#sbatch -o $dirp/work/logs/"$q"surface.out -c 16 --job-name=COEF_"$q" --export=ids=${s} $dirp/scripts/make_coefficient_surface_mesh.sh
sbatch -o $dirp/work/logs/"$q"surface.out -c 16 -q lopri --job-name=COEF_"$q" --export=ids=${q} $dirp/scripts/make_coefficient_surface_mesh.sh

done

#01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 17 20 21 22

	

