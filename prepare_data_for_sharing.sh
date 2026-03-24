#!/bin/bash

dirp=/group/mlr-lab/Saskia/ECoG_LASSO/derivatives/

# make a new directory
mkdir -p $dirp/OSF

# for every participant
for q in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 17 20 21 22; do

# make directory
mkdir $dirp/OSF/sub-"$q"/

# copy data into the OSF folder
# decibel-normalised power
cp $dirp/wavelet/sub-"$q"/dBpower.csv $dirp/OSF/sub-"$q"/dBpower.csv
cp $dirp/wavelet/sub-"$q"/dBpower.mat $dirp/OSF/sub-"$q"/dBpower.mat
# phase
cp $dirp/wavelet/sub-"$q"/phase.csv $dirp/OSF/sub-"$q"/phase.csv
cp $dirp/wavelet/sub-"$q"/phase.mat $dirp/OSF/sub-"$q"/phase.mat
# voltage
cp $dirp/wavelet/sub-"$q"/voltage.csv $dirp/OSF/sub-"$q"/voltage.csv
cp $dirp/wavelet/sub-"$q"/voltage.mat $dirp/OSF/sub-"$q"/voltage.mat

# zip directory
tar -zcvf $dirp/OSF/sub-"$q".tar.gz $dirp/OSF/sub-"$q"/

# delete unzipped
rm -r $dirp/OSF/sub-"$q"/

done

