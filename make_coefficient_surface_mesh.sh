#!/bin/bash

# this script is to be run after make_coefficient_volumes.m and before surface_plot.py.

# load workbench
module load workbench/1.5.0

# set directories
dirp=/group/mlr-lab/Saskia/ECoG_LASSO/
work=/group/mlr-lab/Saskia/ECoG_LASSO/work

# unzip template mesh from fsaverage
if [ ! -f $work/coefficients/pial_right.gii ]; then
cp $dirp/scripts/fsaverage/pial_left.gii.gz $work/coefficients/pial_left.gii.gz
cp $dirp/scripts/fsaverage/pial_right.gii.gz $work/coefficients/pial_right.gii.gz
gunzip $work/coefficients/pial_left.gii.gz
gunzip $work/coefficients/pial_right.gii.gz
fi

# run only if volumes containing coefficients exist
if [ -d $work/sub-"$ids"/coefficients/volume/ ]; then

# make output directory
rm -rf $work/sub-"$ids"/coefficients/surface/
mkdir $work/sub-"$ids"/coefficients/surface

# get a list of the volumes
cd $work/sub-"$ids"/coefficients/volume
dir > $work/sub-"$ids"/coefficients/volume.txt

# for each volume (see "done" statement to see how input is piped in from text)
while read volume; do

# get filename for reading and saving
filename=$(basename $volume .nii)

# if the coefficients are in the left hemisphere
if [[ $filename == L* ]]; then
# project to the left pial surface
wb_command -volume-to-surface-mapping $work/sub-"$ids"/coefficients/volume/"$filename".nii $work/coefficients/pial_left.gii $work/sub-"$ids"/coefficients/surface/"$filename".func.gii -trilinear
# smooth (6mm to match univariate - can be changed!). Treat zeros as missing data
wb_command -metric-smoothing $work/coefficients/pial_left.gii $work/sub-"$ids"/coefficients/surface/"$filename".func.gii 6 $work/sub-"$ids"/coefficients/surface/"$filename".func.gii -fix-zeros

# if they are in the right hemisphere
elif [[ $filename == R* ]];then 
# do the same thing in the right
wb_command -volume-to-surface-mapping $work/sub-"$ids"/coefficients/volume/"$filename".nii $work/coefficients/pial_right.gii $work/sub-"$ids"/coefficients/surface/"$filename".func.gii -trilinear
wb_command -metric-smoothing $work/coefficients/pial_right.gii $work/sub-"$ids"/coefficients/surface/"$filename".func.gii 6 $work/sub-"$ids"/coefficients/surface/"$filename".func.gii -fix-zeros

fi

done < $work/sub-"$ids"/coefficients/volume.txt

fi
 

