#!/bin/bash

# load ffmpeg
module load ffmpeg

dirp=/group/mlr-lab/Saskia/ECoG_LASSO/

mkdir $dirp/derivatives/figures/coefficients/videos/

# loop over analyses
for analysis in dB_power dB_theta_power dB_alpha_power dB_beta_power dB_gamma_power dB_high_gamma_power voltage; do

ffmpeg -r 10 -pattern_type glob -i "$dirp/derivatives/figures/coefficients/surface_plots/L_"$analysis"_timepoint_*_coefficients.jpeg" -filter:v minterpolate $dirp/derivatives/figures/coefficients/videos/L_"$analysis".mp4
ffmpeg -r 10 -pattern_type glob -i "$dirp/derivatives/figures/coefficients/surface_plots/R_"$analysis"_timepoint_*_coefficients.jpeg" -filter:v minterpolate $dirp/derivatives/figures/coefficients/videos/R_"$analysis".mp4

done
