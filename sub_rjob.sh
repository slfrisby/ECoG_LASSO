#!/bin/bash

# load R
module load R/4.3.1
# load compilers
module load gcc/9.5.0

dirp=/group/mlr-lab/Saskia/ECoG_LASSO
#q=01

# fits models for main analysis
#$dirp/scripts/./sub_fit_all_models.r --arg1=$dirp/derivatives/wavelet/sub-$q/

# runs eLife replication
#$dirp/scripts/./sub_elife_replic.r --arg1=$dirp/derivatives/wavelet/sub-$q/

# fits models at multiple preprocessing stages
$dirp/scripts/./sub_decode_during_preprocessing.r --arg1=$dirp/derivatives/wavelet/sub-$q/

# fits models to preprocessed electrooculogram
#$dirp/scripts/./sub_fit_all_eog_models.r --arg1=$dirp/derivatives/wavelet/sub-$q/



