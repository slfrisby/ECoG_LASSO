#!/bin/bash
echo "
++++++++++++++++++++++++" 

#p is set in sub_job.sh 
#p=$1
dirp=/group/mlr-lab/Saskia/ECoG_LASSO
work=/group/mlr-lab/Saskia/ECoG_LASSO/work

# runs preprocessing 
# matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');preprocess('"$p"');exit"

# runs wavelet transform
#matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');wavelet('"$p"');exit"

# extracts waveletted data at multiple preprocessing stages
matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');decode_during_preprocessing('"$p"');exit"

# extracts preprocessed electrooculogram (EOG)
#matlab_r2023b -nodisplay -nodesktop -r "addpath('$dirp/scripts/');preprocess_eyes('"$p"');exit"


