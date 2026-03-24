% decode based on just word frequency, to ensure that this cannot account
% for decoding performance in the main analysis. 

% set target
Y = [zeros(50,1);ones(50,1)];

% set features (taken from Morrison et al., 1997)
load('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/frequency.mat');
X = cell2mat(frequency(:,2));

% set other parameters that the script expects
acrossRun = 0;

% decode
output = get_log_hoerr(X,Y,'acrossRun',0);