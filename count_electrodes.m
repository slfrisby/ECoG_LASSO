% by Saskia. Counts the number of electrodes implanted in each participant
% (for inclusion in Methods section). 

root = '/group/mlr-lab/Saskia/ECoG_LASSO';

subs = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','17','20','21','22'};
lat = {'L','L','L','L','L','R','L','R','L','L','L','R','L','L','L','L','L','L','L'};

% initialise counts
totalelecs = nan(1,length(subs));
vtelecs = nan(1,length(subs));
lvtelecs = nan(1,length(subs));
rvtelecs = nan(1,length(subs));

for q = 1:length(subs)

    p = subs{q};
    
    % read channels file as cell array
    channelsfile = table2cell(readtable([root,'/data/sub-',p,'/ieeg/sub-',p,'_task-naming_run-01_channels.tsv'],'FileType','Text','Delimiter','\t'));
    
    % initialise counters
    te = 0;
    vt = 0;
    lt = 0;
    rt = 0;

    % loop through rows of channels file
    for rowidx = 1:size(channelsfile,1)
        % if the electrode is labelled as 'ECOG', i.e. surface grid
        if strcmp(channelsfile(rowidx,2),'ECOG')
            % count it towards the total
            te = te + 1;
            % ventral temporal electrodes are additionally labelled as
            % 'good' or 'bad' (as opposed to others which are labelled
            % 'n/a'
            if strcmp(channelsfile(rowidx,6),'good') | strcmp(channelsfile(rowidx,6),'bad')
                % count it towards the ventral temporal total
                vt = vt + 1; 
            end
        end
        % store values
        totalelecs(q) = te;
        vtelecs(q) = vt;
        % if the patient is left hemisphere
        if lat{q} == 'L'
            % update the left hemisphere total
            lvtelecs(q) = vt;
        % else if they are right hemisphere
        elseif lat{q} == 'R'
            rvtelecs(q) = vt;
        end
    end
end

% calculate summary statistics
disp(['The mean of the total number of electrodes is ',num2str(mean(totalelecs)),'. The range of the total number of electrodes is ',num2str(min(totalelecs)),' - ',num2str(max(totalelecs)),'. The mean number of ventral temporal electrodes is ',num2str(mean(vtelecs)),'. The range of the number of ventral temporal electrodes is ',num2str(min(vtelecs)),' - ',num2str(max(vtelecs)),'.'])
disp(['For the left-hemisphere patients, the total number of electrodes is ',num2str(nanmean(lvtelecs)),'. The range is ',num2str(min(lvtelecs)),' - ',num2str(max(lvtelecs)),'. For the right-hemisphere patients, the total number of electrodes is ',num2str(nanmean(rvtelecs)),'. The range is ',num2str(min(rvtelecs)),' - ',num2str(max(rvtelecs)),'.'])