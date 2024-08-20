% by Saskia. Sets up metadata input folder, including coordinates, needed
% for calling setup_data() inside preprocess().

subs = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','17','20','21','22'};

% copy folder from ECoG_Data_Prep master. The contents of this folder will
% be updated to suit this dataset. N.B. Running this script on the MRC CBU
% compute cluster will call the CBU wrapper function copyfile, which calls
% the native MATLAB copyfile without throwing errors.
copyfile('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/ECoG_Data_Prep-master/data','/group/mlr-lab/Saskia/ECoG_LASSO/work/metadata_input','f');

% setup outputs
ELECTRODE = {};
subject = {};
xyz = [];

% setup variables
allchannels = {};
counter = 0;

for q = 1:length(subs)

    p = subs{q};

    cd(['/group/mlr-lab/Saskia/ECoG_LASSO/data/sub-',p,'/ieeg/']);
    
    % Selects all channels.tsv files in case the semantic judgement .tsvs
    % are needed later. For now, we use the naming files only.
    chanfiles = dir(['/group/mlr-lab/Saskia/ECoG_LASSO/data/sub-',p,'/ieeg/sub-',p,'_task-*_channels.tsv']);
   
    % Step 1 - get names of ventral temporal electrodes not including those
    % in the seizure onset zone.

    % read naming channels.tsv file
    channels = table2cell(readtable(chanfiles(1).name,'FileType','text','Delimiter','\t'));
    chanidx = [];

    for i = 1:size(channels,1)
        % if the channel is marked as 'good', mark it for keeping
        if strcmp(channels{i,6},'good')
            chanidx(i) = 1;
        % else mark it for chucking (it is either 'bad', i.e. below the
        % seizure onset zone or containing a very large number of
        % artifacts, or 'n/a', i.e. outside the region of interest)
        else
           chanidx(i) = 0;
        end
    end
    % keep only channels marked for keeping
    channels = channels(logical(chanidx'),1);

    % Step 2 - get coordinates of those electrodes

    coords = table2cell(readtable(['sub-',p,'_electrodes.tsv'],'FileType','text','Delimiter','\t'));

    % save the first column as a separate variable
    coords1 = coords(:,1);

    % then remove everything but the separate coordinate values and convert
    % contents of the cell array to double if necessary (coords remains as
    % a cell array).
    if ischar(coords{1,2}) % this asks what sort of object the first x coordinate is
        coords = num2cell(str2double(coords(:,2:4)));
    else
        coords = coords(:,2:4);
    end

    % initialise output
    coords2keep = [];
    
    for i = 1:size(channels,1)
        % find index of electrode in coords
        rowidx = find(ismember(coords1,channels(i)));
        % find coordinates 
        coords2keep(i,:) = [coords{rowidx,1},coords{rowidx,2},coords{rowidx,3}];
    end

    % Step 3 - construct metadata input to setup_data

    skip = size(channels,1);
    xyz(counter+1:counter+skip,:) = coords2keep;
    subject(counter+1:counter+skip,1) = {['sub-',p]};
    ELECTRODE(q,1) = {channels};

    counter = counter + size(channels,1);
    
end

% save variables 
save('/group/mlr-lab/Saskia/ECoG_LASSO/work/metadata_input/coords/electrodes.mat','ELECTRODE','subject','xyz');

% construct excel spreadsheet
spreadsheet = table;
spreadsheet.subject = subject;
spreadsheet.electrode = vertcat(ELECTRODE{:});
spreadsheet.x = xyz(:,1);
spreadsheet.y = xyz(:,2);
spreadsheet.z = xyz(:,3);
writetable(spreadsheet,'/group/mlr-lab/Saskia/ECoG_LASSO/work/metadata_input/coords/MNI_basal_electrodes_sub-01-22_w_label.csv')

