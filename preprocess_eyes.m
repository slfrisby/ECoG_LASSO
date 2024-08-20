function preprocess_eyes(p)

% by Saskia. This function epochs electrooculogram recordings so that they
% can be decoded.
    addpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/')
    addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/eeglab2023.1'))
    addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/ECoG_Data_Prep-master'))
    root = '/group/mlr-lab/Saskia/ECoG_LASSO';
    cd(root);
 
    % make work subdirectory
    mkdir([root,'/work/sub-',p,'/eyes/']);

    % tell us where it is up to
    disp(['Working on subject ',p,' ...']);

    %% Part 1 - epoch
    
    % copy metadata input into work subdirectory
    copyfile('/group/mlr-lab/Saskia/ECoG_LASSO/work/metadata_input',[root,'/work/sub-',p,'/eyes/metadata_input'],'f');

    % make metadata to include just EOG channels

    % locate channels.tsv
    cd([root,'/data/sub-',p,'/ieeg/']);
    % Selects all channels.tsv files in case the semantic judgement .tsvs
    % are needed later. For now, we use the naming files only.
    chanfiles = dir(['/group/mlr-lab/Saskia/ECoG_LASSO/data/sub-',p,'/ieeg/sub-',p,'_task-*_channels.tsv']);
    % read naming channels.tsv file
    channels = table2cell(readtable(chanfiles(1).name,'FileType','text','Delimiter','\t'));
    % initialise
    chanidx = [];
    % select channels
    for i = 1:size(channels,1)
        % if the channel is in the vATL, mark it for keeping
        if strcmp(channels{i,2},'EOG') 
            chanidx(i) = 1;
        % else mark it for chucking
        else
           chanidx(i) = 0;
        end
    end
    channels = channels(logical(chanidx'),1);
    % participant 22 has 6 EOG channels and 4 of these are redundant (i.e.
    % channel 1 is a summary of channels 3 and 4 and channel 2 is a summary
    % of channels 5 and 6). So just keep the first 2
    if strcmp(p,'22') 
        channels = channels(1:2);
    end
    % make metadata. EOG does not have coordinates, so just set these to 0
    xyz = repmat([0 0 0],size(channels,1),1);
    subject(1:size(xyz,1),1) = {['sub-',p]};
    ELECTRODE(str2num(p)) = {channels};
    % save variables 
    save([root,'/work/sub-',p,'/eyes/metadata_input/coords/electrodes.mat'],'ELECTRODE','subject','xyz');
    spreadsheet = table;
    spreadsheet.subject = subject;
    spreadsheet.electrode = vertcat(ELECTRODE{:});
    spreadsheet.x = xyz(:,1);
    spreadsheet.y = xyz(:,2);
    spreadsheet.z = xyz(:,3);
    writetable(spreadsheet,[root,'/work/sub-',p,'/eyes/metadata_input/coords/MNI_basal_electrodes_sub-01-22_w_label.csv'])

    % copy raw data into work subdirectory. The filename is 'filt_data.mat'
    % even though the data are not filtered! This is because, in the main
    % workflow, setup_data is run after filtering, so this is the filename
    % that setup_data is set up to expect.
    copyfile([root,'/data/sub-',p,'/ieeg/sub-',p,'_task-naming_run-01_ieeg.mat'],[root,'/work/sub-',p,'/eyes/filt_data.mat'])
    % epoch data using setup_data
    setup_data( ...
    'WindowStart', -1000, ... 
    'WindowSize', 4000, ...
    'BaselineWindow', 200, ...
    'subjects', str2num(p), ...
    'boxcar', 1, ...
    'slope_interval', 0, ...
    'average', 0, ...
    'datacode', '', ...
    'dataroot', [root,'/work/sub-',p,'/eyes/'], ...
    'metaroot', [root,'/work/sub-',p,'/eyes/metadata_input/'], ...
    'datarootout', ['/group/mlr-lab/Saskia/ECoG_LASSO/work/sub-',p,'/eyes/'], ...
    'cvpath', [root,'/work/metadata_input/cv/cvpartition_10fold_DilkinaSplit.mat'], ...
    'overwrite', 1, ...
    'WriteIndividualMetadata', 1);

    % reorganise data ready for wavelet decomposition 
    a = load([root,'/work/sub-',p,'/eyes/full/BoxCar/001/WindowStart/-1000/WindowSize/4000/s',p,'_.mat']);
    disp(['Organising data for subject ',p,' ...']);
    % find the number of electrodes (the second dimension is spatiotemporal
    % features - value for each timepoint for electrode 1, followed by 
    % data for each timepoint for electrode 2, etc.)
    nelecs = size(a.X,2)/4000;
    % reshape into a matrix of trials x time x electrodes
    X = reshape(a.X,400,4000,nelecs);
    % permute into the format electrodes x time x trials (data are not in perfect Cox
    % format at this point)
    X = permute(X,[3 2 1]);
    % save and clear
    save([root,'/work/sub-',p,'/eyes/EOG.mat'],'X','-v7.3');
    disp('Done organising!')
    
    %% Part 2 - wavelet transform

    % We want to test whether EOG is decodable:
    % - in its raw form
    % - in a warped/smeared form following wavelet decomposition
    % Therefore, we run the data through the same wavelet pipeline as we
    % use for the ECoG electrodes.

    %% Part 2A - power

    power = zeros(size(X,3), size(X,1), 166, 60);
    baselinepower = zeros(size(X,3), size(X,1), 21, 60);

    % we want these results to be as comparable to the main analysis as
    % possible. Load main metadata and set bad trials to nans.
    load([root,'/work/sub-',p,'/metadata.mat']);
    X(:,:,~metadata.filters(10).filter)=NaN;

    % wavelet
    for elec = 1:size(X,1)
        % wavelet main data
        [output, freqs, times] = timefreq(squeeze(X(elec,:,:)), 1000, 'cycles', [5 15], 'freqs', [4 200], 'nfreqs', 60, 'freqscale', 'log', 'timesout',[1001:10:2651]);
        % wavelet baseline
        [baselineoutput, freqs, times] = timefreq(squeeze(X(elec,:,:)), 1000, 'cycles', [5 15], 'freqs', [4 200], 'nfreqs', 60, 'freqscale', 'log', 'timesout',[701:10:901]);
        % calculate power
        pow = output.*conj(output);
        baselinepow = baselineoutput.*conj(baselineoutput);
        % store
        power(:,elec,:,:) = permute(pow,[3,2,1]);
        baselinepower(:,elec,:,:) = permute(baselinepow,[3,2,1]);
    end

    % non-normalised power first. Calculate mean over items, excluding
    % nans if present
    tmp = cat(5,power(1:100,:,:,:),power(101:200,:,:,:),power(201:300,:,:,:),power(301:400,:,:,:));
    power = nanmean(tmp,5);

    % if there are missing trials, NaNs will remain. Replace these with the
    % median power across trials (for that electrode, band and frequency).
    tmp = find(isnan(power));
    [trial,elec,window,freq] = ind2sub(size(power),tmp);
    for i = 1:length(trial)
        power(trial(i),elec(i),window(i),freq(i)) = nanmedian(power(:,elec(i),window(i),freq(i)));
    end
   
    % save the data as a matlab variable
    save([root,'/derivatives/wavelet/sub-',p,'/eyes_power.mat'],'power');

    % reformat into trials x spatiofrequency features
    tmp = permute(power,[1,3,2,4]);
    powersff = reshape(tmp,size(tmp,1),[]);
   
    % write .csv
    writematrix(powersff,[root,'/derivatives/wavelet/sub-',p,'/eyes_power.csv']);

    % then power with decibel normalisation. Calculate condition-average
    % baseline
    baselinepower = squeeze(nanmean(baselinepower,[1,3]));
    % make it the same size as the data
    baselinepower = permute(repmat(baselinepower,1,1,100,166),[3,1,4,2]);
    % decibel-convert
    dBpower = 10*log10(power./baselinepower);
    % save
    save([root,'/derivatives/wavelet/sub-',p,'/eyes_dBpower.mat'],'dBpower');

    % reformat into trials x spatiofrequency features
    tmp = permute(dBpower,[1,3,2,4]);
    dBpowersff = reshape(tmp,size(tmp,1),[]);
    % write .csv
    writematrix(dBpowersff,[root,'/derivatives/wavelet/sub-',p,'/eyes_dBpower.csv']);


     %% Part 2B - Phase

    % initialise output
    phase = zeros(size(X,3)/4, size(X,1), 166, 60);

    % average data
    tmp = cat(4,X(:,:,1:100),X(:,:,101:200),X(:,:,201:300),X(:,:,301:400));
    avX = nanmean(tmp,4);

   % wavelet
    for elec = 1:size(X,1)
        [output, freqs, times] = timefreq(squeeze(avX(elec,:,:)), 1000, 'cycles', [5 15], 'freqs', [4 200], 'nfreqs', 60, 'freqscale', 'log', 'timesout',[1001:10:2651]); 
        % calculate phase
        pha = angle(output);
        % store
        phase(:,elec,:,:) = permute(pha,[3,2,1]);
    end

   % median-impute remaining NaNs
   tmp = find(isnan(phase));
   [trial,elec,window,freq] = ind2sub(size(phase),tmp);
   for i = 1:length(trial)
        phase(trial(i),elec(i),window(i),freq(i)) = nanmedian(phase(:,elec(i),window(i),freq(i)));
   end

   % save
   save([root,'/derivatives/wavelet/sub-',p,'/eyes_phase.mat'],'phase');

   % reformat into trials x spatiofrequency features
   tmp = permute(phase,[1,3,2,4]);
   phasesff = reshape(tmp,size(tmp,1),[]);
   % write .csv
   writematrix(phasesff,[root,'/derivatives/wavelet/sub-',p,'/eyes_phase.csv']);

   %% Part 2C - voltage

   % average data
   tmp = cat(4,X(:,:,1:100),X(:,:,101:200),X(:,:,201:300),X(:,:,301:400));
   voltage = permute(nanmean(tmp,4),[3,1,2]);

   % median-impute remaining NaNs
   tmp = find(isnan(voltage));
   [trial,elec,timepoint] = ind2sub(size(voltage),tmp);
   for i = 1:length(trial)
        voltage(trial(i),elec(i),timepoint(i)) = nanmedian(voltage(:,elec(i),timepoint(i)));
   end

   % save
   save([root,'/derivatives/wavelet/sub-',p,'/eyes_voltage.mat'],'voltage');

   % reformat into trials x spatiofrequency features
   tmp = permute(voltage,[1,3,2]);
   voltagesff = reshape(tmp,size(tmp,1),[]);
   % write .csv
   writematrix(voltagesff,[root,'/derivatives/wavelet/sub-',p,'/eyes_voltage.csv']);
    
end








