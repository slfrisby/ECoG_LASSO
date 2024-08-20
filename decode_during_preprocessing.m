function decode_during_preprocessing(p)

% by Saskia. This function enables testing of the influence of various
% preprocessing steps by producing data ready for decoding at intermediate
% stages of the preprocessing pipeline:
% - raw
% - after filtering 
% - after exlucsion of bad channels
% - after common average referencing
% - after rejection of bad trials
% This function synthesises preprocess, inspect_data, wavelet and
% make_metadata_input. The other functions are more thoroughlmkdiry commented so
% can explain the rationale for processing steps within this function.

    addpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/')
    addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/eeglab2023.1'))
    addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/ECoG_Data_Prep-master'))

    root = '/group/mlr-lab/Saskia/ECoG_LASSO';
    cd(root);

    % make work subdirectory
    mkdir([root,'/work/sub-',p,'/decode_during_preprocessing/']);

    % tell us where it is up to
    disp(['Working on subject ',p,' ...']);

    %% Part 1 - extraction and organisation of preprocessed data
    
    %% Part 1A - raw data

    % copy metadata input into work subdirectory
    copyfile('/group/mlr-lab/Saskia/ECoG_LASSO/work/metadata_input',[root,'/work/sub-',p,'/decode_during_preprocessing/metadata_input'],'f');

    % make metadata to include all ventral temporal channels - none
    % rejected

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
        if strcmp(channels{i,6},'good') | strcmp(channels{i,6},'bad') 
            chanidx(i) = 1;
        % else mark it for chucking
        else
           chanidx(i) = 0;
        end
    end
    % keep only channels marked for keeping
    channels = channels(logical(chanidx'),1);
    % get coordinates
    coords = table2cell(readtable(['sub-',p,'_electrodes.tsv'],'FileType','text','Delimiter','\t'));
    % get first column separately
    coords1 = coords(:,1);
    % then remove everything but the separate coordinate values and convert
    % contents of the cell array to double if necessary (coords remains as
    % a cell array).
    if ischar(coords{1,2}) % this asks what sort of object the first x coordinate is
        coords = num2cell(str2double(coords(:,2:4)));
    else
        coords = coords(:,2:4);
    end
    % initialise
    coords2keep = [];
    for i = 1:size(channels,1)
        % find index of electrode in coords
        rowidx = find(ismember(coords1,channels(i)));
        % find coordinates 
        coords2keep(i,:) = [coords{rowidx,1},coords{rowidx,2},coords{rowidx,3}];
    end
    % make metadata
    xyz = coords2keep;
    subject(1:size(xyz,1),1) = {['sub-',p]};
    ELECTRODE(str2num(p)) = {channels};
    % save variables 
    save([root,'/work/sub-',p,'/decode_during_preprocessing/metadata_input/coords/electrodes.mat'],'ELECTRODE','subject','xyz');
    spreadsheet = table;
    spreadsheet.subject = subject;
    spreadsheet.electrode = vertcat(ELECTRODE{:});
    spreadsheet.x = xyz(:,1);
    spreadsheet.y = xyz(:,2);
    spreadsheet.z = xyz(:,3);
    writetable(spreadsheet,[root,'/work/sub-',p,'/decode_during_preprocessing/metadata_input/coords/MNI_basal_electrodes_sub-01-22_w_label.csv'])

    % copy raw data into work subdirectory. The filename is 'filt_data.mat'
    % even though the data are not filtered! This is because, in the main
    % workflow, setup_data is run after filtering, so this is the filename
    % that setup_data is set up to expect.
    copyfile([root,'/data/sub-',p,'/ieeg/sub-',p,'_task-naming_run-01_ieeg.mat'],[root,'/work/sub-',p,'/decode_during_preprocessing/filt_data.mat'])
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
    'dataroot', [root,'/work/sub-',p,'/decode_during_preprocessing/'], ...
    'metaroot', [root,'/work/sub-',p,'/decode_during_preprocessing/metadata_input/'], ...
    'datarootout', ['/group/mlr-lab/Saskia/ECoG_LASSO/work/sub-',p,'/decode_during_preprocessing/'], ...
    'cvpath', [root,'/work/metadata_input/cv/cvpartition_10fold_DilkinaSplit.mat'], ...
    'overwrite', 1, ...
    'WriteIndividualMetadata', 1);

    % reorganise data ready for wavelet decomposition 
    a = load([root,'/work/sub-',p,'/decode_during_preprocessing/full/BoxCar/001/WindowStart/-1000/WindowSize/4000/s',p,'_.mat']);
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
    save([root,'/work/sub-',p,'/decode_during_preprocessing/data_raw.mat'],'X','-v7.3');
    disp('Done organising!')

    %% Part 1B - filtered data

    % copy across ready-filtered data from work directory (overwriting
    % previous file called filt_data.mat)
    copyfile([root,'/work/sub-',p,'/filt_data.mat'],[root,'/work/sub-',p,'/decode_during_preprocessing/filt_data.mat'])
    
    % using same metadata as before (i.e. retaining all ventral temporal
    % channels), epoch
    setup_data( ...
    'WindowStart', -1000, ... 
    'WindowSize', 4000, ...
    'BaselineWindow', 200, ...
    'subjects', str2num(p), ...
    'boxcar', 1, ...
    'slope_interval', 0, ...
    'average', 0, ...
    'datacode', '', ...
    'dataroot', [root,'/work/sub-',p,'/decode_during_preprocessing/'], ...
    'metaroot', [root,'/work/sub-',p,'/decode_during_preprocessing/metadata_input/'], ...
    'datarootout', ['/group/mlr-lab/Saskia/ECoG_LASSO/work/sub-',p,'/decode_during_preprocessing/'], ...
    'cvpath', [root,'/work/metadata_input/cv/cvpartition_10fold_DilkinaSplit.mat'], ...
    'overwrite', 1, ...
    'WriteIndividualMetadata', 1);

    % reorganise data ready for wavelet decomposition 
    a = load([root,'/work/sub-',p,'/decode_during_preprocessing/full/BoxCar/001/WindowStart/-1000/WindowSize/4000/s',p,'_.mat']);
    disp(['Organising data for subject ',p,' ...']);
    nelecs = size(a.X,2)/4000;
    X = reshape(a.X,400,4000,nelecs);
    X = permute(X,[3 2 1]);
    save([root,'/work/sub-',p,'/decode_during_preprocessing/data_filtered.mat'],'X','-v7.3');
    disp('Done organising!')

    %% Part 1C - data with bad channels rejected

    % using main metadata to automatically reject channels, epoch
    setup_data( ...
    'WindowStart', -1000, ... 
    'WindowSize', 4000, ...
    'BaselineWindow', 200, ...
    'subjects', str2num(p), ...
    'boxcar', 1, ...
    'slope_interval', 0, ...
    'average', 0, ...
    'datacode', '', ...
    'dataroot', [root,'/work/sub-',p,'/decode_during_preprocessing/'], ...
    'metaroot', [root,'/work/metadata_input/'], ...
    'datarootout', ['/group/mlr-lab/Saskia/ECoG_LASSO/work/sub-',p,'/decode_during_preprocessing/'], ...
    'cvpath', [root,'/work/metadata_input/cv/cvpartition_10fold_DilkinaSplit.mat'], ...
    'overwrite', 1, ...
    'WriteIndividualMetadata', 1);
    
    % reorganise data ready for wavelet decomposition 
    a = load([root,'/work/sub-',p,'/decode_during_preprocessing/full/BoxCar/001/WindowStart/-1000/WindowSize/4000/s',p,'_.mat']);
    disp(['Organising data for subject ',p,' ...']);
    nelecs = size(a.X,2)/4000;
    X = reshape(a.X,400,4000,nelecs);
    X = permute(X,[3 2 1]);
    save([root,'/work/sub-',p,'/decode_during_preprocessing/data_channelsrejected.mat'],'X','-v7.3');
    disp('Done organising!')

    %% Part 1D - data with common average referencing applied

    % Load data into EEGlab.
    load([root,'/work/sub-',num2str(p),'/data_reformat.mat'])
    load([root,'/scripts/eeglab_init.mat']); 
    EEG.setname = ['sub-',p];
    EEG.nbchan = size(X,1);
    EEG.data = X;
    EEG.srate = 1000; % data collected at 2000 Hz have now been downsampled to 1000 Hz
    EEG.subject = ['sub-',p];
    EEG = eeg_checkset(EEG);
    % apply common average reference
    EEG.data = reref(EEG.data,[]);

    % data are already organised for wavelet decomposition. Save
    disp(['Organising data for subject ',p,' ...']);
    X = EEG.data;
    save([root,'/work/sub-',p,'/decode_during_preprocessing/data_rereferenced.mat'],'X','-v7.3');
    disp('Done organising!')

   %% Note on part 1E - data with bad trials rejected

   % This will be computed during part 2!

   %% Part 2 - wavelet decomposition
   
   % loop over each dataset. Rereferenced data appears twice because we need 
   % a version before trial rejection and a version after trial rejection.
   filenames = {'raw','filtered','channelsrejected','rereferenced','rereferenced'};

   for f = 1:5
        
       % load data and metadata
       load([root,'/work/sub-',p,'/decode_during_preprocessing/data_',filenames{f},'.mat']);
       % if trials are to be rejected, load metadata
       if f == 5
           load([root,'/work/sub-',p,'/metadata.mat']);
       end
        
       %% Part 2A - Power
       
       % initialise variables
       power = zeros(size(X,3), size(X,1), 166, 60);
       baselinepower = zeros(size(X,3), size(X,1), 21, 60);

       % if trials are to be rejected, set bad trials to nans
       if f == 5
           X(:,:,~metadata.filters(10).filter)=NaN;
       end
       
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
       if f == 5
            save([root,'/derivatives/wavelet/sub-',p,'/trialsrejected_power.mat'],'power');
       else
            save([root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_power.mat'],'power');
       end

       % reformat into trials x spatiofrequency features
       tmp = permute(power,[1,3,2,4]);
       powersff = reshape(tmp,size(tmp,1),[]);
       
       % write .csv
       if f == 5
            writematrix(powersff,[root,'/derivatives/wavelet/sub-',p,'/trialsrejected_power.csv']);
       else
            writematrix(powersff,[root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_power.csv']);
       end

       % then power with decibel normalisation. Calculate condition-average
       % baseline
       baselinepower = squeeze(nanmean(baselinepower,[1,3]));
       % make it the same size as the data
       baselinepower = permute(repmat(baselinepower,1,1,100,166),[3,1,4,2]);
       % decibel-convert
       dBpower = 10*log10(power./baselinepower);
       % save
       if f == 5
            save([root,'/derivatives/wavelet/sub-',p,'/trialsrejected_dBpower.mat'],'dBpower');
       else
            save([root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_dBpower.mat'],'dBpower');
       end
       % reformat into trials x spatiofrequency features
       tmp = permute(dBpower,[1,3,2,4]);
       dBpowersff = reshape(tmp,size(tmp,1),[]);
       % write .csv
       if f == 5
            writematrix(dBpowersff,[root,'/derivatives/wavelet/sub-',p,'/trialsrejected_dBpower.csv']);
       else
            writematrix(dBpowersff,[root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_dBpower.csv']);
       end

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
      if f == 5
           save([root,'/derivatives/wavelet/sub-',p,'/trialsrejected_phase.mat'],'phase');
      else
           save([root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_phase.mat'],'phase');
      end
      % reformat into trials x spatiofrequency features
      tmp = permute(phase,[1,3,2,4]);
      phasesff = reshape(tmp,size(tmp,1),[]);
      % write .csv
      if f == 5
           writematrix(phasesff,[root,'/derivatives/wavelet/sub-',p,'/trialsrejected_phase.csv']);
      else
           writematrix(phasesff,[root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_phase.csv']);
      end

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
      if f == 5
            save([root,'/derivatives/wavelet/sub-',p,'/trialsrejected_voltage.mat'],'voltage');
      else
            save([root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_voltage.mat'],'voltage');
      end
      % reformat into trials x spatiofrequency features
      tmp = permute(voltage,[1,3,2]);
      voltagesff = reshape(tmp,size(tmp,1),[]);
      % write .csv
      if f == 5
           writematrix(voltagesff,[root,'/derivatives/wavelet/sub-',p,'/trialsrejected_voltage.csv']);
      else
           writematrix(voltagesff,[root,'/derivatives/wavelet/sub-',p,'/',filenames{f},'_voltage.csv']);
      end

    
   end

end








