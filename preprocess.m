function preprocess(p)

% by Saskia. This function does the hands-free section of the preprocessing
% pipeline for the ECoG time-frequency analysis. It should be followed by
% the manual section - inspect_data.m . 

    addpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/')
    addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/eeglab2023.1'))
    % eeglab should have the SASICA plugin installed. To adapt this
    % function for ECoG, open eeg_SASICA.m and comment out lines 148-150
    % and 710-818.
    addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/ECoG_Data_Prep-master'))

    root = '/group/mlr-lab/Saskia/ECoG_LASSO';
    cd(root);
 
    % make work directory
    mkdir([root,'/work/sub-',p]);

    % tell us where it is up to
    disp(['Working on subject ',p,' ...']);
    
    %% STEP 1 - Filter data to remove line noise and slow drifts

    % load sidecar .json, which contains sampling rate.
    tmp = fileread([root,'/data/sub-',p,'/ieeg/sub-',p,'_task-naming_run-01_ieeg.json']);
    jsonfile = jsondecode(tmp);

    % Load data. The data comes in the form of a struct containing both data
    % and stimulus onset "tags". Load this whole struct into an arbitrary
    % variable.
    a = load([root,'/data/sub-',p,'/ieeg/sub-',p,'_task-naming_run-01_ieeg.mat']);

    % Find the name(s) of the data structs within this master struct. In
    % some cases there are multiple data structs.
    structs2load = {};
    tmp = fieldnames(a);
    for i = 1:length(tmp)
        if isstruct(a.(tmp{i}))
            structs2load = [structs2load,tmp(i)];       
        end      
    end
    
    % Then load in the structs one at a time.
    for t = 1:length(structs2load)

        % EEGlab opens a GUI on startup. However, on the compute nodes of the high-performance
        % cluster it is not possible to open a GUI. Rather than starting EEGlab
        % in the normal way, we create empty EEGlab structures...
        load([root,'/scripts/eeglab_init.mat']); 
        % ... and populate them...
        EEG.setname = ['sub-',p];
        EEG.nbchan = size(a.(structs2load{t}).DATA,2);
        EEG.data = a.(structs2load{t}).DATA';
        EEG.srate = jsonfile.SamplingFrequency; 
        EEG.subject = ['sub-',p];
        % ...and use the toolbox's built-in function to check that the
        % structures are set up correctly.
        EEG = eeg_checkset(EEG); 
        
        % To bring up a GUI (when working on the login node, rather than
        % the submit node):
        % eeglab redraw

        % plot channel spectra before any filtering. Set frequency range to
        % be the range that we will use for wavelet convolution.
        % EEG.pnts*EEG.srate gives the length of data acquisition in seconds.
        f = figure('visible','off'); 
        f = pop_spectopo(EEG, 1, [0 EEG.pnts/EEG.srate*1000], 'EEG' , 'freqrange',[4 200],'electrodes','off');
        saveas(gcf,[root,'/work/sub-',p,'/raw_',structs2load{t},'.png'])
        close(gcf)

        % 1A - Filter out line noise made by the ECoG recording equipment and other 
        % nearby appliances. Line noise tends to occur at either 50 Hz or 60 Hz.
        % Power line noise in Japan is at 60 Hz, so we will filter at 60 Hz and its 
        % harmonics 120 Hz and 180 Hz (Clarke 2020). Here, this is done using the EEGlab plugin CleanLine.
        % Information about CleanLine can be found here:
        % https://github.com/sccn/cleanline . NOTE: Cleanline is not 100% 
        % effective at removing line noise, but tends not to remove any
        % signal, whereas other methods (notch filtering, Zapline-plus) are
        % more aggressive. Common average referencing, further down the
        % pipeline, removes even more line noise (online help pages suggest
        % that one should do this before applying CleanLine, but my own
        % tests show that the order makes little difference in these data.)
        %%%
        disp('Filtering...')
        EEG = pop_cleanline(EEG,'LineFrequencies',[60 120 180],'ScanForLines',0,'SignalType','Channels','ChanCompIndices',[1:EEG.nbchan],'SlidingWinLength',4,'SlidingWinStep',2,'ComputeSpectralPower',0,'VerboseOutput',0);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');
        
        % plot channel spectra after CleanLine.
        f = figure('visible','off'); 
        f = pop_spectopo(EEG, 1, [0 EEG.pnts/EEG.srate*1000], 'EEG' , 'freqrange',[4 200],'electrodes','off');
        saveas(gcf,[root,'/work/sub-',p,'/cleanline_',structs2load{t},'.png'])
        close(gcf)

        % 1B - Apply a low-pass filter at 0.5 Hz to filter out slow drifts. While this does
        % not make any difference to the time-frequency results themselves, it will
        % make ERP extraction and comparison between frequency results and ERPs
        % easier later on. This must be done before epoching because the edge
        % artifact of this kind of filter can be up to 6 seconds long, which would
        % be very problematic.
        
        % and
        
        % 1C - Apply a high-pass filter at 300 Hz. For participants that are
        % sampled at 1000 Hz, the recording equipment imposes this filter anyway.
        % (For participants sampled at 2000 Hz, the filter is 600 Hz). Very high frequency activity is
        % not of interest to the study and tends to be artifactual.
        EEG = pop_eegfiltnew(EEG,'locutoff',0.5,'hicutoff',300,'plotfreqz',0);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on');

        % We need to remove the data from EEGlab and save it so that Chris's
        % windowing scripts can get to it. We also need to make sure that the other
        % important info that Chris's scripts need are saved where they need to be.
        data = EEG.data';
        a.(structs2load{t}).DATA = data;

    end
    
    % Save and clear
    save([root,'/work/sub-',p,'/filt_data.mat'],'-struct', 'a', '-v7.3')
    clearvars -except jsonfile p root

    %% STEP 2 - Downsample, epoch (into one big epoch, rather than many little windows), reject channels below the seizure onset zone, and baseline-correct

    % These are all achieved using Chris Cox's function setup_data. The version
    % deployed here is edited to suit these data. The original can be
    % downloaded from https://github.com/crcox/ECoG_Data_Prep. 

    % make metadata input if necessary
    if ~exist([root,'/work/metadata_input'])
        make_metadata_input
    end

    % Parameters:

    % WindowStart and WindowSize
    % We wish to generate one big window per trial (rather than a moving 
    % or opening window). The window of interest for later analysis is 
    % 0:1650 ms (where 0 is stimulus onset) but we need a buffer so 
    % that any edge artefacts lie outside the window of interest. We 
    % also want to avoid overlap because we want to run ICA later. We 
    % therefore select -1000:3000 ms.

    % BaselineWindow
    % This is set to match the voltage analysis of the same data. It is
    % not perfect because the previous stimulus was still onscreen.
    % Additionally, some advise the use of a baseline that stops well
    % before stimulus onset to avoid bleeding of activity into the
    % pretrial period. However, it is important that these results are
    % comparable with previous work, so we select a baseline of -200 ms
    % until stimulus onset.

    % subjects
    % We wish to proces one dataset at a time. 

    % boxcar
    % This defines the number of MILISECONDS over which the data is
    % averaged. For example, if the data were sampled at 2000 Hz and we 
    % had the following timeseries:
    % [2, 4, 6, 8];
    % A boxcar average of 1 ms would down-sample to 1000 Hz and yield:
    % [3, 7];
    % Note that, because this parameter is defined in miliseconds
    % rather than in data points, the effect depends on the sampling
    % rate. We set boxcar = 1 so that data collected at 2000 Hz will be
    % downsampled to match data collected at 1000 Hz (which will be
    % unaffected).

    % slope_interval
    % If slope interval > 0, results describe rate of change in voltage rather
    % than voltage value.

    % average
    % 'False' means that data from all 4 presentations of the same
    % stimulus will not be averaged at this stage.

    % datacode
    % Chris's data have a suffix '_raw' or '_ref'.

    % dataroot
    % The data are found in the work directory for each participant.

    % metaroot
    % This points the script towards important metadata.

    % datarootout
    % We wish to output the data to the work folder for each
    % participant.

    % cvpath
    % This points the script towards cross-validation schemes that are
    % stored as part of the metadata (though not used later on).

    % overwrite
    % 1 for true.

    % writeIndividualMetadata
    % Since we are running the analysis individually, we will save
    % metadata individually. (Note that for CHTC we shoud save it all
    % together.)
    
    setup_data( ...
    'WindowStart', -1000, ... 
    'WindowSize', 4000, ...
    'BaselineWindow', 200, ...
    'subjects', str2num(p), ...
    'boxcar', 1, ...
    'slope_interval', 0, ...
    'average', 0, ...
    'datacode', '', ...
    'dataroot', [root,'/work/sub-',p], ...
    'metaroot', [root,'/work/metadata_input/'], ...
    'datarootout', [root,'/work/sub-',p], ...
    'cvpath', [root,'/work/metadata_input/cv/cvpartition_10fold_DilkinaSplit.mat'], ...
    'overwrite', 1, ...
    'WriteIndividualMetadata', 1);

    % Reformat the data so that it can be loaded back into EEGlab
    
    a = load([root,'/work/sub-',p,'/full/BoxCar/001/WindowStart/-1000/WindowSize/4000/s',p,'_.mat']);
    
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
    save([root,'/work/sub-',p,'/data_reformat.mat'],'X','-v7.3');
    disp('Done organising!')
    clearvars -except jsonfile p root

end








