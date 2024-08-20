function wavelet(p)

    % by Saskia. This script takes the preprocessed data (run through both
    % preprocess.m and inspect_data.m) and performs complex Morlet wavelet
    % convolution, extracting both power and phase.
    
    addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/eeglab2023.1'))
    root = '/group/mlr-lab/Saskia/ECoG_LASSO';
    
    % load data and metadata
    load([root,'/work/sub-',p,'/clean.mat']);
    load([root,'/work/sub-',p,'/metadata.mat']);

    disp(['Working on patient ',p,' ...'])
    
    %% Part 1 - Power

    % initialise output - electrodes x timepoints (166 between 0 and 1650) x
    % trials x frequencies (60, logarithmically spaced between 4 and 200
    % Hz)
    power = zeros(size(X,3), size(X,1), 166, 60);

    % also initialise a variable to store baseline power
    baselinepower = zeros(size(X,3), size(X,1), 21, 60);
    
    % set bad trials to nans
    X(:,:,~metadata.filters(10).filter)=NaN;

    % We are interested in the average power over repeated presentations of
    % the same stimulus. Because of this, we conduct wavelet decomposition,
    % extract power, and THEN average.

    % loop over electrodes
    for elec = 1:size(X,1)
        
        % Perform time-frequency analysis with the following parameters:
        % Sampling rate = 1000 Hz
        % Number of wavelets = 60
        % Frequency range = 4 - 200 Hz. In the results, 4 Hz is the first entry and 200 Hz is the 60th entry. (Clarke 2020 uses 190 Hz as the upper
        % cutoff).
        % Frequency spacing = logarithmic
        % Number of wavelet cycles increases with frequency to increase
        % accuracy. 5-cycle wavelet used at 4 Hz, increasing to a 15-cycle
        % wavelet at 15 Hz.
        % Times = 0 : 1650 ms in 50 ms intervals (since the current epoch is -1000:2999, we
        % extract times 1001, 1051, 1100, ..., 2651)

        [output, freqs, times] = timefreq(squeeze(X(elec,:,:)), 1000, 'cycles', [5 15], 'freqs', [4 200], 'nfreqs', 60, 'freqscale', 'log', 'timesout',[1001:10:2651]);
        
        % Also produce output for the baseline period. The
        % baseline-correction applied during preprocessing is -200:-1 ms.
        % Here we apply a baseline of -300:-100 ms - there is a gap before
        % 0 ms to mitigate the effct of temporal leakage from the trial
        % into the baseline (although some of this may still be present),
        % the baseline is as short as possible (which is desirable because
        % there is no gap between stimuli, so we want to avoid capturing
        % signal from the previous stimulus as far as possible), and it is
        % long enough to capture an entire cycle at 4 Hz, the lowest
        % frequency of interest.
        [baselineoutput, freqs, times] = timefreq(squeeze(X(elec,:,:)), 1000, 'cycles', [5 15], 'freqs', [4 200], 'nfreqs', 60, 'freqscale', 'log', 'timesout',[701:10:901]);
    
        % Calculate power. (Multiplying a number by its complex conjugate gives
        % its magnitude squared. Magnitude in this case is amplitude, and
        % amplitude squared is power.)
        pow = output.*conj(output);
        baselinepow = baselineoutput.*conj(baselineoutput);

        % fill in results and baseline 
        power(:,elec,:,:) = permute(pow,[3,2,1]);
        baselinepower(:,elec,:,:) = permute(baselinepow,[3,2,1]);

    end

    %% Part 1A - power without baseline correction

    % calculate mean over items, excluding nans
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
    mkdir([root,'/derivatives/wavelet/sub-',p]);
    save([root,'/derivatives/wavelet/sub-',p,'/power.mat'],'power');

    % also save a log-transformed version for visualisation purposes
    logpower = log(power);
    save([root,'/derivatives/wavelet/sub-',p,'/logpower.mat'],'logpower');

    % reformat the data into a matrix of trials x spatiofrequency features
    % (electrode 1 at timepoint 1 and frequency 1, electrode 1 at timepoint
    % 2 and frequency 1, ... electrode 1 at timepoint 166 and frequency 1,
    % electrode 2 at timepoint 1 and frequency 1, ... electrode n at
    % timepoint 166 and frequency 1, electrode 1 at timepoint 1 and
    % frequency 2, ...)
    tmp = permute(power,[1,3,2,4]);
    powersff = reshape(tmp,size(tmp,1),[]);

    % write .csv
    writematrix(powersff,[root,'/derivatives/wavelet/sub-',p,'/power.csv']);

    %% Part 1B - power with decibel normalisation

    % We will use a condition-average baseline (in this case, this means
    % that we average the baseline across all trials). This increases the
    % signal-to-noise of the baseline power and means that results will not
    % be sensitive to trial-to-trial differences in baseline. The baseline
    % is, however, electrode- and frequency-specific.

    % calculate condition-average baseline
    baselinepower = squeeze(nanmean(baselinepower,[1,3]));

    % make it the same size as the data
    baselinepower = permute(repmat(baselinepower,1,1,100,166),[3,1,4,2]);

    % decibel-convert
    dBpower = 10*log10(power./baselinepower);

    % save
    save([root,'/derivatives/wavelet/sub-',p,'/dBpower.mat'],'dBpower');

    % reformat into trials x spatiofrequency features
    tmp = permute(dBpower,[1,3,2,4]);
    dBpowersff = reshape(tmp,size(tmp,1),[]);

    % write .csv
    writematrix(dBpowersff,[root,'/derivatives/wavelet/sub-',p,'/dBpower.csv']);

   %% Part 2 - Phase

   % For power calculations, we are more interested in the average power
   % than in the power of the average - so we conduct wavelet decomposition
   % and THEN average across trials. However, phase values cannot be
   % averaged, so for phase we MUST calculate an average first and then
   % wavelet. (This approach is also used by Clarke 2020).

   % initialise output
   phase = zeros(size(X,3)/4, size(X,1), 166, 60);

   % average data
   tmp = cat(4,X(:,:,1:100),X(:,:,101:200),X(:,:,201:300),X(:,:,301:400));
   avX = nanmean(tmp,4);

   % loop over electrodes
   for elec = 1:size(X,1)

        % wavelet
        [output, freqs, times] = timefreq(squeeze(avX(elec,:,:)), 1000, 'cycles', [5 15], 'freqs', [4 200], 'nfreqs', 60, 'freqscale', 'log', 'timesout',[1001:10:2651]);
        % Calculate phase. (This is the angle of the complex number). 
        pha = angle(output);
        % fill in results
        phase(:,elec,:,:) = permute(pha,[3,2,1]);

   end
   
    % median-impute remaining NaNs
    tmp = find(isnan(phase));
    [trial,elec,window,freq] = ind2sub(size(phase),tmp);
    for i = 1:length(trial)
        phase(trial(i),elec(i),window(i),freq(i)) = nanmedian(phase(:,elec(i),window(i),freq(i)));
    end
    
    % save
    save([root,'/derivatives/wavelet/sub-',p,'/phase.mat'],'phase');
    
    % reformat into trials x spatiofrequency features
    tmp = permute(phase,[1,3,2,4]);
    phasesff = reshape(tmp,size(tmp,1),[]);
    
    % write .csv
    writematrix(phasesff,[root,'/derivatives/wavelet/sub-',p,'/phase.csv']);

    %% Part 3 - voltage

    % These are simply ERPs - averaged over repeated presentations of the
    % same stimulus. They differ from the voltages in Rogers et al. (2021)
    % because of the differences in preprocessing.

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
    save([root,'/derivatives/wavelet/sub-',p,'/voltage.mat'],'voltage');
    
    % reformat into trials x spatiofrequency features
    tmp = permute(voltage,[1,3,2]);
    voltagesff = reshape(tmp,size(tmp,1),[]);
    
    % write .csv
    writematrix(voltagesff,[root,'/derivatives/wavelet/sub-',p,'/voltage.csv']);

end

