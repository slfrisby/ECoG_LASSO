% by Saskia. Enables interactive channel and trial rejection.
% This script is NOT designed for the cluster compute nodes.

addpath(genpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/eeglab2023.1'))
% eeglab should have the SASICA plugin installed. To adapt this
% function for ECoG, open eeg_SASICA.m and comment out lines 148-150
% and 710-818.
root = '/group/mlr-lab/Saskia/ECoG_LASSO';

subs = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','17','20','21','22'};
badtrials = struct;
badtrials(1).subject = 1;
badtrials(1).trials = [34,98,142,149,179,180,194,222,224,294,308,393];
badtrials(2).subject = 2;
badtrials(2).trials = [24,80,83,113,117,164,165,175,209,225,233,254,261,270,296,307,308,339,340,341,344,354,359,364,365,366,372,380,390,391];
badtrials(3).subject = 3;
badtrials(3).trials = [269,346];
badtrials(4).subject = 4;
badtrials(4).trials = [6,8,22,48,51,75,78,89,96,97,102,118,122,134,143,148,154,157,160,167,169,177,179,184,205,208,217,223,236,238,241,245,251,265,285,305,307,311,315,324,341,348,360,388,400];
badtrials(5).subject = 5;
badtrials(5).trials = [6,29,30,31,32,34,64,91,102,160,222,231,273,380,381,389];
badtrials(6).subject = 6;
badtrials(6).trials = [28,52,65,69,83,94,119,127,130,143,150,155,156,164,166,168,179,185,194,197,207,211,236,246,264,285,288,300,305,309,313,317,319,325,327,329,332,334,343,344,358,365,369,382,386,394,396,399];
badtrials(7).subject = 7;
badtrials(7).trials = [65,72,80,90,105,110,117,146,164,168,170,171,179,181,187,202,203,204,207,210,220,221,222,228,229,233,234,235,247,249,251,253,266,267,273,277,278,279,283,286,288,299,305,306,307,309,310,318,320,321,325,327,333,336,338,342,345,350,358,359,363,371,372,374,379,381,384,385,397,399];
badtrials(8).subject = 8;
badtrials(8).trials = [10,16,19,35,36,37,43,48,53,57,65,66,70,71,73,75,77,81,82,83,84,87,89,91,94,95,98,104,105,107,108,109,114,115,116,117,119,127,128,129,130,131,132,133,134,135,136,139,141,147,148,152,153,158,161,162,164,171,175,177,182,184,186,188,191,192,193,197,199,200,201,203,205,207,208,209,210,211,214,218,220,221,222,224,227,232,233,234,236,239,242,244,245,247,248,249,250,251,253,254,258,263,273,276,279,280,282,283,288,289,291,292,293,296,297,299,301,303,304,305,306,307,310,315,316,320,327,328,329,330,331,332,333,336,337,338,339,341,342,343,346,349,350,353,354,355,357,358,359,362,363,364,366,367,368,369,370,371,371,374,375,381,383,384,385,388,389,391,393,395,399,400];
badtrials(9).subject = 9;
badtrials(9).trials = [1,40,46,50,60,82,108,114,124,144,155,174,199,208,212,218,331,337,343,353,357,371,378,382,388,396,399];
badtrials(10).subject = 10;
badtrials(10).trials = [27,45,64,137,231,235,329,363,371];
badtrials(11).subject = 11;
badtrials(11).trials = [6,17,19,21,32,43,44,48,55,56,58,59,63,64,73,81,86,86,90,92,99,104,111,120,133,145,152,163,168,188,197,200,250,263,281,307,310,371,398];
badtrials(12).subject = 12;
badtrials(12).trials = [5,7,14,18,19,23,25,26,34,35,36,41,45,47,50,59,63,67,70,73,74,75,77,81,82,86,91,93,98,100,103,107,108,110,111,125,131,146,149,153,162,164,165,171,186,197,212,227,230,236,237,240,242,247,249,254,256,261,265,266,267,278,281,282,291,292,294,295,301,303,306,307,316,317,318,323,328,334,335,339,344,348,365,374,380,384,386,395,398,399];
badtrials(13).subject = 13;
badtrials(13).trials = [4,11,12,26,29,32,36,45,54,63,64,77,88,92,93,102,105,106,117,133,138,141,149,150,154,156,165,167,168,174,187,191,199,200,201,202,209,210,217,218,219,239,243,250,253,276,280,288,289,303,306,307,308,309,310,313,314,324,329,339,351,356,359,367,372,376,378];
badtrials(14).subject = 14;
badtrials(14).trials = [8,9,20,39,58,65,73,74,95,106,107,110,115,125,126,127,148,160,171,177,182,184,202,216,220,221,227,230,234,258,261,268,271,275,276,308,312,317,318,320,321,375,376,377,380,388];
badtrials(15).subject = 15;
badtrials(15).trials = [6,9,21,33,36,53,61,63,80,81,88,118,120,122,123,128,133,142,167,182,184,192,194,197,207,208,216,225,248,255,273,274,276,280,307,308,313,320,323,324,331,356,359,373,378,391,397,398];
badtrials(16).subject = 17;
badtrials(16).trials = [3,5,7,8,11,12,13,14,16,19,23,26,27,28,29,30,31,33,34,36,37,41,42,45,46,47,48,49,51,53,55,57,58,60,61,63,64,65,67,68,69,70,71,72,73,74,79,82,83,84,85,89,90,91,92,93,95,97,98,99,101,102,103,104,106,107,108,109,110,112,113,114,115,117,118,121,124,125,126,127,128,129,130,132,135,136,138,140,141,142,143,144,150,151,152,153,154,158,160,161,162,164,165,166,169,170,171,174,176,177,179,181,182,186,189,190,191,193,194,197,199,200,202,203,204,205,206,207,209,211,213,215,216,218,219,220,221,222,223,225,226,227,230,231,233,234,235,236,237,238,240,242,243,246,247,249,250,251,252,253,254,255,257,261,262,263,264,265,266,267,270,271,272,273,274,275,276,277,278,279,281,286,287,289,290,293,294,295,296,297,299,300,301,303,304,305,307,308,309,310,311,312,313,315,316,317,318,319,320,321,322,323,324,325,326,327,330,331,332,333,334,335,336,337,339,340,343,344,345,347,348,349,351,353,354,355,356,357,359,360,364,366,367,368,369,371,372,373,374,376,378,379,380,381,388,389,390,391,392,393,394,395,397,399];
badtrials(17).subject = 20;
badtrials(17).trials = [89,152,177,206,211,212,230,234,236,245,249,251,259,282,297,326,327,328,330,339,343,369,371,393,396,400];
badtrials(18).subject = 21;
badtrials(18).trials = [7,8,9,20,32,34,39,43,47,49,50,53,62,64,65,69,76,78,79,82,86,94,99,100,105,108,110,113,118,125,130,132,135,137,141,147,149,151,152,153,165,173,174,175,177,179,188,197,203,207,208,212,213,217,222,223,227,230,237,238,242,243,247,254,256,257,259,261,264,265,269,271,272,276,282,283,284,285,291,299,302,304,306,310,313,315,316,317,319,322,326,324,327,332,334,337,339,345,346,348,350,353,354,363,368,380,386,399];
badtrials(19).subject = 22;
badtrials(19).trials = [110,125,126,215,229,318,387];
% Participants 7 and 8 have obviously autocorrelated components, but these
% have been retained because the component clearly contains signal as well
% as noise. 

for q = 1:length(subs)

    p = subs{q};

    % Load data into EEGlab.
    load([root,'/work/sub-',num2str(p),'/data_reformat.mat'])
    load([root,'/scripts/eeglab_init.mat']); 
    EEG.setname = ['sub-',p];
    EEG.nbchan = size(X,1);
    EEG.data = X;
    EEG.srate = 1000; % data collected at 2000 Hz have now been downsampled to 1000 Hz
    EEG.subject = ['sub-',p];
    EEG = eeg_checkset(EEG);
    % eeglab redraw

    % Create a filter for trials.
    trialfilter = ones(400,1);

    %% Stage 1 - common average referencing
    
    % retain old data in case it is useful to plot
    mastoiddata = EEG.data;
    % compute and apply common average reference
    EEG.data = reref(EEG.data,[]);
    % eegplot(EEG.data)

    % to plot new data on top of old
    % eegplot(mastoiddata,'data2',EEG.data)

    %% Stage 2 - automatic rejection with a liberal threshold
    % Auto-reject any trials containing a deviation of over 10 standard
    % deviations from the channel mean (over all timepoints and trials). 
    % This accounts for the fact that some participants have more variable
    % data than others and this variation does not necessarily reflect
    % abnormal activity. (This is different to the threshold in the eLife
    % paper, but for most participants the threshold ends up similar in
    % practice.)

    % initialise variables
    trials2rej = [];
    channelspikecount = [];

    % for each channel
    for c = 1:size(EEG.data,1)
        % get data and calculate the mean and standard deviation for that
        % channel
        channeldata = squeeze(EEG.data(c,:,:));
        m = mean(channeldata,'all');
        s = std(channeldata,0,'all');
        % find the indices of the trials where there are spikes over 10
        % standard deviations above or below the mean
        [~,spike] = ind2sub(size(channeldata),find(channeldata > m+10*s | channeldata < m-10*s));
        % store these in a running total for that participant
        trials2rej = vertcat(trials2rej,spike);
        % store the number of spikes per participant.
        channelspikecount(c) = length(unique(spike));
    end

    % Some trials contain spikes in multiple channels. We need only count
    % each trial once.
    trials2rej = unique(trials2rej);

    % Plot a figure to assess the distribution of spikes across channels.
    % Sometimes there is one channel that has an atypically large number of
    % spikes. If this is the case, it can be better to reject that channel
    % rather than all those trials.
    f = figure;
    f = bar(1:size(EEG.data,1),squeeze(channelspikecount));
    xticks(1:size(EEG.data,1))
    xlabel('Electrode')
    ylabel('Number of outliers')
    title(['Sub-',p])
    saveas(gcf,[root,'/work/sub-',p,'/spike_distribution.png'])
    % close(gcf)   

    % mark those trials for rejection
    trialfilter(trials2rej) = 0;

    %% Stage 3 - manual trial rejection
    
    % eegplot(EEG.data)

    % Set the y-axis to 500. Inspect the trials that are auto-rejected, to ensure that this is sensible.
    % Then make a note of any trials that:
         % - appear to contain ictal or interictal activity
         %   (https://www.learningeeg.com/epileptiform-activity). Be careful of wicket spikes, which look
         %   periodic but are actually a feature of the healthy temporal lobe.
         % - appear to contain artifacts, especially muscle activity. N.B. do not panic about
         %   muscle activity, as much of this is removed by common average referencing. Also note 
         %   electrode pop.

    % Fill in the variable badtrials by hand. 
    trialfilter(badtrials(q).trials) = 0;

    %% Stage 4 - Independent components analysis (ICA)
    
    % retain a copy of the EEG data
    completedata = EEG.data;
    % retain data shape to help with data reshaping later
    [s1, s2, s3] = size(EEG.data);
    % reshape data ready for ICA (i.e. concatenate trials to make one big
    % timeseries).
    icadat = EEG.data;
    % before ICA we reject trials that are flagged for rejection. However,
    % we do not reject these trials from the main dataset yet because it will put the trial labels out of sync. 
    icadat(:,:,trialfilter==0) = [];
    icadat = reshape(icadat, size(icadat,1), []);
    % Load ICA data into EEGlab
    EEG.data = icadat; 

    % Conduct ICA. The number of ICs is 75% the number of good electrodes
    % (Clarke 2020).
    [weights,sphere,compvars,bias,signs,lrates,ICs] = runica(icadat,'extended',1,'PCA',round(size(icadat,1)*0.75));

    save([root,'/work/sub-',p,'/ICA.mat'],'weights','sphere','compvars','bias','signs','lrates','ICs','-v7.3')

    EEG.icaweights=weights;
    EEG.icasphere=sphere;
    EEG.icawinv=pinv(EEG.icaweights*EEG.icasphere); % see https://sccn.ucsd.edu/pipermail/eeglablist/2009/002907.html to explain why this is so
    EEG.icaact = reshape(ICs,round(size(icadat,1)*0.75),s2,[]);
    EEG = eeg_checkset(EEG);
    % eeglab redraw
   
    %% 4.1 - SASICA

    % We use SASICA to flag components with low autocorrelation
    % (usually represents muscle activity) and/or high focal trial acrivity
    % (usually represents muscle activity, electrode pop, or similar).

    % configure SASICA
    cfg = SASICA('getdefs');
    cfg.autocorr.enable = logical(1);
    cfg.trialfoc.enable= logical(1);
    % apply SASICA
    [EEG,cfg] = eeg_SASICA(EEG,cfg);
    % list trials with low autocorrelation or high focal trial activity
    flagcomps = find(EEG.reject.gcompreject);

    % Scroll component activations (view 20 trials at once for ease). See
    % whether you agree that the component looks noisy or has some unusual
    % trials. If you don't agree, don't reject the component. For focal
    % trial activity, consider rejecting the focal trials rather than the
    % whole component. NOTE: because only the (assumed) good trials have
    % been loaded into EEGlab, the trial indices displayed by EEGlab will
    % not be the true trial indices in the complete data. To find the true
    % trial indices to store in badtrials, use, for example:
    % tmp = EEG.data(1,1,345)
    % [x,y,z] = ind2sub(size(completedata),find(completedata==tmp))
    % z is the trial index in the complete data. If there are multiple
    % values of z, the correct value is the one where x = 1 and y = 1.

    % eegplot(EEG.data)

    %% 4.2 - Microsaccade filtering

    % This is based on the method of Clarke (2020). It works by:
        % 1. Performing ICA
        % 2. Filtering independent components for gamma/high gamma (the range in which we
        % expect to see microsaccadic activity)
        % 3. Convolving a template of saccade-related activity with the ICA
        % components
        % 4. Plotting a graph showing which components have the highest
        % number of saccade events. If it is clear that one or two
        % components have MANY more saccade events than others, these
        % components should be rejected.

    % Load the independent components into the data space of eeglab in order to use the filter  
    EEG.data = ICs;
    EEG = eeg_checkset(EEG);
    % filter the components for gamma and high-gamma activity
    EEG = pop_eegfiltnew(EEG,'locutoff',20,'hicutoff',190,'plotfreqz',0);
    fICs = EEG.data;
    fICs = reshape(fICs, size(fICs,1), []);

    % initialise variables
    mpd=20; % Minimum peak distance for saccades - originally mpd=round(20/(1000/srate)).
    locs = []; % this will contain the locations of the saccades within the IC
    z=[]; % this will contain the saccade-template-filtered ICs

    % for each component
    for i = 1:size(fICs,1)
        % filter the component with a saccade template
        z(i,:) = filtSRP(double(fICs(i,:))',1000);
        [pks,locs] = findpeaks(z(i,:),'minpeakheight',2*mean(abs(z(i,:)),2),'minpeakdistance',mpd);
        sac_rate(i) = length(locs)/(length(z)*1000); % number of events per second
    end
    % sort components by number of saccades
    sortsacs = sortrows([1:length(sac_rate); sac_rate]',2,'descend')';
    % plot a figure to assess whether some components have dramatically
    % more microsaccades than others
    f = figure;
    f = bar(sortsacs(2,:));
    set(gca,'xtick',[])
    xlabel('Component')
    ylabel('Number of saccades')
    title(['Sub-',p])
    saveas(gcf,[root,'/work/sub-',p,'/microsaccade_count.png'])
    close(gcf)  

    %% IF YOU AGREE THAT COMPONENTS SHOULD BE REJECTED - reject components

    % % list components
    % comps2keep = 1:size(ICs,1);
    % % mark flagged components for removal
    % comps2keep(ismember(comps2keep,flagcomps))=[];
    % icadat = reshape(completedata,size(completedata,1),[]);
    % % reject and reconstruct
    % [icaprojdata] = icaproj(icadat,weights,comps2keep);
    % % reshape
    % X = reshape(icaprojdata,s1,s2,s3);

    %% IF NO COMPONENTS NEED REJECTING
    X = completedata;


    %% Save
    save([root,'/work/sub-',p,'/clean.mat'],'X');
    mkdir([root,'/inspected/sub-',p,'/'])
    save([root,'/inspected/sub-',p,'/clean.mat'],'X');

    %% Save filter as part of metadata
    load([root,'/work/sub-',p,'/full/BoxCar/001/WindowStart/-1000/WindowSize/4000/metadata__',p,'.mat'])
    metadata.filters(10).label = 'badtrials';
    metadata.filters(10).dimension = 1;
    metadata.filters(10).filter = logical(trialfilter);
    save([root,'/work/sub-',p,'/metadata.mat'],'metadata');

    % Throw a warning if there is not at least one good trial for each
    % stimulus.
    detector = trialfilter(1:100)+trialfilter(101:200)+trialfilter(201:300)+trialfilter(301:400);
    livingdetector = trialfilter(1:50)+trialfilter(101:150)+trialfilter(201:250)+trialfilter(301:350);
    nonlivingdetector = trialfilter(51:100)+trialfilter(151:200)+trialfilter(251:300)+trialfilter(351:400);
    if ~isempty(find(detector==0))
        warning(['participant ',p,' is missing ',num2str(length(find(livingdetector==0))),' living trials and ',num2str(length(find(nonlivingdetector==0))),' nonliving trials.'])
    end

    %% Clear
    clearvars -except root subs badtrials
end
