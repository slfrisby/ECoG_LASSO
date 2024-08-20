function make_coefficient_volumes()
% Takes .csv files of coefficients and constructs volumes. These volumes
% can then be projected onto the surface and used to make plots and
% videos. This is the first stage in the coefficient plotting pipeline (the
% next is make_coefficient_surface_mesh.sh).

    addpath('/group/mlr-lab/AH/Projects/spm12/')
    root = '/group/mlr-lab/Saskia/ECoG_LASSO';
    cd(root);
    
    bands = {'dB_power','dB_theta_power','dB_alpha_power','dB_beta_power','dB_gamma_power','dB_high_gamma_power','voltage'};

    % load mni template provided by Riki and Akihiro
    % gunzip([root,'/scripts/data_info/MNI152_T1_2mm_brain_nocereb.nii.gz']);
    % vol = spm_vol([root,'/scripts/data_info/MNI152_T1_2mm_brain_nocereb.nii']);
    % bug fix - use spare functional volume in same MNI space as template.
    % For whatever reason, this works while the template does not
    vol = spm_vol([root,'/scripts/tmpT1.nii.gz']);
    template = spm_read_vols(vol);
    % also get header information
    T = vol.mat;

    % load mni coordinates
    coords = readtable([root,'/scripts/mni_coordinates.csv']);
    coords = table2array(coords);
    % delete rows with nans
    % coords = rmmissing(coords,1);
    % create index of which coordinates are LH (0) and which are RH (1)
    hindex = coords(:,2) > 0;
    % transform coordinates into matrix space
    coords(:,2:4) = mni2cor(coords(:,2:4),T);

    % get participant index. Participant 20 is excluded because we have no
    % coordinates for them
    subs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 22];

    % for each band
    for b = 1:length(bands)

        % load coefficients as double
        coefs = table2array(readtable([root,'/work/coefficients/R/',bands{b},'_coefficients.csv']));

        % for every participant
        for s = 1:length(subs)
            
            if ~exist([root,'/work/sub-',sprintf('%02d',subs(s)),'/coefficients/volume/'])
                mkdir([root,'/work/sub-',sprintf('%02d',subs(s)),'/coefficients/volume/']);
            end
            
            
            % get that participant's coordinates
            scoords = coords(coords(:,1)==subs(s),:);
            scoefs = coefs(coords(:,1)==subs(s),:);
            shindex = hindex(coords(:,1)==subs(s),:);

            % remove nans
            nanidx = isnan(scoords);
            nanidx = logical(any(nanidx,2));
            scoords(nanidx,:) = [];
            scoefs(nanidx,:) = [];
            shindex(nanidx,:) = [];
        
            % loop over timepoints
            for tp = 1:size(scoefs,2) 
               % initialise volumes
               LH = zeros(size(template));
               RH = zeros(size(template));
               % loop over electrodes
               for elec = 1:size(scoords,1) 
                   % if the coordinates are in the left hemisphere
                   if ~(shindex(elec))
                        % add coefficients to LH image
                        LH(scoords(elec,2),scoords(elec,3),scoords(elec,4)) = scoefs(elec,tp);
                   else
                       % otherwise add them to the RH image
                       RH(scoords(elec,2),scoords(elec,3),scoords(elec,4)) = scoefs(elec,tp);
                   end
               end
            
               % save volume
               % if the electrodes are only in the left hemisphere
               if sum(shindex)==0
                   % save the left hemisphere
                   tmp = vol;
                   tmp.fname = [root,'/work/sub-',sprintf('%02d',subs(s)),'/coefficients/volume/L_sub-',sprintf('%02d',subs(s)),'_',bands{b},'_timepoint_',sprintf('%03d',tp),'_coefficients.nii'];
                   spm_write_vol(tmp,LH);
               % if they are only in the right hemisphere
               elseif sum(shindex)>0
                   % save the right hemisphere
                   tmp = vol;
                   tmp.fname = [root,'/work/sub-',sprintf('%02d',subs(s)),'/coefficients/volume/R_sub-',sprintf('%02d',subs(s)),'_',bands{b},'_timepoint_',sprintf('%03d',tp),'_coefficients.nii'];
                   spm_write_vol(tmp,RH);
               end
            end
        end
    end
end