% adapted for this dataset by Saskia

%% Template Matlab script to create an BIDS compatible _ieeg.json file
% This example lists all required and optional fields.
% When adding additional metadata please use CamelCase
%
% DHermes, 2017
% modified Jaap van der Aar & Giulio Castegnaro 30.11.18

% Writing json files relies on the JSONio library
% https://github.com/bids-standard/bids-matlab
%
% Make sure it is in the matab/octave path

addpath('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/bids-matlab-main')

try
    bids.bids_matlab_version;
catch
    warning('%s\n%s\n%s\n%s', ...
            'Writing the JSON file seems to have failed.', ...
            'Make sure that the following library is in the matlab/octave path:', ...
            'https://github.com/bids-standard/bids-matlab');
end

%% 

% load struct of data that is specific to each participant
load('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/json_details.mat')

% loop through patients
for q = 1:size(json_details,2)
    
    clearvars -except q json_details

    % this_dir = fileparts(mfilename('fullpath'));
    root_dir = '/group/mlr-lab/Saskia';

    project = 'ECoG_LASSO/data';

    sub_label = json_details(q).sub_label;
    task_label = json_details(q).task_label;
    run_label = '01';

    name_spec.modality = 'ieeg';
    name_spec.suffix = 'ieeg';
    name_spec.ext = '.json';
    name_spec.entities = struct('sub', sub_label, ...
                                'task', task_label, ...
                                'run', run_label);

    % using the 'use_schema', true
    % ensures that the entities will be in the correct order
    bids_file = bids.File(name_spec, 'use_schema', true);

    % Contrust the fullpath version of the filename
    json_name = fullfile(root_dir, project, bids_file.bids_path, bids_file.filename);

    %% General fields, shared with MRI BIDS and MEG BIDS

    % to get the definition of each column,
    % you can use the bids.Schema class from bids matlab
    % For example
    schema = bids.Schema;
    def = schema.get_definition('iEEGCoordinateUnits');
    fprintf(def.description);

    %% Required fields:
    json.TaskName = task_label;
    json.SamplingFrequency = json_details(q).SamplingFrequency;
    json.PowerLineFrequency = 60;
    json.SoftwareFilters = '';

    %% Recommended fields:
    HardwareFilters.HighpassFilter.CutoffFrequency = [];
    HardwareFilters.LowpassFilter.CutoffFrequency = [];
    json.HardwareFilters = HardwareFilters; %

    json.Manufacturer = 'Ad-Tech Medical';
    json.ManufacturersModelName = '';
    json.TaskDescription = '';
    json.Instructions = '';
    json.CogAtlasID = 'https://www.cognitiveatlas.org/FIXME';
    json.CogPOID = 'http://www.cogpo.org/ontologies/CogPOver1.owl#FIXME';
    json.InstitutionName = '';
    json.InstitutionAddress = '';
    json.DeviceSerialNumber = '';
    json.ECOGChannelCount = json_details(q).ECOGChannelCount;
    json.SEEGChannelCount = json_details(q).SEEGChannelCount;
    json.EEGChannelCount = [];
    json.EOGChannelCount = json_details(q).EOGChannelCount;
    json.ECGChannelCount = json_details(q).ECGChannelCount;
    json.EMGChannelCount = [];
    json.MiscChannelCount = json_details(q).MiscChannelCount;
    json.TriggerChannelCount = json_details(q).TriggerChannelCount;
    json.RecordingDuration = [];
    json.RecordingType = '';
    json.EpochLength = [];
    json.SubjectArtefactDescription = '';
    json.SoftwareVersions = '';

    %% Specific iEEG fields:

    % If mixed types of references, manufacturers or electrodes are used, please
    % specify in the corresponding table in the _electrodes.tsv file

    %% Required fields:
    json.iEEGReference = json_details(q).iEEGReference;

    %% Recommended fields:
    json.ElectrodeManufacturer = 'Ad-Tech Medical';
    json.ElectrodeManufacturersModelName = '';
    json.iEEGGround = '';
    json.iEEGPlacementScheme = '';
    json.iEEGElectrodeGroups = '';

    %% Optional fields:
    json.ElectricalStimulation = '';
    json.ElectricalStimulationParameters = '';

    %% Write JSON
    % Make sure the directory exists
    bids.util.mkdir(fileparts(json_name));
    bids.util.jsonencode(json_name, json);

end