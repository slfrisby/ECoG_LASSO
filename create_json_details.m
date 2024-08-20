% by Saskia. Creates participant-specific details for .json files made by
% createBIDS_ieeg_json.m .

subs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,21,22];

% the first 19 fields in the struct are naming data

json_details(1).sub_label = '01';
json_details(1).task_label = 'naming';
json_details(1).SamplingFrequency = 2000;
json_details(1).ECOGChannelCount = 61;
json_details(1).SEEGChannelCount = [];
json_details(1).EOGChannelCount = 1;
json_details(1).ECGChannelCount = [];
json_details(1).MiscChannelCount = [];
json_details(1).TriggerChannelCount = 3;
json_details(1).iEEGReference = 'mastoid';
json_details(1).SubjectArtefactDescription = '';

json_details(2).sub_label = '02';
json_details(2).task_label = 'naming';
json_details(2).SamplingFrequency = 1000;
json_details(2).ECOGChannelCount = 102;
json_details(2).SEEGChannelCount = [];
json_details(2).EOGChannelCount = 1;
json_details(2).ECGChannelCount = [];
json_details(2).MiscChannelCount = [];
json_details(2).TriggerChannelCount = 3;
json_details(2).iEEGReference = 'mastoid';
json_details(2).SubjectArtefactDescription = '';

json_details(3).sub_label = '03';
json_details(3).task_label = 'naming';
json_details(3).SamplingFrequency = 2000;
json_details(3).ECOGChannelCount = 61;
json_details(3).SEEGChannelCount = [];
json_details(3).EOGChannelCount = 1;
json_details(3).ECGChannelCount = [];
json_details(3).MiscChannelCount = [];
json_details(3).TriggerChannelCount = 3;
json_details(3).iEEGReference = 'mastoid';
json_details(3).SubjectArtefactDescription = '';

json_details(4).sub_label = '04';
json_details(4).task_label = 'naming';
json_details(4).SamplingFrequency = 1000;
json_details(4).ECOGChannelCount = 86;
json_details(4).SEEGChannelCount = [];
json_details(4).EOGChannelCount = [];
json_details(4).ECGChannelCount = [];
json_details(4).MiscChannelCount = 6;
json_details(4).TriggerChannelCount = 3;
json_details(4).iEEGReference = 'mastoid';
json_details(4).SubjectArtefactDescription = '';

json_details(5).sub_label = '05';
json_details(5).task_label = 'naming';
json_details(5).SamplingFrequency = 2000;
json_details(5).ECOGChannelCount = 56;
json_details(5).SEEGChannelCount = [];
json_details(5).EOGChannelCount = [];
json_details(5).ECGChannelCount = [];
json_details(5).MiscChannelCount = 5;
json_details(5).TriggerChannelCount = 2;
json_details(5).iEEGReference = 'mastoid';
json_details(5).SubjectArtefactDescription = '';

json_details(6).sub_label = '06';
json_details(6).task_label = 'naming';
json_details(6).SamplingFrequency = 1000;
json_details(6).ECOGChannelCount = 82;
json_details(6).SEEGChannelCount = [];
json_details(6).EOGChannelCount = 1;
json_details(6).ECGChannelCount = [];
json_details(6).MiscChannelCount = 1;
json_details(6).TriggerChannelCount = 2;
json_details(6).iEEGReference = 'mastoid';
json_details(6).SubjectArtefactDescription = '';

json_details(7).sub_label = '07';
json_details(7).task_label = 'naming';
json_details(7).SamplingFrequency = 1000;
json_details(7).ECOGChannelCount = 114;
json_details(7).SEEGChannelCount = [];
json_details(7).EOGChannelCount = [];
json_details(7).ECGChannelCount = [];
json_details(7).MiscChannelCount = 1;
json_details(7).TriggerChannelCount = 4;
json_details(7).iEEGReference = 'mastoid';
json_details(7).SubjectArtefactDescription = '';

json_details(8).sub_label = '08';
json_details(8).task_label = 'naming';
json_details(8).SamplingFrequency = 1000;
json_details(8).ECOGChannelCount = 94;
json_details(8).SEEGChannelCount = [];
json_details(8).EOGChannelCount = [];
json_details(8).ECGChannelCount = [];
json_details(8).MiscChannelCount = 4;
json_details(8).TriggerChannelCount = 3;
json_details(8).iEEGReference = 'mastoid';
json_details(8).SubjectArtefactDescription = '';

json_details(9).sub_label = '09';
json_details(9).task_label = 'naming';
json_details(9).SamplingFrequency = 1000;
json_details(9).ECOGChannelCount = 123;
json_details(9).SEEGChannelCount = [];
json_details(9).EOGChannelCount = [];
json_details(9).ECGChannelCount = [];
json_details(9).MiscChannelCount = 6;
json_details(9).TriggerChannelCount = 3;
json_details(9).iEEGReference = 'mastoid';
json_details(9).SubjectArtefactDescription =  '';

json_details(10).sub_label = '10';
json_details(10).task_label = 'naming';
json_details(10).SamplingFrequency = 1000;
json_details(10).ECOGChannelCount = 106;
json_details(10).SEEGChannelCount = [];
json_details(10).EOGChannelCount = [];
json_details(10).ECGChannelCount = 1;
json_details(10).MiscChannelCount = 8;
json_details(10).TriggerChannelCount = 3;
json_details(10).iEEGReference = 'mastoid';
json_details(10).SubjectArtefactDescription =  '';

json_details(11).sub_label = '11';
json_details(11).task_label = 'naming';
json_details(11).SamplingFrequency = 1000;
json_details(11).ECOGChannelCount = 106;
json_details(11).SEEGChannelCount = [];
json_details(11).EOGChannelCount = [];
json_details(11).ECGChannelCount = 1;
json_details(11).MiscChannelCount = 6;
json_details(11).TriggerChannelCount = 3;
json_details(11).iEEGReference = 'mastoid';
json_details(11).SubjectArtefactDescription =  '';

json_details(12).sub_label = '12';
json_details(12).task_label = 'naming';
json_details(12).SamplingFrequency = 1000;
json_details(12).ECOGChannelCount = 90;
json_details(12).SEEGChannelCount = [];
json_details(12).EOGChannelCount = [];
json_details(12).ECGChannelCount = [];
json_details(12).MiscChannelCount = 6;
json_details(12).TriggerChannelCount = 3;
json_details(12).iEEGReference = 'mastoid';
json_details(12).SubjectArtefactDescription =  '';

json_details(13).sub_label = '13';
json_details(13).task_label = 'naming';
json_details(13).SamplingFrequency = 1000; 
json_details(13).ECOGChannelCount = 82;
json_details(13).SEEGChannelCount = [];
json_details(13).EOGChannelCount = [];
json_details(13).ECGChannelCount = 1;
json_details(13).MiscChannelCount = 5;
json_details(13).TriggerChannelCount = 5;
json_details(13).iEEGReference = 'mastoid'
json_details(13).SubjectArtefactDescription =  '';

json_details(14).sub_label = '14';
json_details(14).task_label = 'naming';
json_details(14).SamplingFrequency = 2000;
json_details(14).ECOGChannelCount = 92;
json_details(14).SEEGChannelCount = [];
json_details(14).EOGChannelCount = [];
json_details(14).ECGChannelCount = 1;
json_details(14).MiscChannelCount = 4;
json_details(14).TriggerChannelCount = 4;
json_details(14).iEEGReference = 'mastoid';
json_details(14).SubjectArtefactDescription =  '';

json_details(15).sub_label = '15';
json_details(15).task_label = 'naming';
json_details(15).SamplingFrequency = 2000;
json_details(15).ECOGChannelCount = 102;
json_details(15).SEEGChannelCount = [];
json_details(15).EOGChannelCount = [];
json_details(15).ECGChannelCount = 1;
json_details(15).MiscChannelCount = 4;
json_details(15).TriggerChannelCount = 3;
json_details(15).iEEGReference = 'mastoid';
json_details(15).SubjectArtefactDescription =  '';

json_details(16).sub_label = '17';
json_details(16).task_label = 'naming';
json_details(16).SamplingFrequency = 2000;
json_details(16).ECOGChannelCount = 128;
json_details(16).SEEGChannelCount = [];
json_details(16).EOGChannelCount = [];
json_details(16).ECGChannelCount = 1;
json_details(16).MiscChannelCount = 6;
json_details(16).TriggerChannelCount = 3;
json_details(16).iEEGReference = 'mastoid';
json_details(16).SubjectArtefactDescription =  '';

json_details(17).sub_label = '20';
json_details(17).task_label = 'naming';
json_details(17).SamplingFrequency = 2000;
json_details(17).ECOGChannelCount = 106;
json_details(17).SEEGChannelCount = [];
json_details(17).EOGChannelCount = [];
json_details(17).ECGChannelCount = 1;
json_details(17).MiscChannelCount = [];
json_details(17).TriggerChannelCount = 4;
json_details(17).iEEGReference = 'mastoid';
json_details(17).SubjectArtefactDescription =  '';

json_details(18).sub_label = '21'; 
json_details(18).task_label = 'naming';
json_details(18).SamplingFrequency = 2000;
json_details(18).ECOGChannelCount = 98; 
json_details(18).SEEGChannelCount = [];
json_details(18).EOGChannelCount = [];
json_details(18).ECGChannelCount = 1;
json_details(18).MiscChannelCount = 4;
json_details(18).TriggerChannelCount = 4;
json_details(18).iEEGReference = 'mastoid';
json_details(18).SubjectArtefactDescription =  '';

json_details(19).sub_label = '22';
json_details(19).task_label = 'naming';
json_details(19).SamplingFrequency = 2000;
json_details(19).ECOGChannelCount = 92;
json_details(19).SEEGChannelCount = 24;
json_details(19).EOGChannelCount = 6;
json_details(19).ECGChannelCount = 1;
json_details(19).MiscChannelCount = 4;
json_details(19).TriggerChannelCount = 4;
json_details(19).iEEGReference = 'mastoid'
json_details(19).SubjectArtefactDescription =  '';

% the next 3 fields are semantic judgement

json_details(20).sub_label = '13';
json_details(20).task_label = 'semanticjudgement';
json_details(20).SamplingFrequency = 1000; 
json_details(20).ECOGChannelCount = 82;
json_details(20).SEEGChannelCount = [];
json_details(20).EOGChannelCount = [];
json_details(20).ECGChannelCount = 1;
json_details(20).MiscChannelCount = 5;
json_details(20).TriggerChannelCount = 5;
json_details(20).iEEGReference = 'mastoid'
json_details(20).SubjectArtefactDescription =  '';

json_details(21).sub_label = '14';
json_details(21).task_label = 'semanticjudgement';
json_details(21).SamplingFrequency = 2000;
json_details(21).ECOGChannelCount = 92;
json_details(21).SEEGChannelCount = [];
json_details(21).EOGChannelCount = [];
json_details(21).ECGChannelCount = 1;
json_details(21).MiscChannelCount = 4;
json_details(21).TriggerChannelCount = 4;
json_details(21).iEEGReference = 'mastoid';
json_details(21).SubjectArtefactDescription =  '';

json_details(22).sub_label = '15';
json_details(22).task_label = 'semanticjudgement';
json_details(22).SamplingFrequency = 2000;
json_details(22).ECOGChannelCount = 102;
json_details(22).SEEGChannelCount = [];
json_details(22).EOGChannelCount = [];
json_details(22).ECGChannelCount = 1;
json_details(22).MiscChannelCount = 4;
json_details(22).TriggerChannelCount = 3;
json_details(22).iEEGReference = 'mastoid';
json_details(22).SubjectArtefactDescription =  '';

% save

save('/group/mlr-lab/Saskia/ECoG_LASSO/scripts/json_details.mat','json_details');