clear;close all;clc;
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\eeglab2021.1_old')
inputpath  = 'D:\Jorge\Data\1.Raw_data\';
cd(inputpath)

 sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};
 subnames = lower(sublists);

 binlabels = {'B1(condition1)','B2(condition2)','B3(condition3)','B4(condition4)'};
 resultpath='D:\Jorge\Data\Source_localisation_data\Pain\';

COMP=[3,1;3,4;2,1];
PE={'sPE','cPE','pPE'};

 brainstorm
 eeglab
 close all
 
% For indexing in brainstorm
for isub=1:length(sublists)
    sublist=sublists{isub};
    subname=subnames{isub};

    RawFiles{isub}=strcat(inputpath,char(sublist),'\',char(subname),'_PE_interp.set');%Create path to rawfile data


    % Run analysis using brainstorm
    % Input files
    sFiles = [];
    SubjectNames = subnames;
    
    % Process: Import MEG/EEG: Events
        sFiles=[];
        sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
            'subjectname',   SubjectNames{isub}, ...
            'condition',     '', ...
            'datafile',      {{RawFiles{isub}}, 'EEG-EEGLAB'}, ...
            'eventname',     'B1(condition1), B2(condition2), B3(condition3),B4(condition4)', ...
            'timewindow',    [-0.1, 1.848], ...
            'epochtime',     [-0.1, 1.848], ...
            'split',         0, ...
            'createcond',    1, ...
            'ignoreshort',   1, ...
            'channelalign',  1, ...
            'usectfcomp',    1, ...
            'usessp',        1, ...
            'freq',          [], ...
            'baseline',      [], ...
            'blsensortypes', 'MEG, EEG');


    % Process: Refine registration
    sFiles = bst_process('CallProcess', 'process_headpoints_refine', sFiles, [], ...
        'tolerance', 0);
    
    % Process: Project electrodes on scalp
    sFiles = bst_process('CallProcess', 'process_channel_project', sFiles, [], ...
        'sensortypes', 'EEG');
        
    % Process: Average: By trial group (folder average)
    sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
        'avgtype',       5, ...  % By trial group (folder average)
        'avg_func',      1, ...  % Arithmetic average:  mean(x)
        'weighted',      0, ...
        'keepevents',    0);

    sFiles_avg=sFiles;
    
    % Input files
    sFiles = sFiles_avg(3);%Use condition with fewest trials to compute a headmodel
  
    % Process: Compute head model
    sFiles = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
        'Comment',     '', ...
        'sourcespace', 1, ...  % Cortex surface
        'meg',         1, ...  % 
        'eeg',         3, ...  % OpenMEEG BEM
        'ecog',        1, ...  % 
        'seeg',        1, ...  % 
        'openmeeg',    struct(...
             'BemSelect',    [1, 1, 1], ...
             'BemCond',      [1, 0.0125, 1], ...
             'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
             'BemFiles',     {{}}, ...
             'isAdjoint',    0, ...
             'isAdaptative', 1, ...
             'isSplit',      0, ...
             'SplitLength',  4000), ...
        'channelfile', '');

    temp = bst_get('Study'); 
    headmodel = temp.HeadModel.FileName;

    %Copy this head model to other condition folders
    db_set_headmodel(headmodel,'AllConditions');

    temp=[];
    headmodel=[];

    %Change sFiles back to average of conditions for the covariance matrices
    sFiles=sFiles_avg;

    % Process: Compute covariance (noise or data)
    sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
        'baseline',       [-0.1, 0], ...
        'datatimewindow', [0, 1.848], ...
        'sensortypes',    ' EEG', ...
        'target',         1, ...  % Noise covariance     (covariance over baseline time window)
        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);  % Replace
    
    % Process: Compute sources [2018]
    sFiles = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
        'output',  1, ...  % Kernel only: shared
        'inverse', struct(...
             'Comment',        'sLORETA: EEG', ...
             'InverseMethod',  'minnorm', ...
             'InverseMeasure', 'sloreta', ...
             'SourceOrient',   {{'fixed'}}, ...
             'Loose',          0.2, ...
             'UseDepth',       0, ...
             'WeightExp',      0.5, ...
             'WeightLimit',    10, ...
             'NoiseMethod',    'reg', ...
             'NoiseReg',       0.1, ...
             'SnrMethod',      'fixed', ...
             'SnrRms',         1e-06, ...
             'SnrFixed',       3, ...
             'ComputeKernel',  1, ...
             'DataTypes',      {{'EEG'}}));

        %Find file source filenames for difference ERPs
        for icon=1:4
            sFiles_source(icon).file=sFiles(icon).FileName;
        end
    
        % Obtain source localisation differences between our comparisons of interest
        for icomp=1:3
            sFiles=[];
            sFiles{icomp}= sFiles_source(COMP(icomp,1)).file;
            sFiles2{icomp}=sFiles_source(COMP(icomp,2)).file;
        
    
            % Process: Difference: A-B
            sFiles = bst_process('CallProcess', 'process_diff_ab', sFiles{icomp}, sFiles2{icomp}, ...
                'source_abs', 0);
            
            % Process: Set name: deviant - standard
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           PE{icomp}, ...
                'isindex',       1);
            
            % Process: Low-pass:40Hz
            sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
                'highpass',    0, ...
                'lowpass',     40, ...
                'tranband',    0, ...
                'attenuation', 'strict', ...  % 60dB
                'ver',         '2019', ...  % 2019
                'mirror',      0, ...
                'overwrite',   1);
            save(strcat(resultpath,sublists{isub},'\',subnames{isub},'_',PE{icomp}),'sFiles')
            
            % Process: Z-score transformation: [-100ms,0ms]
            sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
                'baseline',   [-0.1, 0], ...
                'source_abs', 0, ...
                'method',     'zscore', ...  % Z-score transformation:    x_std = (x - &mu;) / &sigma;
                'overwrite',  0);
        end

end
%% Perform grand-averaging

% Obtain file link to be used for grand-averaging
sFiles=[];
GA_input=[];
for isub=1:length(sublists)
    for icomp=1:3
        temp=[];
        temp=load(strcat(resultpath,sublists{isub},'\',subnames{isub},'_',PE{icomp},'.mat'));

        sFiles{icomp,isub}=temp.sFiles;
        GA_input{icomp,isub}=sFiles{icomp,isub}.FileName;%links for individual files to be used for grand-averaging
    end
end

% Perform grand-averaging
for icomp=1:3
    % Process: Average: By trial group (grand average)
    sFiles = bst_process('CallProcess', 'process_average', {GA_input{icomp,:}}, [], ...
    'avgtype',         7, ...  % By trial group (grand average)
    'avg_func',        2, ...  % Average absolute values:  mean(abs(x))
    'weighted',        0, ...
    'scalenormalized', 0);
end


panel_time('SetCurrentTime',1.5)%useful for visualisation
%% Group difference
% Input files
sFiles=[];
sFiles = {...
    'Group_analysis/Pain_difference/results_average_230717_1039.mat', ...
    'Group_analysis/Pain_difference/results_average_230717_1041.mat', ...
    'Group_analysis/Pain_difference/results_average_230717_1040.mat'};
sFiles2 = {...
    'Group_analysis/Controls_difference/results_average_230717_1041.mat', ...
    'Group_analysis/Controls_difference/results_average_230717_1043.mat', ...
    'Group_analysis/Controls_difference/results_average_230717_1042.mat'};

% Start a new report
bst_report('Start', sFiles);

% Process: Difference: A-B, abs
sFiles = bst_process('CallProcess', 'process_diff_ab', sFiles, sFiles2, ...
    'source_abs', 1);

% Process: Low-pass:40Hz
sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
    'highpass',    0, ...
    'lowpass',     40, ...
    'tranband',    0, ...
    'attenuation', 'strict', ...  % 60dB
    'ver',         '2019', ...  % 2019
    'mirror',      0, ...
    'overwrite',   1);
%% Run permutation test. All time samples (avg time = 0)
sFiles=[];
sFiles2=[];
TOI=[1.432 1.488; 1.730 1.840; 1.432 1.510];

cd(inputpath)
load("GA_input_pain.mat")

cd('D:\Jorge\Data\Fan_healthy_data\All_paradigm_data\1 Preprocessing_eeglab\');
load("GA_input_controls.mat")

for icomp=1:3

    % Input files
    sFiles = {GA_input_pain{icomp,:}};
    sFiles2 = {GA_input_controls{icomp,:}};
    
    % Process: Perm absmean test [1432ms,1488ms]          H0:(|mean(A)|=|mean(B)|), H1:(|mean(A)|<>|mean(B)|)
    sFiles = bst_process('CallProcess', 'process_test_permutation2', sFiles, sFiles2, ...
        'timewindow',     [TOI(icomp,1) TOI(icomp,2)], ...
        'scoutsel',       {'Desikan-Killiany', {'bankssts L', 'bankssts R', 'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'frontalpole L', 'frontalpole R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'temporalpole L', 'temporalpole R', 'transversetemporal L', 'transversetemporal R'}}, ...
        'scoutfunc',      1, ...  % Mean
        'isnorm',         0, ...
        'avgtime',        1, ...
        'iszerobad',      1, ...
        'Comment',        '', ...
        'test_type',      'absmean', ...  % Absolute mean test:    (works with unconstrained sources)T = (|mean(A)|-|mean(B)|) / sqrt(|var(A)|/nA + |var(B)|/nB)
        'randomizations', 1000, ...
        'tail',           'two');  % Two-tailed
end
%% Extract one value from each participant averaged over TOI for each comparison
GA_all=[GA_input_pain GA_input_controls];
temp=erase(GA_all,'.mat');
GA_all_z=strcat(temp,'_zscore.mat');


for icomp=1:3
    sFiles=[];
    sFiles={GA_all_z{icomp,:}};
    % Process: Extract values: [1432ms,1488ms] 68 scouts abs
    sFiles = bst_process('CallProcess', 'process_extract_values', sFiles, [], ...
    'timewindow', [TOI(icomp,1) TOI(icomp,2)], ...
    'scoutsel',   {'Desikan-Killiany', {'bankssts L', 'bankssts R', 'caudalanteriorcingulate L', 'caudalanteriorcingulate R', 'caudalmiddlefrontal L', 'caudalmiddlefrontal R', 'cuneus L', 'cuneus R', 'entorhinal L', 'entorhinal R', 'frontalpole L', 'frontalpole R', 'fusiform L', 'fusiform R', 'inferiorparietal L', 'inferiorparietal R', 'inferiortemporal L', 'inferiortemporal R', 'insula L', 'insula R', 'isthmuscingulate L', 'isthmuscingulate R', 'lateraloccipital L', 'lateraloccipital R', 'lateralorbitofrontal L', 'lateralorbitofrontal R', 'lingual L', 'lingual R', 'medialorbitofrontal L', 'medialorbitofrontal R', 'middletemporal L', 'middletemporal R', 'paracentral L', 'paracentral R', 'parahippocampal L', 'parahippocampal R', 'parsopercularis L', 'parsopercularis R', 'parsorbitalis L', 'parsorbitalis R', 'parstriangularis L', 'parstriangularis R', 'pericalcarine L', 'pericalcarine R', 'postcentral L', 'postcentral R', 'posteriorcingulate L', 'posteriorcingulate R', 'precentral L', 'precentral R', 'precuneus L', 'precuneus R', 'rostralanteriorcingulate L', 'rostralanteriorcingulate R', 'rostralmiddlefrontal L', 'rostralmiddlefrontal R', 'superiorfrontal L', 'superiorfrontal R', 'superiorparietal L', 'superiorparietal R', 'superiortemporal L', 'superiortemporal R', 'supramarginal L', 'supramarginal R', 'temporalpole L', 'temporalpole R', 'transversetemporal L', 'transversetemporal R'}}, ...
    'scoutfunc',  1, ...  % Mean
    'isnorm',     1, ...
    'avgtime',    1, ...
    'dim',        2, ...  % Concatenate time (dimension 2)
    'Comment',    '');

end
%%
fields= cellstr(sPE_TOI.Description);

for i=1:68
    fields{i}=strrep(fields{i},' ','');
end
Values=sPE_TOI.Value';
fields=fields';




