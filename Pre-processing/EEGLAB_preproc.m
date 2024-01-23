
%% Preprocess - EEGLAB
% This is the FIRST section of data processing. 
%
% INPUT data for this script:
%   * raw EGG files collected by Biosemi device: 'xxx.bdf'
%
% TASKS this processing script should accomplish:
%   * Preprocess EEG data:
%       * Remove unsed channels + add location files 
%       * Resampling
%       * Re-reference
%       * Filting
%       * Epoch
%       * ICA
%       * Manual artifact rejection
%
% Main OUTPUT data for this script:
%   * Intermediate data in the intermidate dictionary: 'xxx.set'
%   * Final data saved in the output dictionary: 'xxx_final.set'
%
% Next step:
%   * Option1 - Preprocess for further statistical analysis by Fieldtrip:
%       * \1 ProcessingScripts\Script2_Preprocess_Fieldtrip.mat
%   * Option2 - Primary ERP visualization by EEGLAB:
%       * 


cd('D:\Jorge\Data')

% Set up global parameters

sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};

subnames = lower(sublists);

interpath  = 'D:\Jorge\Data\1.Raw_data\';%load raw data
outputpath = 'D:\Jorge\Data\EEGLAB_clean\';%store pre-proc data

addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\eeglab2021.1');
addpath('D:\Jorge\Data\1.Raw_data')
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\functions')
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\1 ProcessingScripts');

eeglab
close all;
%% Step 1 - Precleaning
session = 'PE';
for i = 1:25
    display(strcat({'load raw EEG in the '},session,{' file: '},sublists(i)));
    sublist = sublists(i);
    subname = subnames(i);
    
    EEGOUT = f_preclean(session,sublist,subname,interpath);
end
%% Step 2 - Run ICA

session = 'PE';
for i = 25
    % 2.1 Determine which subject
    display(strcat({'Start running ICA:'},sublists(i)));
    sublist = sublists(i);
    subname = subnames(i);
    
    % 2.2 Load EEG set file
    EEG = pop_loadset(strcat(interpath,char(sublist),'\', ...
        char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean.set'));
    
    % 2.3 Run ICA 
    [ALLEEG,EEGOUT_ICA,CURRENTSET] = processMARA([],EEG,[],[0,1,0,0,0]);
    
    EEG = pop_saveset (EEGOUT_ICA, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean_ICA.set'), ...
        'filepath',  strcat(interpath,char(sublist)));
    
end

%% Step 3 - Inspect ICA components

session = 'PE';

for i = 25

    % 3.1 Determine which subject
    display(strcat({'Inspect ICA artifacts:'},sublists(i)));
    sublist = sublists(i);
    subname = subnames(i);
    
    % 3.2 Load EEG file
    EEG = pop_loadset(strcat(interpath,char(sublist),'\',char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean_ICA.set'));

    % 3.3 Inspect ICA components
    EEG = pop_selectcomps(EEG, [1:35]);
    uiwait;
    EEG = pop_selectcomps(EEG, [36:size(EEG.icaact,1)]);
    uiwait;

    % 3.4 Remove ICA components
    EEGOUT_ICA_clean = pop_subcomp(EEG,[],1);

    % 3.6 Save data
     EEGOUT = pop_saveset (EEGOUT_ICA_clean, 'filename', ...
         strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean_ICA_clean.set'), 'filepath', strcat(interpath,char(sublist)));
end

%% Step 4 - Postcleaning
% add lines to post clean
session = 'PE';

for i = 25

    % 4.1 Determine which subject
    display(strcat({'Maunal artifact rejection:'},sublists(i)));
    sublist = sublists(i);
    subname = subnames(i);

    % 4.2 Manual artifact rejection
    EEGOUT = f_postclean(session,sublist,subname,interpath,outputpath);%add variables outputpath
end

%% Step 5 - ERP bin operation for PE and Pr session (by ERPLAB)

session   = 'PE';
inputpath = 'D:\1 Projects (offline)\1 Auditory Illusion Project\3 Data\3 AnalysisData\1 Preprocessing_eeglab\';

for i = 1:length(sublists)

    % Determine which subject
    display(strcat({'Operate Bins for:'},sublists(i)));
    sublist = sublists(i);
    subname = subnames(i);

    % Bin operation
    f_binoperateERP(session,sublist,subname,inputpath,outputpath)
end

%% Step 6 - Extract Mean ERP for each session per subject 

session = 'PE';
inputpath = 'D:\1 Projects (offline)\1 Auditory Illusion Project\3 Data\3 AnalysisData\1 Preprocessing_eeglab\';

for i = 1:length(sublists)
    % Determine which subject
    display(strcat({'Extract Mean ERP for:'},sublists(i)));
    sublist = sublists(i);
    subname = subnames(i);
    
    % 5.2 Extract mean ERP
    f_MeanERP(session,sublist,subname,inputpath,outputpath)
end