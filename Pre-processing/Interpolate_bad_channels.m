cd('D:\Jorge\Data')
% eeglab
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
%%
session = 'PE';
for i = 25
    % 2.1 Determine which subject
    display(strcat({'Start interpolating bad channels:'},sublists(i)));
    sublist = sublists(i);
    subname = subnames(i);
 % 2.2 Load EEG set file
    EEG_all_chans = pop_loadset(strcat(interpath,char(sublist),'\', ...
        char(subname),'_',session,'_chanloc.set'));
    
   all_chans=EEG_all_chans.chanlocs;

   EEG_final= pop_loadset(strcat(interpath,char(sublist),'\', ...
        char(subname),'_',session,'_final.set'));

 % 12.1 (IF BAD CHANNELS) Interpolate bad channels
    % if isempty(badchan_index) == 0 

        EEGOUT_interp = interpol(EEG_final, EEG_all_chans.chanlocs, 'spherical');
        EEGOUT = pop_saveset (EEGOUT_interp, 'filename', ...
            strcat(char(subname),'_',session,'_interp.set'), 'filepath', strcat(interpath,char(sublist)));
 end
