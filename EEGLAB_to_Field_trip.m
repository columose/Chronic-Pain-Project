% Set up global parameters

sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};


subnames = lower(sublists);
inputpath= 'D:\Jorge\Data\1.Raw_data\';
outputpath='D:\Jorge\Data\Field_trip_pre_process\';

addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\eeglab2021.1');
addpath('D:\Jorge\Data\1.Raw_data')
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\functions')
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\1 ProcessingScripts');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\fieldtrip-20221126')
% Determine the session and conditions
session    = 'PE';
conditions = {'std','dev','omi','omistd'};
%% Step 1 - Read data and Generate full ERP for specific condition per subject
for icon = 1
    for i = 1
        
        disp(strcat({'Working in subject '} ,subnames(i), {' in the '},conditions(icon),{' condition.'}));

        % 1. Read data from EEGLAB to fieldtrip format
        cd(strcat(inputpath,sublists{i}));
        cfg = [];
        cfg.dataset = strcat(subnames{i},'_',session,'_interp.set');
        
        %determine which paradigm
        strcmp(session,'PE')
        cfg = f_trldef_pe(cfg,conditions{icon});
      
        data_all = ft_preprocessing(cfg);
        % save(strcat(outputpath,sublists{i},'\',subnames{i},'_',conditions{icon},'_preproc'),'data_all');
        
        % 2. Generation of full ERP ST
        cfg.latency = 'all';
        cfg.keeptrials = 'yes';
        data_ERP_ST = ft_timelockanalysis(cfg, data_all);
        % save(strcat(outputpath,sublists{i},'\',subnames{i},'_',conditions{icon},'_preproc_ERP_ST'),'data_ERP_ST');
        
        % 3. Generation of full ERP
        cfg.latency = 'all';
        cfg.keeptrials = 'no';
        data_ERP = ft_timelockanalysis(cfg, data_all);
        % save(strcat(outputpath,sublists{i},'\',subnames{i},'_',conditions{icon},'_preproc_ERP'),'data_ERP');
        
    end
end