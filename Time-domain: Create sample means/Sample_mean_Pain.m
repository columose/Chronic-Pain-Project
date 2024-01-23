cd('D:\Jorge\Data\Field_trip_pre_process');

% 0. Define which condition
conditions = {'omi','omistd'} ;%or{'omi','std'}

sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};
subnames = lower(sublists);

inputpath  = 'D:\Jorge\Data\Field_trip_pre_process\';
outputpath = 'D:\Jorge\Data\mean_sample_data\';

addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\fieldtrip-20221126');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\functions');
%% 1. Load averaged ERP of individuals in two gruops

% 1.1 Create variables
% group set for further test
Omission_PE  = [];
Standard_PE = [];
for isub = 25
    display(['Loading omission condition ',sublists{isub}])

    % average data
    load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{1},'_preproc_ERP.mat'));
    
    %  data with all omission trials
    Omi_ERP_ST = load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{1},'_preproc_ERP_ST.mat'));
    Omi_ERP_ST = Omi_ERP_ST.data_ERP_ST;

    % display(['Loading standard condition ',sublists{isub}])
    % 
    % % data with all standard trials
    % std_ERP_ST = load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{2},'_preproc_ERP_ST.mat'));
    % std_ERP_ST = std_ERP_ST.data_ERP_ST;

    display(['Loading omission standard condition ',sublists{isub}])

    % data with all standard trials
    omi_std_ERP_ST = load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{2},'_preproc_ERP_ST.mat'));
    omi_std_ERP_ST = omi_std_ERP_ST.data_ERP_ST;

    % Calculate the mean of the sample mean for the standard condition
    omi_std_ERP = f_timelock_sampleMean(Omi_ERP_ST,omi_std_ERP_ST);
    omi_std_ERP.cfg = [];
    save(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{2},'_preproc_ERP_sampleMean'),'omi_std_ERP');

    % Omission_PE{isub}  = data_ERP;
    % Standard_PE{isub} = std_ERP;
end