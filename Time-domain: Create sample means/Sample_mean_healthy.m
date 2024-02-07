cd('D:\Jorge\Data\Fan_healthy_data')

% 0. Define which condition
conditions = {'omi','std'};  % {'omi','omistd'} or{'omi','std'}


sublists  = {'SDT001','AAT002','CAS003','MUS004','OET005','EAW006', ...
            'NIS007','CUD008','PEB009','EOS012','RNZ014','AMY015', ...
            'PBC016','CHY017','CRY018','TAA019','CAS020','EUS021','GYW022', ...
            'ROR024','LAS026','JUT028','DIT029','AOR030','SOS031', ...
            'LAS032','SRS033','OOS034','ELR035','SES036','DUS037','MRV038', ...
            'LAS039','SAW040','MIR041','CEX043','SOS044','REZ045','LAS046'};
subnames=lower(sublists);


inputpath  = 'D:\Jorge\Data\Fan_healthy_data\';
outputpath = 'D:\Jorge\Data\Fan_healthy_data\';

addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\fieldtrip-20221126');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\functions');

%% 1. Load averaged ERP of individuals in two gruops

% 1.1 Create variables
% group set for further test
load('Deviant_PE.mat')
load('Standard_PE.mat')
%%
for iSub=2:length(sublists)
    std_ERP_healthy=[];
    std_ERP_healthy{iSub,1}=f_timelock_sampleMean(Deviant_PE{1,iSub},Standard_PE{1,iSub});
    std_ERP_healthy{iSub,1}.cfg=[];
    save(strcat(inputpath,sublists{iSub},'\',subnames{iSub},'_',conditions{2},'_sampleMean_healthy'),'std_ERP_healthy');
end

%%

for isub = 1:length(sublists)
    display(['Loading omission condition ',sublists{isub}])

    % average data
    load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{1},'_preproc_ERP.mat'));
    %  data with all omission trials
    Omi_ERP_ST = load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{1},'_preproc_ERP_ST.mat'));
    Omi_ERP_ST = Omi_ERP_ST.data_ERP_ST;

    display(['Loading standard condition ',sublists{isub}])

    % data with all standard trials
    std_ERP_ST = load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{2},'_preproc_ERP_ST.mat'));
    std_ERP_ST = std_ERP_ST.data_ERP_ST;

    % display(['Loading omission standard condition ',sublists{isub}])
    % 
    % % data with all standard trials
    % omi_std_ERP_ST = load(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{3},'_preproc_ERP_ST.mat'));
    % omi_std_ERP_ST = omi_std_ERP_ST.data_ERP_ST;

    % Calculate the mean of the sample mean for the standard condition
    std_ERP = f_timelock_sampleMean(Omi_ERP_ST,std_ERP_ST);
    std_ERP.cfg = [];
    save(strcat(inputpath,sublists{isub},'\',subnames{isub},'_',conditions{2},'_preproc_ERP_sampleMean'),'std_ERP');

    % Omission_PE{isub}  = data_ERP;
    % Standard_PE{isub} = std_ERP;
end
