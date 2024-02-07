% This script is for comparing sPE and cPE ERPS between the chronic pain
% and healthy participant groups
datapath='D:\Jorge\Data\Field_trip_pre_process\';
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\fieldtrip-20221126');

sublists_chronic = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};
subnames_chronic = lower(sublists_chronic);

conditions = {'std','dev','omi','omistd'};

%Load single trial data for condition
for i=1:length(sublists_chronic)
    cd(strcat(datapath,sublists_chronic{i}));
    temp_std=load(strcat('eeg',num2str(i),'_',char(conditions(1)),'_preproc_ERP_sampleMean.mat'));
    temp_dev=load(strcat('eeg',num2str(i),'_',char(conditions(2)),'_preproc_ERP.mat'));
    temp_omi=load(strcat('eeg',num2str(i),'_',char(conditions(3)),'_preproc_ERP.mat'));
    temp_omistd=load(strcat('eeg',num2str(i),'_',char(conditions(4)),'_preproc_ERP_sampleMean.mat'));
   
    STD_all{i,1}=temp_std.std_ERP;
    DEV_all{i,1}=temp_dev.data_ERP;
    OMI_all{i,1}=temp_omi.data_ERP;
    OMI_STD_all{i,1}=temp_omistd.omi_std_ERP;

    temp_std=[];
    temp_dev=[];
    temp_omi=[];
    temp_omistd=[];

    % Then take the difference of the averages using ft_math
    cfg           = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    sPE_ERP_chronic{i,1} = ft_math(cfg, OMI_all{i,1}, STD_all{i,1});%Stimulus-driven Prediction error variable 
    cPE_ERP_chronic{i,1} = ft_math(cfg, OMI_all{i,1}, OMI_STD_all{i,1});%Context-driven Prediction error variable 
    pPE_ERP_chronic{i,1} = ft_math(cfg, DEV_all{i,:},STD_all{i,:});%Probability-driven variable is deviant- standard 
end

%Load Healthy participant data
sublists_healthy  = {'SDT001','AAT002','CAS003','MUS004','OET005','EAW006', ...
            'NIS007','CUD008','PEB009','EOS012','RNZ014','AMY015', ...
            'PBC016','CHY017','CRY018','TAA019','CAS020','EUS021','GYW022', ...
            'ROR024','LAS026','JUT028','DIT029','AOR030','SOS031', ...
            'LAS032','SRS033','OOS034','ELR035','SES036','DUS037','MRV038', ...
            'LAS039','SAW040','MIR041','CEX043','SOS044','REZ045','LAS046'};

subnames_healthy=lower(sublists_healthy);
addpath('D:\Jorge\Data\Fan_healthy_data')
inputpath  = 'D:\Jorge\Data\Fan_healthy_data\';
sPE_ERP_healthy_temp=load('sPE_ERP.mat');
cPE_ERP_healthy_temp=load('cPE_ERP.mat');

for isub=1:length(cPE_ERP_healthy_temp.cPE_avg)

    temp_std_healthy=[];
    temp_std_healthy=load(strcat(inputpath,sublists_healthy{isub},'\',subnames_healthy{isub},'_std_sampleMean_healthy'));%load standard condition sample mean data
    STD_ERP_healthy{isub,1}=temp_std_healthy.std_ERP_healthy{isub,1};

    temp_dev_healthy=[];
    temp_dev_healthy=load(strcat(inputpath,sublists_healthy{isub},'\',subnames_healthy{isub},'_deviant_timelock'));%load deviant condition data
    DEV_ERP_healthy{isub,1}=temp_dev_healthy.DEV_ERP_healthy;

    sPE_ERP_healthy{isub,1}=sPE_ERP_healthy_temp.sPE_avg{1,isub};
    cPE_ERP_healthy{isub,1}=cPE_ERP_healthy_temp.cPE_avg{1,isub};
    pPE_ERP_healthy{isub,1} = ft_math(cfg, DEV_ERP_healthy{isub,1},STD_ERP_healthy{isub,1});%LC variable is deviant- standard 
end
%% Perform statistical analyses- Independent samples t-test
load('biosemi64_neighb.mat');
cfg = [];
cfg.latency = 'all';
cfg.channel='all';
cfg.avgoverchan='no';
cfg.avgovertime='no';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'nonparametric_individual';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;  
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.correcttail      = 'alpha';
cfg.numrandomization = 1000; 
design=[ones(1,length(DEV_all)) ones(1,length(sPE_ERP_healthy))*2];

cfg.design           = design;
cfg.ivar             = 1;

Tstat(1)=ft_timelockstatistics(cfg,sPE_ERP_chronic{:},sPE_ERP_healthy{:});
Tstat(2)=ft_timelockstatistics(cfg,cPE_ERP_chronic{:},cPE_ERP_healthy{:});%no significant difference
Tstat(3)=ft_timelockstatistics(cfg,pPE_ERP_chronic{:},pPE_ERP_healthy{:});
%%
condition={'sPE', 'cPE','pPE'};
clusterpath='D:\Jorge\Scripts\Processing\Cluster_stat_values';%we save the cluster probability values to this directory for later use if need be
comp='chronic healthy comparison';
for icon=1:length(condition)
    for ipval=1:5
        pval(ipval,1)=Tstat(icon).posclusters(ipval).prob; pval(ipval,2)=Tstat(icon).negclusters(ipval).prob;
        save(strcat(clusterpath,'\',condition{icon},'_',comp),'pval')
    end
end
tsamp=DEV_all{1,1}.time;
label=DEV_all{1,1}.label;
% Create nice masks around our time of interest. More efficient for
% scrolling
for icon=1:3
    Mask(icon).condition=condition(icon);
    Mask(icon).pos(1,:)=tsamp;
    Mask(icon).pos(2:65,:)=Tstat(icon).posclusterslabelmat==1;
    Mask(icon).neg(1,:)=tsamp;
    Mask(icon).neg(2:65,:)=Tstat(icon).negclusterslabelmat==1;

    [~,Mask(icon).time_pos]=find(Mask(icon).pos(2:65,:)==1);
    [~,Mask(icon).time_neg]=find(Mask(icon).neg(2:65,:)==1);
    Mask(icon).TOI_pos=unique(Mask(icon).time_pos);
    Mask(icon).TOI_neg=unique(Mask(icon).time_neg);

    Mask(icon).Mask_TOI_pos=Mask(icon).pos(:,Mask(icon).TOI_pos(1):Mask(icon).TOI_pos(end));
    Mask(icon).Mask_TOI_neg=Mask(icon).neg(:,Mask(icon).TOI_neg(1):Mask(icon).TOI_neg(end));
end
%% Compute the grand averaged ERP without keeping subjects
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'no';
GA_sPE_chronic= ft_timelockgrandaverage(cfg, sPE_ERP_chronic{:});
GA_cPE_chronic= ft_timelockgrandaverage(cfg, cPE_ERP_chronic{:});%chronic pain participants
GA_pPE_chronic= ft_timelockgrandaverage(cfg, pPE_ERP_chronic{:});

GA_sPE_healthy= ft_timelockgrandaverage(cfg, sPE_ERP_healthy{:});
GA_cPE_healthy= ft_timelockgrandaverage(cfg, cPE_ERP_healthy{:});
GA_pPE_healthy= ft_timelockgrandaverage(cfg, pPE_ERP_healthy{:});%healthy participants
%% Plot significant clusters for stimulus-driven PE

    cfg=[];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    Diff_ERP(1) = ft_math(cfg, GA_sPE_chronic,GA_sPE_healthy);
    Diff_ERP(2) = ft_math(cfg, GA_cPE_chronic,GA_cPE_healthy);
    Diff_ERP(3) = ft_math(cfg, GA_pPE_chronic,GA_pPE_healthy);

for icon=3
    % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
    pos_cluster_pvals = [Tstat(icon).posclusters(:).prob];
    
    % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
    % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
    % respectively
    % pos_clust = find(pos_cluster_pvals < 0.025);
    % pos       = ismember(Tstat_sPE.posclusterslabelmat, pos_clust);
    % 
    % % and now for the negative clusters...
    % neg_cluster_pvals = [Tstat_sPE.negclusters(:).prob];
    % neg_clust         = find(neg_cluster_pvals < 0.025);
    % neg               = ismember(Tstat_sPE.negclusterslabelmat, neg_clust);
    
    pos = Tstat(icon).posclusterslabelmat == 1; % or == 2, or 3, etc.
    neg = Tstat(icon).negclusterslabelmat == 1;
    
    timestep      = 0.025; % timestep between time windows for each subplot (in seconds)
    sampling_rate = 500; % Data has a temporal resolution of 300 Hz
    sample_count  = length(Tstat(icon).time);
    % number of temporal samples in the statistics object
    j = [-.1:timestep:1.85]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*sampling_rate:sample_count]; % temporal endpoints in M/EEG samples
    
    % First ensure the channels to have the same order in the average and in the statistical output.
    % This might not be the case, because ft_math might shuffle the order
    [i1,i2] = match_str(Diff_ERP(icon).label, label);
    
    n = 1;
    for k = 58:78
       subplot(4,5,n);
       cfg = [];
       cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
       cfg.zlim = [-1.5 1.5];
       % If a channel is in a to-be-plotted cluster, then
       % the element of pos_int with an index equal to that channel
       % number will be set to 1 (otherwise 0).
    
       % Next, check which channels are in the clusters over the
       % entire time interval of interest.
       pos_int = zeros(numel(Diff_ERP(icon).label),1);
       neg_int = zeros(numel(Diff_ERP(icon).label),1);
       pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
       neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
    
       cfg.highlight   = 'on';
       % Get the index of the to-be-highlighted channel
       cfg.highlightchannel = find(pos_int | neg_int);%mess around with this
       cfg.comment     = 'xlim';
       cfg.commentpos  = 'title';
       cfg.layout      = 'biosemi64.lay';
       cfg.interactive = 'no';
       cfg.figure      = 'gca'; % plots in the current axes, here in a subplot
       ft_topoplotER(cfg, Diff_ERP(icon));
       n = n +1;
    end
end

% sPE pos cluster time is 1.432 to 1.488. ROI=[5 12 13 20 32 33 39:41 46:51 56 57]-1
%sPE neg cluster time is 1.556 to 1.616. ROI = [2:7 10:12 14 15 34 39 49]-1

%pPE pos cluster time is 1.432 to 1.51. ROI=[5 6 11:14 20 33 39:41 46:51 57]-1
%pPE neg cluster time is 1.25 to 1.33. ROI =[5 6 11:14 33 37 39:41 45:51 57]-1
