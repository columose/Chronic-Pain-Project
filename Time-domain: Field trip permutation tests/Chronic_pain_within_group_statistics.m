% This script is for grand averaging the data, and determining regions of
% interets within the chronic pain group
datapath='D:\Jorge\Data\Field_trip_pre_process\';
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\fieldtrip-20221126');

sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};
subnames = lower(sublists);

conditions = {'std','dev','omi','omistd'};

%Load single trial data for condition
for i=1:length(sublists)
    cd(strcat(datapath,sublists{i}));
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
end

%% Check singleplot data. Change cell values for different participants
cfg=[];
cfg.channel='all';
cfg.layout = 'biosemi64.lay';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.baseline=[-0.1 -0.05];
ft_multiplotER(cfg, OMI_all{25,1})

%% Perform statistical analyses- Dependent samples t-test
load('biosemi64_neighb.mat');
cfg = [];
cfg.latency = 'all';
cfg.channel='all';
cfg.avgoverchan='no';
cfg.avgovertime='no';
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
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

design(1,:)=[ones(1,length(sublists)),(ones(1,length(sublists))*2)];
design(2,:) =[1:length(sublists),1:length(sublists)];

cfg.design           = design;
cfg.ivar             = 1;

%Note we always do deviant-standard in the comparisons
Tstat(1)=ft_timelockstatistics(cfg,OMI_all{:},STD_all{:});
Tstat(2)=ft_timelockstatistics(cfg,OMI_all{:},OMI_STD_all{:});%cPE is not significant
Tstat(3)=ft_timelockstatistics(cfg,DEV_all{:},STD_all{:});
%%
condition={'sPE', 'cPE','pPE'};
clusterpath='D:\Jorge\Scripts\Processing\Cluster_stat_values';%we save the cluster probability values to this directory for later use if need be
comp='chronic within';
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
%% Grand average the data
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'no';
GA_Std= ft_timelockgrandaverage(cfg, STD_all{:});
GA_Dev= ft_timelockgrandaverage(cfg, DEV_all{:});
GA_Omi= ft_timelockgrandaverage(cfg, OMI_all{:});
GA_Omi_std= ft_timelockgrandaverage(cfg, OMI_STD_all{:});
%% Plot significant clusters for various conditions
cfg=[];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
Diff_ERP(1) = ft_math(cfg, GA_Omi,GA_Std);
Diff_ERP(2) = ft_math(cfg, GA_Omi,GA_Omi_std);
Diff_ERP(3) = ft_math(cfg, GA_Dev,GA_Std);

for icon=3
    Diff_ERP(icon).condition=condition(icon);

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
figure(icon)
for k = 58:78
   subplot(5,4,n);
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
% Time is 1.27 - 1.37 for sPE negcluser 
%ROI = [12 47:51]-1. Neg and pos are the same time frame

%cPE is not significant


%Time is 1.464 to 1.558 for neglcuster pPE
% ROI = [4:6 11:14 20 33 37 39:41 46:51 57]-1


%% Ensure that omission and omission standard data are different
plot(tsamp,mean(GA_Omi.avg(46:51,:),1))
hold on
plot(tsamp,mean(GA_Omi_std.avg(46:51,:),1))
legend('Omi','Omi Std')


