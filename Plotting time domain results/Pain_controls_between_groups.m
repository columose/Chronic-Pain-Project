datapath='D:\Jorge\Data\Field_trip_pre_process\';
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\fieldtrip-20221126');
addpath('D:\Jorge\Scripts\Processing\functions');

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
inputpath  = 'D:\Jorge\Data\Fan_healthy_data\ERP data\';
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
% %% Apply hilbert transform
% for isub=1:25
%     cfg = [];
%     cfg.channel    = 'all';
%     cfg.bpfilter   = 'yes';
%     cfg.bpfreq     = [4 8];
%     cfg.hilbert    = 'yes';
%     cfg.keeptrials = 'no';
%     temp{isub,1} = ft_preprocessing(cfg,pPE_ERP_chronic{isub,1});
%     pain_hilbert_indiv(isub,:,:)=temp{isub,1}.avg;
% end
% 
% pain_hilbert_GA=squeeze(mean(pain_hilbert_indiv(:,48,:),1));% electrode Cz
% tsamp=-0.1:.002:1.848;
% plot(tsamp,pain_hilbert_GA')
% %%
% for isub=1:39
%     cfg = [];
%     cfg.channel    = 'all';
%     cfg.bpfilter   = 'yes';
%     cfg.bpfreq     = [4 8];
%     cfg.hilbert    = 'yes';
%     cfg.keeptrials = 'no';
%     temp{isub,1} = ft_preprocessing(cfg,pPE_ERP_healthy{isub,1});
%     controls_hilbert_indiv(isub,:,:)=temp{isub,1}.avg;
% end
% 
% controls_hilbert_GA=squeeze(mean(controls_hilbert_indiv(:,48,:),1));% electrode Cz
% 
% tsamp=-0.1:.002:1.848;
% plot(tsamp,pain_hilbert_GA')
% hold on
% plot(tsamp,controls_hilbert_GA')
% hold on
% legend('Pain','Controls')
%% Compute the grand averaged ERP without keeping subjects
variable={'sPE_healthy','sPE_chronic','cPE_healthy','cPE_chronic','pPE_healthy','pPE_chronic'};
All_EEG(1).GA=sPE_ERP_healthy;
All_EEG(2).GA=sPE_ERP_chronic;
All_EEG(3).GA=cPE_ERP_healthy;
All_EEG(4).GA=cPE_ERP_chronic;
All_EEG(5).GA=pPE_ERP_healthy;
All_EEG(6).GA=pPE_ERP_chronic;
for ivar=1:6
    EEG_GA(ivar).variable=variable{ivar};% Create a Core structure for our eeg data. This makes processing more efficient

    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    cfg.keepindividual = 'no';
    EEG_GA(ivar).GA=ft_timelockgrandaverage(cfg, All_EEG(ivar).GA{:});
end

for i=1:length(sPE_ERP_healthy)%grand averaged structure 
    EEG_GA(1).individual(i,:,:)=sPE_ERP_healthy{i,1}.avg;
    EEG_GA(3).individual(i,:,:)=cPE_ERP_healthy{i,1}.avg;
    EEG_GA(5).individual(i,:,:)=pPE_ERP_healthy{i,1}.avg;
end

for i=1:length(cPE_ERP_chronic)
    EEG_GA(2).individual(i,:,:)=sPE_ERP_chronic{i,1}.avg;
    EEG_GA(4).individual(i,:,:)=cPE_ERP_chronic{i,1}.avg;
    EEG_GA(6).individual(i,:,:)=pPE_ERP_chronic{i,1}.avg;
end

%% Create a Core eeg structure for efficient looping
tsamp=OMI_all{1,1}.time;
label=OMI_all{1,1}.label;

EEG_Plot(1).Var='sPE positive group difference';
EEG_Plot(2).Var='sPE negative group difference';
EEG_Plot(3).Var='pPE positive group difference';
EEG_Plot(4).Var='pPE negative group difference';
EEG_Plot(5).Var='cPE non-significant group difference';


EEG_Plot(1).ROI=[5 12 13 20 32 33 39:41 46:51 56 57]-1;
% EEG_Plot(1).ROI=find(ismember(label,{'Cz'}));%artificial
EEG_Plot(2).ROI=[2:7 10:12 14 15 34 39 49]-1;%provide region of interest for each significant cluster
EEG_Plot(3).ROI=[5 6 11:14 20 33 39:41 46:51 57]-1;
EEG_Plot(4).ROI=[5 6 11:14 33 37 39:41 45:51 57]-1;
EEG_Plot(5).ROI=find(ismember(label,{'P6','P8','P10','CP6','TP8','T8'}));

startpoint=526;
endpoint=975;
timerange=startpoint:endpoint;%for plotting purposes. Index of the beginning of our final tone
zeropoint=1.2;
time=tsamp(timerange)-zeropoint;


EEG_Plot(1).TOI(1,1)=find(tsamp==1.432);EEG_Plot(1).TOI(1,2)=795;
EEG_Plot(2).TOI(1,1)=find(tsamp==1.556);EEG_Plot(2).TOI(1,2)=find(tsamp==1.616);%provide time of interest for each significant cluster
EEG_Plot(3).TOI(1,1)=find(tsamp==1.432);EEG_Plot(3).TOI(1,2)=find(tsamp==1.51);
EEG_Plot(4).TOI(1,1)=find(tsamp==1.25);EEG_Plot(4).TOI(1,2)=find(tsamp==1.33);
EEG_Plot(5).TOI(1,1)=find(tsamp==1.73);EEG_Plot(5).TOI(1,2)=find(tsamp==1.84);
%% Create core clusters for plotting
% Deviant-standard
outfolder={'sPE','sPE','pPE','pPE','cPE'};
COMP=[1,2;1,2;5,6;5,6;3,4];
for icon=1:5 %for all conditions. Efficient way of calculating all data
    EEG_Plot(icon).Topo_healthy=mean(EEG_GA(COMP(icon,1)).GA.avg(:,EEG_Plot(icon).TOI(1,1):EEG_Plot(icon).TOI(1,2)),2);%Create scalpmap and ERP data for each cluster
    EEG_Plot(icon).Topo_chronic=mean(EEG_GA(COMP(icon,2)).GA.avg(:,EEG_Plot(icon).TOI(1,1):EEG_Plot(icon).TOI(1,2)),2);
    EEG_Plot(icon).Difference=EEG_Plot(icon).Topo_chronic-EEG_Plot(icon).Topo_healthy;
    EEG_Plot(icon).mean_erp_healthy=mean(EEG_GA(COMP(icon,1)).GA.avg(EEG_Plot(icon).ROI,:),1);
    EEG_Plot(icon).mean_erp_chronic=mean(EEG_GA(COMP(icon,2)).GA.avg(EEG_Plot(icon).ROI,:),1);
    EEG_Plot(icon).outfolder=outfolder(icon);
end

EEG_Plot(1).minmax=[-2.4 2.4];
EEG_Plot(1).minmax_difference=[-2.4 2.4];

EEG_Plot(2).minmax=[-1.6 1.6];
EEG_Plot(2).minmax_difference=[-1.3 1.3];

EEG_Plot(3).minmax=[-4.2 4.2];
EEG_Plot(3).minmax_difference=[-2.5 2.5];

EEG_Plot(4).minmax=[-1.6 1.6];
EEG_Plot(4).minmax_difference=[-2.2 2.2];

EEG_Plot(5).minmax=[-1.1 1.1];
EEG_Plot(5).minmax_difference=[-0.7 0.7];

gray=[.5 .5 .5];
black=[0 0 0];
color=[black; black; black; black; gray];
% color=[black; black; black; black; black];
%% Create subplots for scalpmaps and ERP difference
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\eeglab2021.1_old');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\Variables');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Data Analysis\1 Process code');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\functions')
eeglab
close all
customcolormap = customcolormap_preset('red-white-blue');
load('ALLEEG.mat')
%% 
% outputpath='D:\Jorge\Figures\PE';
% scalppath='D:\Jorge\Figures\Time_frequency\Phase_single_plots';
outputpath='D:\Jorge\Figures\PE\Local_global';
for icon=1:length(EEG_Plot)
    fh(icon)=figure('Position',[200 500 800 250]);
    figure(icon)
    t=tiledlayout(1,4);
    t.TileSpacing='none';
    % nexttile(1);topoplot(EEG_Plot(icon).Topo_chronic,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax, 'style', 'map', 'shading', 'interp','emarker',{'.','k',5,1}); colormap(customcolormap);
    % if icon==1
    %     title('Pain')
    % end
    % 
    % hold on
    % nexttile(2);topoplot(EEG_Plot(icon).Topo_healthy,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax, 'style', 'map', 'shading', 'interp','emarker',{'.','k',5,1}); colormap(customcolormap);
    % colorbar('Ticks',EEG_Plot(icon).minmax)
    % if icon==1
    %     title('Controls')
    % end
    % hold on

    ax=gca;
    ax.FontWeight='bold';
    % nexttile(4);topoplot(EEG_Plot(icon).Difference,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax_difference, 'style', 'map', 'shading', 'interp','emarker',{'.','k',0,1},'emarker2',{[EEG_Plot(icon).ROI],'o',color(icon,:),4,1}); colormap(customcolormap);
    nexttile(4);topoplot(EEG_Plot(icon).Difference,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax_difference, 'style', 'map', 'shading', 'interp','electrodes','off','emarker2',{[EEG_Plot(icon).ROI],'o',color(icon,:),4,1}); colormap(customcolormap);
    % colorbar('Ticks',EEG_Plot(icon).minmax_difference)
    if icon==1
        title('Difference')
    end
    hold on

    ax=gca;
    ax.FontWeight='bold';
    fontsize(16,'points')

    % exportgraphics(fh(icon),strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'.png'),'Resolution',1200)
    % exportgraphics(fh(icon),strcat(scalppath,'\',EEG_Plot(icon).Var,'.png'),'Resolution',1200)
end
close all
%% Plot ERP
for icon=1:length(EEG_Plot)
    fh(icon) = figure('Position',[200 400 400 150]);
    EEG_Plot(icon).TOI_ERP=EEG_Plot(icon).TOI-startpoint;

    plot(time,EEG_Plot(icon).mean_erp_chronic(timerange),'Color',[.65 0 .65],'LineWidth',1.25)
    hold on
    plot(time,EEG_Plot(icon).mean_erp_healthy(timerange),'Color',[0 .8 0],'LineWidth',1.25)
    hold on

    gray=[.7 .7 .7];
    xline(0,'--',Color=gray)
    yline(0,'--',Color=gray)
    ylimit=[-3 3];

    v = [time(EEG_Plot(icon).TOI_ERP(1)) -ylimit; time(EEG_Plot(icon).TOI_ERP(1)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) -ylimit];
    f = [1 2 3 4];
    if icon==5
        patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'EdgeColor','none','FaceAlpha',0.2)
    else
        patch('Faces',f,'Vertices',v,'FaceColor','#93C47D','EdgeColor','none','FaceAlpha',0.2)
    end

    xlim([time(1) time(end)]);
    ylim(ylimit)
    yticks(ylimit)

    if icon==1
            leg=legend('CP','HC','Fontsize',6,'Box','off','Location','southwest');
            leg.ItemTokenSize=[8,8];
    end
    
    if icon<3 || icon==5
        xlabel('Time from omission(s)','FontWeight','bold')
    else
        xlabel('Time from 5th stimulus(s)','FontWeight','bold')
    end

    hold on
    box off
    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=7;
    ylabel('Amplitude(\muV)','FontWeight','bold')
    cd('D:\Jorge\Figures\PE')
    % exportgraphics(fh(icon),strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_ERP.jpg'),'Resolution',1200)
end
% close all
%% Create a variable with one data point per participant squeezed over ROI
% time
COMP=[1,2;1,2;5,6;5,6;3,4];
for icon=1:length(EEG_Plot)
    Bar_data(icon).Var=EEG_Plot(icon).Var;
    Bar_data(icon).over_ROI_1=squeeze(mean(EEG_GA(COMP(icon,1)).individual(:,EEG_Plot(icon).ROI,:),2));
    Bar_data(icon).data_points_1=squeeze(mean(Bar_data(icon).over_ROI_1(:,EEG_Plot(icon).TOI),2));%for squeezing data over relevant dimensions
    Bar_data(icon).over_ROI_2=squeeze(mean(EEG_GA(COMP(icon,2)).individual(:,EEG_Plot(icon).ROI,:),2));
    Bar_data(icon).data_points_2=squeeze(mean(Bar_data(icon).over_ROI_2(:,EEG_Plot(icon).TOI),2));

    Bar_data(icon).SEM_1=std(Bar_data(icon).data_points_1)./sqrt(size(sublists_healthy,2));
    Bar_data(icon).SEM_2=std(Bar_data(icon).data_points_2)./sqrt(size(sublists,2));
    Bar_data(icon).err=[Bar_data(icon).SEM_2 Bar_data(icon).SEM_1]';%for standard error of the mean and data points
    Bar_data(icon).amp_data_1=mean(Bar_data(icon).data_points_1);
    Bar_data(icon).amp_data_2=mean(Bar_data(icon).data_points_2);
    Bar_data(icon).amp_data=[mean(Bar_data(icon).data_points_2) mean(Bar_data(icon).data_points_1)]';
end

for icon=1:5
    fdr_vars{icon,1}=EEG_Plot(icon).Var;
end

% These stats were obtained from SPSS ancova with age as covariate
pval=[0.002; 0.0001; 0.002; 0.003; .055];
Bonferroni=.05/5;
Bonferroni_p=pval<Bonferroni;

FDR_corrected=f_fdr_correct(fdr_vars,pval);


%% Create a barchart of the mean amplitude squeezed over significant channel-time points in chronic and healthy groups
outputpath='D:\Jorge\Figures\PE';
for icon=1:length(EEG_Plot)
    fh(icon) = figure('Position',[200 500 400 300]);
    hb=bar(Bar_data(icon).amp_data);
    hb.FaceColor="flat";
    hb.CData(1,:)=[.65 0 .65];
    hb.CData(2,:)=[0 .8 0];
    hold on
    for k=1:size(Bar_data(icon).amp_data,2)
        % get x positions per group
        xpos= hb(k).XData + hb(k).XOffset;
        % draw errorbar
        errorbar(xpos, Bar_data(icon).amp_data(:,k), Bar_data(icon).err(:,k), 'LineStyle', 'none', ... 
            'Color', 'k', 'LineWidth', 1);
    end
    hold on
    box off
    ylim([-1.5 1.5])
    yticks(ylim)
    % ylabel('Amplitude(\muV)')
    ax=gca;
    ax.FontSize=13;
    ax.YLim=ylim;
    ax.XTickLabel={'CP','HC'};
    ax.FontWeight='bold';

    if icon<5
        sigstar([1,2],pval(icon))
    end

    if icon==5
        title('Corrected for age')
    end

    % exportgraphics(fh(icon),strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_Barchart.jpg'),'Resolution',1200)
end
% close all
%% Look at the MMN:P300 ratio in deviant - standard comparison
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\3 Data Collection')
healthy_age=table2array(readtable('2 Personal data.xlsx','Sheet','Raw data','Range','G2:G40'));
chronic_age=[18	19	19	46	29	62	43	68	28	62	33	58	61	24	38	67	22	42	41	30	37	42	25	32 67];

% average over region and time of interest for MMN and P00 differences to obtain ratio values
for isub=1:length(sublists)
    MMN.chronic(isub,1)=squeeze(mean(pPE_ERP_chronic{isub,1}.avg(EEG_Plot(4).ROI,EEG_Plot(4).TOI(1):EEG_Plot(4).TOI(2)),'all'));
    P300.chronic(isub,1)=squeeze(mean(pPE_ERP_chronic{isub,1}.avg(EEG_Plot(3).ROI,EEG_Plot(3).TOI(1):EEG_Plot(3).TOI(2)),'all'));
end
ratio.chronic=MMN.chronic./P300.chronic;

%repeat the same process for healthy subjects
for isub=1:length(sublists_healthy)
    MMN.healthy(isub,1)=squeeze(mean(pPE_ERP_healthy{isub,1}.avg(EEG_Plot(4).ROI,EEG_Plot(4).TOI(1):EEG_Plot(4).TOI(2)),'all'));%average over region and time of interest for MMN difference
    P300.healthy(isub,1)=squeeze(mean(pPE_ERP_healthy{isub,1}.avg(EEG_Plot(3).ROI,EEG_Plot(3).TOI(1):EEG_Plot(3).TOI(2)),'all'));
end
ratio.healthy=MMN.healthy./P300.healthy;

% Correlate VAS score with MMN:P300 ratio
cd('D:\Jorge\Data')
VAS_pain=table2array(readtable("AGES chronic pain eeg.xlsx",'Range','k2:k26'));
[RHO,pval]=corr(VAS_pain,ratio.chronic,'type','Spearman');

[RHO,pval]=corr(MMN.chronic,P300.chronic,'type','Spearman')
[RHO,pval]=corr(MMN.healthy,P300.healthy,'type','Spearman')



%% Correlate BAI and BDI with mean ERP
cd('D:\Jorge\Data')
% BDI minimal 0-13 mild 14-19 moderate 20-28 severe 29-63
% BAI minimal 0-7 mild 8-15 moderate 16-25 severe 26-63

BAI.chronic=table2array(readtable("AGES chronic pain eeg.xlsx",'Range','G2:G26'));
BDI.chronic=table2array(readtable("AGES chronic pain eeg.xlsx",'Range','I2:I26'));
VAS_pain=table2array(readtable("AGES chronic pain eeg.xlsx",'Range','k2:k26'));


[RHO,pval]=corr(VAS_pain,Bar_data(3).data_points_2,'type','Spearman')

cd('D:\Anusa_fan\For fan')
BAI.healthy=table2array(readtable("2 Personal data.xlsx",'Range','X2:X40'));
BDI.healthy=table2array(readtable("2 Personal data.xlsx",'Range','W2:W40'));

%Independent sampels t-tests between BDI and BAI scores
[~,BAI.pval]=ttest2(BAI.chronic,BAI.healthy,'Alpha',0.05/2);
[~,BDI.pval]=ttest2(BDI.chronic,BDI.healthy,'Alpha',0.05/2);

BDI.avg_chronic=mean(BDI.chronic);
BAI.avg_chronic=mean(BAI.chronic);

BDI.avg_healthy=mean(BDI.healthy);
BAI.avg_healthy=mean(BAI.healthy);

% Correlate mean amp data with BAI and BDI scores
pval=[];
RHO=[];
for icon=1:length(Bar_data)
    [RHO(icon).healthy_BDI, pval(icon).healthy_BDI]=corr(BDI.healthy,Bar_data(icon).data_points_1,'type','Spearman');
    [RHO(icon).chronic_BDI, pval(icon).chronic_BDI]=corr(BDI.chronic,Bar_data(icon).data_points_2,'type','Spearman');
end

for icon=1:length(Bar_data)
    [RHO(icon).healthy_BAI, pval(icon).healthy_BAI]=corr(BAI.healthy,Bar_data(icon).data_points_1,'type','Spearman');
    [RHO(icon).chronic_BAI, pval(icon).chronic_BAI]=corr(BAI.chronic,Bar_data(icon).data_points_2,'type','Spearman');
end

% ttest on age
% ttest2(BAI.chronic,BAI.healthy,'Alpha',0.05/3);
%% Plot nexttile ERP figures and save to folders
outputpath='D:\Jorge\Figures\PE\Nexttile test\';
fh=figure('Position',[200 500 400 450]);
t=tiledlayout('vertical');
t.TileSpacing='compact';

for icon=[5 4 3]
    EEG_Plot(icon).TOI_ERP=EEG_Plot(icon).TOI-startpoint;
    nexttile

    plot(time,EEG_Plot(icon).mean_erp_chronic(timerange),'Color',[.65 0 .65],'LineWidth',1.25)
    hold on
    plot(time,EEG_Plot(icon).mean_erp_healthy(timerange),'Color',[0 .8 0],'LineWidth',1.25)
    hold on

    gray=[.7 .7 .7];
    xline(0,'--',Color=gray)
    yline(0,'--',Color=gray)
    ylimit=[-3 3];

    v = [time(EEG_Plot(icon).TOI_ERP(1)) -ylimit; time(EEG_Plot(icon).TOI_ERP(1)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) -ylimit];
    f = [1 2 3 4];
    if icon==5
        patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'EdgeColor','none','FaceAlpha',0.2)
    else
        patch('Faces',f,'Vertices',v,'FaceColor','#93C47D','EdgeColor','none','FaceAlpha',0.2)
    end

    xlim([time(1) time(end)]);
    ylim(ylimit)
    yticks(ylimit)


    leg=legend('CP','HC','Fontsize',6,'Box','off','Location','southwest');
    leg.ItemTokenSize=[8,8];

    
    if icon==3
        xlabel('Time relative to event onset(s)','FontWeight','bold')
    end

    hold on
    box off
    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=7;
    ylabel('Amplitude(\muV)','FontWeight','bold')
    cd('D:\Jorge\Figures\PE')
    if icon==3
        % exportgraphics(fh,strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_ERP.png'),'Resolution',1200)
    end
end
% close all