% This script is for plotting the significant differences between deviant
% and standard in healthy participants 

%Load Healthy participant data
inputpath  = 'D:\Jorge\Data\Fan_healthy_data\ERP data\';
sublists_healthy  = {'SDT001','AAT002','CAS003','MUS004','OET005','EAW006', ...
            'NIS007','CUD008','PEB009','EOS012','RNZ014','AMY015', ...
            'PBC016','CHY017','CRY018','TAA019','CAS020','EUS021','GYW022', ...
            'ROR024','LAS026','JUT028','DIT029','AOR030','SOS031', ...
            'LAS032','SRS033','OOS034','ELR035','SES036','DUS037','MRV038', ...
            'LAS039','SAW040','MIR041','CEX043','SOS044','REZ045','LAS046'};

subnames_healthy=lower(sublists_healthy);
addpath('D:\Jorge\Data\Fan_healthy_data')
for isub=1:length(subnames_healthy)
    cd(strcat(inputpath,sublists_healthy{isub}));
    temp_std_healthy=[];
    temp_std_healthy=load(strcat(char(subnames_healthy(isub)),'_std_sampleMean_healthy'));%load standard condition sample mean data
    STD_ERP_healthy{isub,1}=temp_std_healthy.std_ERP_healthy{isub,1};

    temp_dev_healthy=[];
    temp_dev_healthy=load(strcat(inputpath,sublists_healthy{isub},'\',subnames_healthy{isub}));%load deviant condition data
    DEV_ERP_healthy{isub,1}=temp_dev_healthy.DEV_ERP_healthy;
end

OMI_ERP_healthy=load('ERP_omi.mat');
OMISTD_ERP_healthy=load('ERP_omistd.mat');
temp_OMI_ERP_healthy=OMI_ERP_healthy.ERP_omi_TOI;
temp_OMISTD_ERP_healthy=OMISTD_ERP_healthy.ERP_omistd_TOI;
clear OMI_ERP_healthy OMISTD_ERP_healthy

for isub = 1:length(temp_OMI_ERP_healthy)
    OMI_ERP_healthy{isub,1} = temp_OMI_ERP_healthy(isub);
    OMISTD_ERP_healthy{isub,1} = temp_OMISTD_ERP_healthy(isub);
end


%% Compute the grand averaged ERP without keeping subjects
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\fieldtrip-20221126')
conditions = {'std','dev','omi','omistd'};

All_EEG(1).GA=STD_ERP_healthy;
All_EEG(2).GA=DEV_ERP_healthy;
All_EEG(3).GA=OMI_ERP_healthy;
All_EEG(4).GA=OMISTD_ERP_healthy;

for icon=1:4
    EEG_GA(icon).condition=conditions(icon);
    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    cfg.keepindividual = 'no';
    EEG_GA(icon).GA_avg= ft_timelockgrandaverage(cfg, All_EEG(icon).GA{:});
end

for icon=1:4
    for isub=1:length(sublists_healthy)
        EEG_GA(icon).individual(isub,:,:)=All_EEG(icon).GA{isub,1}.avg(:,:);
    end
end
%% Create a Core EEG structure for efficient looping
tsamp=DEV_ERP_healthy{1,1}.time;
label=DEV_ERP_healthy{1,1}.label;

% EEG_Plot(1).Var='sPE healthy';EEG_Plot(1).condition_1='Omi';EEG_Plot(1).condition_2='Std';
% EEG_Plot(2).Var='cPE healthy';EEG_Plot(2).condition_1='Omi';EEG_Plot(2).condition_2='Omi Std';
% EEG_Plot(3).Var='pPE healthy early';EEG_Plot(3).condition_1='Dev';EEG_Plot(3).condition_2='Std';
% EEG_Plot(4).Var='pPE healthy mid';EEG_Plot(4).condition_1='Dev';EEG_Plot(4).condition_2='Std';
% EEG_Plot(5).Var='pPE healthy late';EEG_Plot(5).condition_1='Dev';EEG_Plot(5).condition_2='Std';

EEG_Plot(1).Var='sPE healthy';EEG_Plot(1).condition_1='Omi XY';EEG_Plot(1).condition_2='Loc Dev';
EEG_Plot(2).Var='cPE healthy';EEG_Plot(2).condition_1='Omi XY';EEG_Plot(2).condition_2='Omi Std';
EEG_Plot(3).Var='pPE healthy early';EEG_Plot(3).condition_1='Glob Dev';EEG_Plot(3).condition_2='Loc Dev';
EEG_Plot(4).Var='pPE healthy mid';EEG_Plot(4).condition_1='Glob Dev';EEG_Plot(4).condition_2='Loc Dev';
EEG_Plot(5).Var='pPE healthy late';EEG_Plot(5).condition_1='Glob Dev';EEG_Plot(5).condition_2='Loc Dev';

EEG_Plot(1).ROI=find(ismember(label,{'Fpz','Fp1','Fp2','AF3','AF4','AF7','AF8','Fz','F1','F2','F4','FC2','FC1','FCz'}));
EEG_Plot(2).ROI=find(ismember(label,{'P6','P8','P10','CP6','TP8','T8'}));
EEG_Plot(3).ROI=[4:7 11:14 19 20 32 33 37 39:41 46:51 56 57]-1;
EEG_Plot(4).ROI=[19:21 51 52 56:59];
EEG_Plot(5).ROI=[12:14 19 46 47 50:52];


startpoint=526;
endpoint=975;
timerange=startpoint:endpoint;%for plotting purposes. Index of the beginning of our final tone
zeropoint=1.2;
time=tsamp(timerange)-zeropoint;


EEG_Plot(1).TOI(1,1)=find(tsamp==1.55);EEG_Plot(1).TOI(1,2)=find(tsamp==1.75);
EEG_Plot(2).TOI(1,1)=find(tsamp==1.73);EEG_Plot(2).TOI(1,2)=find(tsamp==1.84);%provide time of interest for each significant cluster
EEG_Plot(3).TOI(1,1)=770;EEG_Plot(3).TOI(1,2)=817;
EEG_Plot(4).TOI(1,1)=find(tsamp==1.552);EEG_Plot(4).TOI(1,2)=865;
EEG_Plot(5).TOI(1,1)=870;EEG_Plot(5).TOI(1,2)=918;
%%
% Create scalpmap ERP variables for the icons.
% Topo 1 refers to the first condition and topo 2 refers to the second
% condition
COMP=[3,1;3,4;2,1;2,1;2,1];%with this variable we can loop through different condition comparisons
for icon=1:length(EEG_Plot)
    EEG_Plot(icon).Topo_1=mean(EEG_GA(COMP(icon,1)).GA_avg.avg(:,EEG_Plot(icon).TOI(1,1):EEG_Plot(icon).TOI(1,2)),2);%Create scalpmap and ERP data for each cluster
    EEG_Plot(icon).Topo_2=mean(EEG_GA(COMP(icon,2)).GA_avg.avg(:,EEG_Plot(icon).TOI(1,1):EEG_Plot(icon).TOI(1,2)),2);
    EEG_Plot(icon).Difference=EEG_Plot(icon).Topo_1-EEG_Plot(icon).Topo_2;
    EEG_Plot(icon).mean_erp_1=mean(EEG_GA(COMP(icon,1)).GA_avg.avg(EEG_Plot(icon).ROI,:),1);
    EEG_Plot(icon).mean_erp_2=mean(EEG_GA(COMP(icon,2)).GA_avg.avg(EEG_Plot(icon).ROI,:),1);

    % EEG_Plot(icon).minmax='minmax';
    % EEG_Plot(icon).minmax_difference='minmax';
end

EEG_Plot(1).minmax=[-0.6 0.6];
EEG_Plot(1).minmax_difference=[-1.1 1.1];

EEG_Plot(2).minmax=[-0.9 0.9];
EEG_Plot(2).minmax_difference=[-1.1 1.1];

EEG_Plot(3).minmax=[-2.3 2.3];
EEG_Plot(3).minmax_difference=[-4.1 4.1];

EEG_Plot(4).minmax=[-1.4 1.4];
EEG_Plot(4).minmax_difference=[-1.4 1.4];

EEG_Plot(5).minmax=[-1.7 1.7];
EEG_Plot(5).minmax_difference=[-1.6 1.6];

EEG_Plot(1).outfolder='sPE';
EEG_Plot(2).outfolder='cPE';
for icon=3:5
    EEG_Plot(icon).outfolder='pPE';
end
%% Create subplots for scalpmaps and ERP difference
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\eeglab2023.0');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\Variables');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Data Analysis\1 Process code');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\functions')
eeglab 
close all;
customcolormap = customcolormap_preset('red-white-blue');
load('ALLEEG.mat')
%%
% outputpath='D:\Jorge\Figures\PE';
outputpath='D:\Jorge\Figures\PE\Local_global';

for icon=2:3
    EEG_Plot(icon).TOI_ERP=EEG_Plot(icon).TOI-startpoint;
    fh(icon)=figure('Position',[200 500 800 250]);
    t=tiledlayout(1,4);
    t.TileSpacing='none';
    % nexttile(1);topoplot(EEG_Plot(icon).Topo_1,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax, 'style', 'map', 'shading', 'interp','emarker',{'.','k',5,1}); colormap(customcolormap);
    % title(EEG_Plot(icon).condition_1)
    % 
    % hold on
    % nexttile(2);topoplot(EEG_Plot(icon).Topo_2,ALLEEG.chanlocs,'maplimits',EEG_Plot(icon).minmax, 'style', 'map', 'shading', 'interp','emarker',{'.','k',5,1}); colormap(customcolormap);
    % colorbar('Ticks',EEG_Plot(icon).minmax)
    % title(EEG_Plot(icon).condition_2)

    hold on
    ax=gca;
    ax.FontWeight='bold';
    nexttile(4);topoplot(EEG_Plot(icon).Difference,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax_difference, 'style', 'map', 'shading', 'interp','electrodes','off','emarker2',{[EEG_Plot(icon).ROI],'o','k',4,1}); colormap(customcolormap);
    % colorbar('Ticks',EEG_Plot(icon).minmax_difference)


    hold on
    ax=gca;
    ax.FontWeight='bold';
    fontsize(16,'points')


    exportgraphics(fh(icon),strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'.png'),'Resolution',1200)
end
close all
%% Plot ERP figures and save to folder
% outputpath='D:\Jorge\Figures\PE';
for icon=1:3
    fh(icon) = figure('Position',[200 400 350 150]);

    EEG_Plot(icon).TOI_ERP=EEG_Plot(icon).TOI-startpoint;

    plot(time,EEG_Plot(icon).mean_erp_1(timerange),'r','LineWidth',1.25)
    hold on
    plot(time,EEG_Plot(icon).mean_erp_2(timerange),'b','LineWidth',1.25)
    hold on
    gray=[.7 .7 .7];
    xline(0,'--',Color=gray)
    yline(0,'--',Color=gray)
    ylimit=[-1.5 1.5];

    v = [time(EEG_Plot(icon).TOI_ERP(1)) -ylimit; time(EEG_Plot(icon).TOI_ERP(1)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) -ylimit];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor','#93C47D','EdgeColor','none','FaceAlpha',0.2)
    xlim([time(1) time(end)]);
    ylim(ylimit)
    yticks(ylimit)
    
    leg=legend(EEG_Plot(icon).condition_1,EEG_Plot(icon).condition_2,'Fontsize',6,'Box','off','Location','southwest');
    leg.ItemTokenSize=[8,8];

    if icon==3
        xlabel('Time from 5th stimulus(s)',FontWeight='bold')
    else
        xlabel('Time from omission(s)','FontWeight','bold')
    end

    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=7;
    ylabel('Amplitude(\muV)','FontWeight','bold')
    
    hold on
    box off

    
    % exportgraphics(fh(icon),strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_ERP.jpg'),'Resolution',1200)
end
% close all
%% Create a variable with one data point per participant squeezed over ROI
% time
for icon=1:5
    Bar_data(icon).Var_1=EEG_Plot(icon).condition_1;
    Bar_data(icon).Var_2=EEG_Plot(icon).condition_2;

    Bar_data(icon).over_ROI_1=squeeze(mean(EEG_GA(COMP(icon,1)).individual(:,EEG_Plot(icon).ROI,:),2));
    Bar_data(icon).data_points_1=squeeze(mean(Bar_data(icon).over_ROI_1(:,EEG_Plot(icon).TOI),2));%for squeezing data over relevant dimensions
    Bar_data(icon).over_ROI_2=squeeze(mean(EEG_GA(COMP(icon,2)).individual(:,EEG_Plot(icon).ROI,:),2));
    Bar_data(icon).data_points_2=squeeze(mean(Bar_data(icon).over_ROI_2(:,EEG_Plot(icon).TOI),2));

    Bar_data(icon).SEM_1=std(Bar_data(icon).data_points_1)./sqrt(size(sublists_healthy,2));
    Bar_data(icon).SEM_2=std(Bar_data(icon).data_points_2)./sqrt(size(sublists_healthy,2));
    Bar_data(icon).err=[Bar_data(icon).SEM_1 Bar_data(icon).SEM_2]';%for standard error of the mean and data points
    Bar_data(icon).amp_data_1=mean(Bar_data(icon).data_points_1);
    Bar_data(icon).amp_data_2=mean(Bar_data(icon).data_points_2);
    Bar_data(icon).amp_data=[mean(Bar_data(icon).data_points_1) mean(Bar_data(icon).data_points_2)]';
end

%Perform dependent samples ttest
for icon=1:5
    [t_h(icon),t_p(icon)]=ttest(Bar_data(icon).data_points_1,Bar_data(icon).data_points_2,'Alpha',0.05/8);%8 is total number of dependent samples ttests
end
%% Plot barchart of squeezed data
CNDs={'Std','Dev','Omi','Omi Std'};
for icon=1:5
    fh(icon) = figure('Position',[200 500 400 300]);
    hb=bar(Bar_data(icon).amp_data);
    hb.FaceColor="flat";
    hb.CData(1,:)=[1 0 0];
    hb.CData(2,:)=[0 0 1];
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
    ylim([-1 1])
    ylabel('Amplitude(\muV)')
    ax=gca;
    ax.FontWeight='bold';
    ax.XTickLabel=[CNDs(COMP(icon,1)) CNDs(COMP(icon,2))];
    fontsize(14,'points')
    if icon==1
        sigstar([1,2],0.01)
    else
        sigstar([1,2],0.001)
    end
    cd('D:\Jorge\Figures\PE')
    outputpath='D:\Jorge\Figures\PE';
    % exportgraphics(fh(icon),strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_Barchart.jpg'),'Resolution',1200)
end

%% Plot ERP figures and save to folders
outputpath='D:\Jorge\Figures\PE\Nexttile test\';
fh=figure('Position',[200 500 350 350]);
t=tiledlayout('vertical');
t.TileSpacing='compact';
for icon=2:3
    EEG_Plot(icon).TOI_ERP=EEG_Plot(icon).TOI-startpoint;% This new variable provides the start and end points of our significant time of interest, with respect to the final stimulus offset
    nexttile


    plot(time,EEG_Plot(icon).mean_erp_1(timerange),'r','LineWidth',1.25) % timerange provides indices for time samples after final tone offset
    hold on
    plot(time,EEG_Plot(icon).mean_erp_2(timerange),'b','LineWidth',1.25)
    hold on

    grey=[.5 .5 .5];
    xline(0,'--',Color=grey)
    yline(0,'--',Color=grey)
    ylimit=[-1.5 1.5];

    v = [time(EEG_Plot(icon).TOI_ERP(1)) -ylimit; time(EEG_Plot(icon).TOI_ERP(1)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) ylimit; time(EEG_Plot(icon).TOI_ERP(2)) -ylimit];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,'FaceColor','#93C47D','EdgeColor','none','FaceAlpha',0.2)% v and f just provide coordinates to plot the green shading of our figures

    xlim([time(1) time(end)]);
    ylim(ylimit)
    yticks(ylimit)
    
    leg=legend(EEG_Plot(icon).condition_1,EEG_Plot(icon).condition_2,'Fontsize',6,'Box','off','Location','southwest');
    leg.ItemTokenSize=[8,8];%reduces size of legend line

    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=7;

    box off

    if icon==3
        xlabel('Time relative to event onset(s)','FontWeight','bold')
        % exportgraphics(fh,strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_ERP.png'),'Resolution',1200)
    end
    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=7;

    
    % exportgraphics(fh,strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_ERP.png'),'Resolution',1200)
end
% close all

