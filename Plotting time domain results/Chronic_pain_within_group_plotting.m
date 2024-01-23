%This script is for plotting scalpmaps and ERP differences within the chronic pain group

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
%% Compute the grand averaged ERP with/without keeping subjects

All_EEG(1).GA=STD_all;
All_EEG(2).GA=DEV_all;
All_EEG(3).GA=OMI_all;
All_EEG(4).GA=OMI_STD_all;

for icon=1:4
    EEG_GA(icon).condition=conditions(icon);
    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    cfg.keepindividual = 'no';
    EEG_GA(icon).GA_avg= ft_timelockgrandaverage(cfg, All_EEG(icon).GA{:});%grand averaged ERP across participants
end

for icon=1:4
    for isub=1:length(sublists)
        EEG_GA(icon).individual(isub,:,:)=All_EEG(icon).GA{isub,1}.avg(:,:);%timelocked data per individual participant
    end
end
%% Create a Core EEG structure for efficient looping
tsamp=OMI_all{1,1}.time;
label=OMI_all{1,1}.label;

EEG_Plot(1).Var='sPE chronic negative';% We input data relevant for the sPE cluster in this structure
EEG_Plot(2).Var='cPE chronic non_significant';
EEG_Plot(3).Var='pPE chronic negative';

EEG_Plot(1).ROI= [12 47:51]-1;% -1 because the first row of the mask was for the time samples
EEG_Plot(2).ROI=find(ismember(label,{'P6','P8','P10','CP6','TP8','T8'}));
EEG_Plot(3).ROI=[4:6 11:14 20 33 37 39:41 46:51 57]-1;%provide region of interest for each significant cluster

startpoint=526;% index of zeropoint
endpoint=975;
timerange=startpoint:endpoint;%for plotting purposes. Indexes the time series from final tone offset
zeropoint=1.2;
time=tsamp(timerange)-zeropoint;


EEG_Plot(1).TOI(1,1)=find(tsamp==1.27);EEG_Plot(1).TOI(1,2)=find(tsamp==1.37);%provide time of interest for each significant cluster
EEG_Plot(2).TOI(1,1)=find(tsamp==1.73);EEG_Plot(2).TOI(1,2)=find(tsamp==1.84);
EEG_Plot(3).TOI(1,1)=783;EEG_Plot(3).TOI(1,2)=find(tsamp==1.558);

%Name conditions in each loop
% EEG_Plot(1).condition_1='Omi';EEG_Plot(1).condition_2='Std';
% EEG_Plot(2).condition_1='Omi';EEG_Plot(2).condition_2='Omi Std';
% EEG_Plot(3).condition_1='Dev';EEG_Plot(3).condition_2='Std';

EEG_Plot(1).condition_1='Omi XY';EEG_Plot(1).condition_2='Loc Dev';
EEG_Plot(2).condition_1='Omi XY';EEG_Plot(2).condition_2='Omi Std';
EEG_Plot(3).condition_1='Glob Dev';EEG_Plot(3).condition_2='Loc Dev';

%%
% Create scalpmap ERP variables for the conditions.
% Topo 1 refers to the first condition and topo 2 refers to the second condition. Must always be deviant-standard

COMP=[3,1;3,4;2,1];% This variable makes it easy to index different comparison conditions in each loop interation. COMP(1)= Omission, standard

for icon=1:length(EEG_Plot)
    EEG_Plot(icon).Topo_1=mean(EEG_GA(COMP(icon,1)).GA_avg.avg(:,EEG_Plot(icon).TOI(1,1):EEG_Plot(icon).TOI(1,2)),2);%Create scalpmap and ERP data for each cluster
    EEG_Plot(icon).Topo_2=mean(EEG_GA(COMP(icon,2)).GA_avg.avg(:,EEG_Plot(icon).TOI(1,1):EEG_Plot(icon).TOI(1,2)),2);
    EEG_Plot(icon).Difference=EEG_Plot(icon).Topo_1-EEG_Plot(icon).Topo_2;
    EEG_Plot(icon).mean_erp_1=mean(EEG_GA(COMP(icon,1)).GA_avg.avg(EEG_Plot(icon).ROI,:),1);
    EEG_Plot(icon).mean_erp_2=mean(EEG_GA(COMP(icon,2)).GA_avg.avg(EEG_Plot(icon).ROI,:),1);
    % 
    % EEG_Plot(icon).minmax='minmax';
    % EEG_Plot(icon).minmax_difference='minmax'; % I use these minmax variables before I've seen the topooplots themselves
end

EEG_Plot(1).minmax=[-0.8 0.8];
EEG_Plot(1).minmax_difference=[-1.2 1.2];

EEG_Plot(2).minmax=[-0.5 0.5];
EEG_Plot(2).minmax_difference=[-0.6 0.6];


EEG_Plot(3).minmax=[-1.2 1.2];
EEG_Plot(3).minmax_difference=[-1.7 1.7];

EEG_Plot(1).outfolder='sPE';% These folders are useful for organisation. Save figures to a PE folder with the following subfolders
EEG_Plot(2).outfolder='cPE';
EEG_Plot(3).outfolder='pPE';

gray=[.5 .5 .5];
black=[0 0 0];
color=[black; gray; black]; % Highlight channels in condition two with gray because they were non-significant
%% Load EEGLAB, colormap etc so we can plot figures
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\eeglab2023.0');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\Variables');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Data Analysis\1 Process code');
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\functions')
eeglab 
close all

customcolormap = customcolormap_preset('red-white-blue');
load('ALLEEG.mat')
%% Create subplots for scalpmaps and ERP difference
% Note, I use nexttile because it is preferred over subplot since Matlab 2019
% outputpath='D:\Jorge\Figures\PE';
outputpath='D:\Jorge\Figures\PE\Local_global';

for icon=1:length(EEG_Plot)
    fh(icon)=figure('Position',[200 500 800 250]); %sets the position of our figure. With these inputs, we get a rectangle with a short vertical and long horizontal
    t=tiledlayout(1,4);
    t.TileSpacing='none';% change input this to compress, and insert colorbar for nexttile 1 if you want to see its minmax values
    % nexttile(1);topoplot(EEG_Plot(icon).Topo_1,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax, 'style', 'map', 'shading', 'interp','emarker',{'.','k',5,1}); colormap(customcolormap);
    % title(EEG_Plot(icon).condition_1)
    % 
    % hold on
    % nexttile(2);topoplot(EEG_Plot(icon).Topo_2,ALLEEG.chanlocs,'maplimits',EEG_Plot(icon).minmax, 'style', 'map', 'shading', 'interp','emarker',{'.','k',5,1}); colormap(customcolormap);
    % colorbar('Ticks',EEG_Plot(icon).minmax)
    % title(EEG_Plot(icon).condition_2)
    % 
    % hold on
    ax=gca;
    ax.FontWeight='bold';
    nexttile(4);topoplot(EEG_Plot(icon).Difference,ALLEEG.chanlocs,'maplimits', EEG_Plot(icon).minmax_difference, 'style', 'map', 'shading', 'interp','electrodes','off','emarker2',{[EEG_Plot(icon).ROI],'o',color(icon,:),4,1}); colormap(customcolormap);
    % colorbar('Ticks',EEG_Plot(icon).minmax_difference)
    hold on

    ax=gca;
    ax.FontWeight='bold';
    fontsize(16,'points')
 
    exportgraphics(fh(icon),strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'.png'),'Resolution',1200)
end
 close all
%% Plot ERP figures and save to folders
for icon=2:3
    EEG_Plot(icon).TOI_ERP=EEG_Plot(icon).TOI-startpoint;% This new variable provides the start and end points of our significant time of interest, with respect to the final stimulus offset
    fh(icon) = figure('Position',[200 400 350 150]);


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
    if icon==2
        patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'EdgeColor','none','FaceAlpha',0.2)% I plot significant shading in green and non significant in white
    else
        patch('Faces',f,'Vertices',v,'FaceColor','#93C47D','EdgeColor','none','FaceAlpha',0.2)% v and f just provide coordinates to plot the green shading of our figures
    end
    xlim([time(1) time(end)]);
    ylim(ylimit)
    yticks(ylimit)
    
    leg=legend(EEG_Plot(icon).condition_1,EEG_Plot(icon).condition_2,'Fontsize',6,'Box','off','Location','southwest');
    leg.ItemTokenSize=[8,8];%reduces size of legend line

    if icon==3
        xlabel('Time relative to event onset(s)','FontWeight','bold')
    % else
    %     xlabel('Time from 5th stimulus(s)','FontWeight','bold')
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
for icon=1:3
    Bar_data(icon).Var_1=EEG_Plot(icon).condition_1;
    Bar_data(icon).Var_2=EEG_Plot(icon).condition_2;

    Bar_data(icon).over_ROI_1=squeeze(mean(EEG_GA(COMP(icon,1)).individual(:,EEG_Plot(icon).ROI,:),2));
    Bar_data(icon).data_points_1=squeeze(mean(Bar_data(icon).over_ROI_1(:,EEG_Plot(icon).TOI),2));%for squeezing data over relevant dimensions
    Bar_data(icon).over_ROI_2=squeeze(mean(EEG_GA(COMP(icon,2)).individual(:,EEG_Plot(icon).ROI,:),2));
    Bar_data(icon).data_points_2=squeeze(mean(Bar_data(icon).over_ROI_2(:,EEG_Plot(icon).TOI),2));

    Bar_data(icon).SEM_1=std(Bar_data(icon).data_points_1)./sqrt(size(sublists,2));%compute standard error from error bars
    Bar_data(icon).SEM_2=std(Bar_data(icon).data_points_2)./sqrt(size(sublists,2));
    Bar_data(icon).err=[Bar_data(icon).SEM_1 Bar_data(icon).SEM_2]';
    Bar_data(icon).amp_data_1=mean(Bar_data(icon).data_points_1);
    Bar_data(icon).amp_data_2=mean(Bar_data(icon).data_points_2);
    Bar_data(icon).amp_data=[mean(Bar_data(icon).data_points_1) mean(Bar_data(icon).data_points_2)]';
end

%Perform dependent samples ttest
for icon=1:3
    [t_h(icon),t_p(icon)]=ttest(Bar_data(icon).data_points_1,Bar_data(icon).data_points_2,'Alpha',0.05/8);% 8 is total number of dependent samples t tests
end


%% Plot barchart of squeezed data
CNDs={'Std','Dev','Omi','Omi Std'};
for icon=1:length(EEG_Plot)
    fh(icon) = figure('Position',[200 500 400 300]);
    hb=bar(Bar_data(icon).amp_data);
    hb.FaceColor="flat";
    hb.CData(1,:)=[1 0 0];% color of bar 1
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
    % ax.XTickLabel={EEG_Plot(icon).condition_1,EEG_Plot(icon).condition_2};
    ax.XTickLabel=[CNDs(COMP(icon,1)) CNDs(COMP(icon,2))];
    fontsize(14,'points')
    if icon==1
        sigstar([1,2],0.001)
    elseif icon==3
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
    if icon==2
        patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'EdgeColor','none','FaceAlpha',0.2)% I plot significant shading in green and non significant in white
    else
        patch('Faces',f,'Vertices',v,'FaceColor','#93C47D','EdgeColor','none','FaceAlpha',0.2)% v and f just provide coordinates to plot the green shading of our figures
    end
    xlim([time(1) time(end)]);
    ylim(ylimit)
    yticks(ylimit)
    
    leg=legend(EEG_Plot(icon).condition_1,EEG_Plot(icon).condition_2,'Fontsize',6,'Box','off','Location','southwest');
    leg.ItemTokenSize=[8,8];%reduces size of legend line

    if icon==3
        xlabel('Time relative to event onset(s)','FontWeight','bold')
    % else
    %     xlabel('Time from 5th stimulus(s)','FontWeight','bold')
    end
    ax=gca;
    ax.FontWeight='bold';
    ax.FontSize=7;
    ylabel('Amplitude(\muV)','FontWeight','bold')
    
    hold on
    box off
    % exportgraphics(fh,strcat(outputpath,'\',char(EEG_Plot(icon).outfolder),'\',EEG_Plot(icon).Var,'_ERP.png'),'Resolution',1200)
end
% close all