%% The purpose of this script is to perform analysis and plotting on ITPC data

% analytic signal values = freq x time  (only real phase-locking values, conversion from complex to real in the decomposition script
% Note that we only look at data from electrode Cz because different channels won't have the exact same phase angle synchrony
clc;close all; clear

%Select which group to analyse
prompt='Select group (CP/HC):';
group=input(prompt,'s');

if strcmp(group,'CP')
    datapath='D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Pain';
    cd(datapath)

    sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};
    subnames = lower(sublists);

else
    datapath='D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Controls';
    cd(datapath)


    sublists= {'SDT001','AAT002','CAS003','MUS004','OET005','EAW006', ...
            'NIS007','CUD008','PEB009','EOS012','RNZ014','AMY015', ...
            'PBC016','CHY017','CRY018','TAA019','CAS020','EUS021','GYW022', ...
            'ROR024','LAS026','JUT028','DIT029','AOR030','SOS031', ...
            'LAS032','SRS033','OOS034','ELR035','SES036','DUS037','MRV038', ...
            'LAS039','SAW040','MIR041','CEX043','SOS044','REZ045','LAS046'};

    subnames  = {'sdt_001','aat_002','cas_003','mus_004','oet_005','eaw_006',...
                'nis_007','cud_008','peb_009','eos_012','rnz_014','amy_015',...
                'pbc_016','chy_017','cry_018','taa_019','cas_020','eus_021','gyw_022', ...
                'ror_024','las_026','jut_028','dit_029','aor_030','sos_031', ...
                'las_032','srs_033','oos_034','elr_035','ses_036','dus_037','mrv_038', ...
                'las_039','saw_040','mir_041','cex_043','sos_044','rez_045','las_046'};
end

% Define conditions 
conditions = {'std','dev','omi','omistd'};
conditions_full_name={'Standard','Deviant','Omission','Omission Standard'};
PE_new_name={'sPE','cPE','sE'};
% Each participants has decomposed data for the 4 conditions
for isub=1:length(sublists)
        for icon=1:4
            cd(strcat(datapath,'\',sublists{isub}))
            data=[];
            temp=[];

            data=load(strcat(subnames{isub},'_',conditions{icon},'_TF_phase.mat'));
            temp=data.as_complex;

            ITPC(icon).tf(isub,:,:)=temp;

        end
end

% Perform grand-averaging
for icon=1:4
      ITPC_GA(icon).tf=squeeze(mean(ITPC(icon).tf,1));% subject x freq x time (average over subject dimension
      ITPC_GA(icon).individual=ITPC(icon).tf;
end

%% Define time-frequency parameters for plotting
freq = logspace(log10(1),log10(44),44);
new_freq=freq(9:end);
new_freq_idx=9:44;
time= -.100:.002:1.848;
start=526;
omi_onset=651;
new_times=time(start:end);
zeropoint=1.2;
timerange=time(start:end)-zeropoint;


% Define limit boundaries
for icon=1:4
        lower.ITPC(icon)=min(ITPC_GA(icon).tf(new_freq_idx,start:end),[],'all');%determine min and max percent changes
        upper.ITPC(icon)=max(ITPC_GA(icon).tf(new_freq_idx,start:end),[],'all');
end
%% Plot within-group IPTC results


xlab={'Time from 5th stimulus (s)','Time from 5th stimulus (s)','Time from omission (s)','Time from omission (s)'};

addpath('D:\Jorge')
outputpath='D:\Jorge\Figures\Time_frequency\Phase_reflected';

for icon=1:4
      
      fh(icon)=figure;
      tiledlayout(4,5)
      colormap(brewermap([],"-RdBu"))

      nexttile
      contourf(timerange,new_freq,ITPC_GA(icon).tf(new_freq_idx,start:end),40,'linecolor','none')
      % contourf(time,new_freq,ITPC_GA(icon).tf(new_freq_idx,:),40,'linecolor','none')
      hold on
    
      climdb=[];
      climdb  = round([lower.ITPC(icon) upper.ITPC(icon)],2); %minmax is obtained from matrix

      cbar = colorbar; % auto colors
      % cbar.Label.String = 'ITPC';
      cbar.Label.Rotation = 270;
      cbar.Label.FontWeight = 'b';
      cbar.Limits=climdb;
      cbar.Ticks=climdb;
      set(gcf,'color','w');

      xline(0,'--k')

      box off
      ax=gca;
      ax.FontWeight='bold';
      fh(icon).WindowState = 'maximized';
      fontsize(12,"points");%for single images

      if icon==1
        title(group)
      elseif icon==4
          xlabel('Time relative to event onset(s)')
    
      end
      % xlabel(xlab{icon})
      if strcmp(group,'CP')
        ylabel('Frequency (Hz)')
      end
       
          % exportgraphics(fh(icon),strcat(outputpath,'\',group,'\',group,'_',conditions{icon},'.jpg'),'Resolution',1200)
          % close all
end
%% Load the ITPC results for both groups
datapath='D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Pain';
cd(datapath)
temp=load('ITPC_GA_pain.mat');
ITPC_GA_pain=temp.ITPC_GA;

temp=[];
datapath='D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Controls';
cd(datapath)
temp=load('ITPC_GA_controls.mat');
ITPC_GA_controls=temp.ITPC_GA;
temp=[];
%% Perform the permutation test on single subject data
for icon=1:4
    real_diff_avg(icon,1).data=ITPC_GA_pain(icon).tf - ITPC_GA_controls(icon).tf;
end

permpath='D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Permutation results';

% Read more about this on page 243 of Mike Cohen's matlab book
% for icon=1:4
%     disp(strcat('working on condition:',conditions{icon}))
%     group_power_avg=[];
%     group_power_avg = cat(1,ITPC_GA_pain(icon).individual, ITPC_GA_controls(icon).individual);%concatenate the data from two conditions
% 
%     for iperm=1:1000
%         disp(strcat('permutation:',num2str(iperm)))
% 
%         fakeord=[];
%         fakeord=randperm(size(group_power_avg,1));
% 
%         fake_diff=[];
%         fake_diff=squeeze(mean(group_power_avg(fakeord(1:25),:,:),1)-mean(group_power_avg(fakeord(26:64),:,:),1));% Null hypothesis is that there should be no group difference
%         permdiff{1,icon}.diff(iperm,:,:)=fake_diff;
%     end
% 
%     save(strcat(permpath,'\',conditions{icon},'_permdiff_','between_groups.mat'),"permdiff")
%     permdiff=[];
% end
%% 2.4 Calcluate the significant results
% pval = .025;
% 
% for icon=1:4
%     load(strcat(permpath,'\',conditions{icon},'_permdiff_between_groups.mat'))
%     permdiff{1,icon}.mean=squeeze(mean(permdiff{1,icon}.diff,1));
%     permdiff{1,icon}.std=squeeze(std(permdiff{1,icon}.diff,1));
%     permdiff{1,icon}.zmap=(real_diff_avg(icon).data-permdiff{1,icon}.mean)./permdiff{1,icon}.std; %Z-score evaluation of difference
% 
%     zthresh=permdiff{1,icon}.zmap;
%     zthresh(abs(zthresh)<norminv(1-pval) & abs(zthresh)>norminv(pval)) = 0;% Change the mask such that 
% 
%     Mask_zthresh=logical(zthresh);
%     save(strcat(permpath,'\',conditions{icon},'_perm_results.mat'),"Mask_zthresh")
%     permdiff=[];
% end
%% Create structure for indexing with the purpose of extracting significant data points
for icon=1:4
    Index(icon).conditions=conditions_full_name{icon};
end

Index(1).bands={'Theta','Delta'};
Index(2).bands={'Delta-Alpha'};
Index(3).bands={'Delta-pre','Delta-post'};
Index(4).bands={'Delta'};

Index(1).band_idx=[16 25; 9 13];
Index(2).band_idx=[9 30];%double check band lims here 
Index(3).band_idx=[11 14; 11 16];
Index(4).band_idx=[9 22];

Index(1).time=[1 226; 216 426];%time indices. Find function is uselesss
Index(2).time=[216 450];
Index(3).time=[1 226; 276 376];
Index(4).time=[1 240];

% Create a matrix with significant data points for all participants by applying mask to data
for icon=1:4
    
    Mask_zthresh=[];
    temp=[];
    load(strcat(permpath,'\',conditions{icon},'_perm_results.mat'));%load mask z thresh

     for isub=1:length(sublists)

         if strcmp(group,'Pain')
            temp(isub).data=squeeze(ITPC_GA_pain(icon).individual(isub,:,:)); 
         else
            temp(isub).data=squeeze(ITPC_GA_controls(icon).individual(isub,:,:));
         end

         temp(isub).data=temp(isub).data.*Mask_zthresh; %apply significance mask to single subject data
    
         Sig(isub).all_data=temp(isub).data(:,start:end);
    
         for iband=1:size(Index(icon).band_idx,1)

             Sig(isub).data{icon,iband}=Sig(isub).all_data(Index(icon).band_idx(iband,1):Index(icon).band_idx(iband,2),Index(icon).time(iband,1):Index(icon).time(iband,2));
             Sig(isub).mn_data{icon,iband}=mean(nonzeros(Sig(isub).data{icon,iband}),'all'); %mean over significant data points
                
             mn_data{icon}.band(isub,iband)=Sig(isub).mn_data{icon,iband}; % This data can now be analysed in SPSS with MANCOVA (age)
         end
     end
end

%Load age data for ancova
chronic_age=[18	19	19	46	29	62	43	68	28	62	33	58	61	24	38	67	22	42	41	30	37	42	25	32 67]';
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\3 Data Collection')
healthy_age=table2array(readtable('2 Personal data.xlsx','Sheet','Raw data','Range','G2:G40'));
%%
p_vals=[0.004; 0.17; 0.001; 0.227; 0.073; 0.001];% in DV order, stats obtained from SPSS MANCOVA(age)


n=1;
for icon=1:4
    for iband=1:size(Index(icon).bands,2)
        vars{n,1}=strcat(Index(icon).conditions,{' '},Index(icon).bands{iband});
        n=n+1;
    end
end

addpath('D:\Jorge\Scripts\Processing\functions')
FDR_corrected=f_fdr_correct(vars,p_vals);


%% Prepare bar chart data
temp=[];
cd('D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Controls')
temp=load('sig_data_controls.mat');
sig_data_controls=temp.mn_data;

temp=[];
cd('D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Pain')
temp=load('sig_data_pain.mat');
sig_data_pain=temp.mn_data;

temmp=[];
n=1;
for icon=1:4
    for iband=1:size(Index(icon).band_idx,1)

        Bar_data{n}.mn(1,:)=mean(sig_data_pain{1,icon}.band(:,iband),1);
        Bar_data{n}.SEM(1,:)=std(sig_data_pain{1,icon}.band(:,iband)./sqrt(25));

        Bar_data{n}.mn(2,:)=mean(sig_data_controls{1,icon}.band(:,iband),1);
        Bar_data{n}.SEM(2,:)=std(sig_data_controls{1,icon}.band(:,iband)./sqrt(39));

        Bar_data{n}.title=Index(icon).bands(iband);
        Bar_data{n}.cond=Index(icon).conditions;
        Bar_data{n}.pvals=p_vals(n,1);

        n=n+1;
    end
end
%% Create a barchart of the mean amplitude squeezed over signifcant data points in both groups
outputpath='D:\Jorge\Figures\Time_frequency\Phase_reflected\Barcharts\';
addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB')

for n=1:6
            fh(n)= figure('Position',[200 500 400 300]);
            hb=bar(Bar_data{1,n}.mn);
            hb.FaceColor="flat";
            hb.CData(1,:)=[.65 0 .65];
            hb.CData(2,:)=[0 .8 0];
            hold on
            for k=1
                % get x positions per group
                xpos= hb(k).XData + hb(k).XOffset;
                % draw errorbar
                errorbar(xpos, Bar_data{1,n}.mn, Bar_data{1,n}.SEM, 'LineStyle', 'none', ... 
                    'Color', 'k', 'LineWidth', 1);
            end
            hold on

            box off
            ylim([0 0.5])
            % ylabel('ITPC')
            ax=gca;
            ax.FontWeight='bold';
            ax.XTickLabel={'CP','HC'};
            fontsize(18,'points')
            title(Bar_data{1,n}.title)

            non_sig_conds=[2 4 5];

            if n~=non_sig_conds
                 sigstar([1,2],Bar_data{1,n}.pvals)
            end

            % exportgraphics(fh(n),strcat(outputpath,Bar_data{1,n}.cond,'_',char(Bar_data{1,n}.title),'_phase.png'),'Resolution',1200)
            close all
end


%% Determine difference figure limits for difference figure
for icon=1:4
    lim(1,icon)=min(real_diff_avg(icon).data(new_freq_idx,start:end),[],'all');
    lim(2,icon)=max(real_diff_avg(icon).data(new_freq_idx,start:end),[],'all');
    lim(3,:)=max(abs(lim)); % provides figure lims
end



%% Plot difference figures
figurepath='D:\Jorge\Figures\Time_frequency\Phase_reflected\Between groups';
addpath('D:\Jorge')
for icon=1:4

    fh=figure(icon);
    tiledlayout(4,5)

    colormap(brewermap([],"-RdBu"))
    load(strcat(permpath,'\',conditions{icon},'_perm_results.mat'));

    nexttile
    contourf(timerange,new_freq,real_diff_avg(icon,1).data(new_freq_idx,start:end),40,'linecolor','none');

    climdb=[];
    climdb  = round([-lim(3,icon) lim(3,icon)],2); %minmax is obtained from matrix
    set(gca,'clim',climdb,'ydir','normal','xlim',[timerange(1) timerange(end)])
    cbar = colorbar; % auto colors
    % cbar.Label.String = 'ITPC difference';
    cbar.Label.Rotation = 270;
    cbar.Label.FontWeight = 'b';
    cbar.Limits=climdb;
    cbar.Ticks=climdb;
    set(gcf,'color','w');

    if icon==1
        title('Group Difference')
    elseif icon==4
        xlabel('Time relative to event onset(s)')
    end

    % xlabel(xlab{icon})
    % ylabel('Frequency (Hz)')
    xline(0,'--k')

     box off
     ax=gca;
     ax.FontWeight='bold';
     ax.FontSize=12;
     fh.WindowState='maximized';
     hold on

     contour(timerange,new_freq,Mask_zthresh(new_freq_idx,start:end),1,'k');
     Mask_zthresh=[];

%     rectangles to shade significant differences after ancova. Sig differences obtained in the previous section
     if icon==1
        rectangle('Position',[-0.25 2 0.45 6.2],'LineStyle','--','EdgeColor','k','LineWidth',1.5)% bottom left corner, width, height
        rectangle('Position',[0.18 2 0.425 2],'LineStyle','--','EdgeColor','k','LineWidth',1.5)

     elseif icon==2
        rectangle('Position',[0.17 2 0.47 10],'LineStyle','--','EdgeColor','k','LineWidth',1.5)
     elseif icon==3
        rectangle('Position',[-0.225 2 0.35 2],'LineStyle','--','EdgeColor','k','LineWidth',1.5)% bottom left corner, width, height
        rectangle('Position',[0.3 2 0.25 2],'LineStyle','--','EdgeColor','k','LineWidth',1.5)
     else
         rectangle('Position',[-0.25 2 0.475 3],'LineStyle','--','EdgeColor','k','LineWidth',1.5)
     end

     % exportgraphics(fh,strcat(figurepath,'\',conditions{icon},'_between_groups.jpg'),'Resolution',1200)
     % close all
end

