% This script is for performing time frequency decomposition using EEGLAB data structures

clear
cd('D:\Jorge\Data\1.Raw_data')

for igroup=1

        if igroup==1
            inputpath  = 'D:\Jorge\Data\1.Raw_data\';
            outputpath='D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Pain\';
        
            sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
            'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};
            subnames = lower(sublists);
        
        else
        
            inputpath='D:\Jorge\Data\Fan_healthy_data\All_paradigm_data\1 Preprocessing_eeglab\';
            outputpath='D:\Jorge\Data\Time_freq_decomposition_EEGLAB\TF_phase_reflection\Controls\';
        
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
        
        
        addpath('C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\MATLAB\eeglab2021.1_old')
        addpath('D:\Jorge\Scripts\Processing\Time_freq\EEGLAB')
        
        eeglab 
        close all
        
        %% 0. Determine conditions and clusters
        session = 'PE';
        timerange = 1:975;
        conditions = {'std','dev','omi','omistd'};
        clusters_PE = {{'Cz'},{'Cz'},{'Cz'},{'Cz'}};
        binlabels = {'B1(condition1)','B2(condition2)','B3(condition3)','B4(condition4)'}; 
        
        %% 1. PE session
        for isub=1
        
            % 1.1 Determine which subject
                display(['Import EEG data: ',sublists{isub}]);
                sublist = sublists{isub};
                subname = subnames{isub};
                
                %%%%%%%%% If the first-time processing %%%%%%%%%
                %%%%%%%%%%%%%%%%%% RUN  BELOW %%%%%%%%%%%%%%%%%%
                % 1.2 Load EEG set file
        
                if igroup==1
                    EEG=[];
                    EEG = pop_loadset(strcat(inputpath,char(sublist),'\',char(subname),'_',session,'_interp.set'));
                else
                    EEG=[];
                    EEG=pop_loadset(strcat(inputpath,char(sublist),'\',char(subname),'_',session,'_final.set'));
                end
        
        
                    for icon=1
        
                        % 1.1 Determine which condition and clusters
                        condition = conditions{icon};
                        clusters  = clusters_PE{icon};
                        binlabel = binlabels{icon};
        
                       display(['Start running the TF decomposition for the condition: ',condition]);
        
                     % 1.3 Create index and empty datasets 
                        %index
                        index = 0; 
                        %dataset
                        dataset = []; 
                        EEGOUT_PE=[];
        
                        % 1.4 Append epoches to each bin dataset
                        for n = 1:(length(EEG.epoch))
                            if strcmp(EEG.epoch(n).eventbinlabel,binlabel) % count trials number of this bin
                                index = index + 1;
                                dataset(:,:,index) = EEG.data(:,:,n);
                            end
                        end
        
                        extra_points=750;
                        reflex_dataset=[];
                        reflex_dataset=[dataset(:,extra_points:-1:1,:) dataset dataset(:,extra_points:-1:1,:)];
                        
                        % add 1500ms of reflected data
                        times_pad= -1600:2:1848+1500; 
                        
                        % 1.5 Create EEG data for the condition
                        EEGOUT_PE.data     = reflex_dataset;
                        EEGOUT_PE.pnts     = EEG.pnts+extra_points*2; % multiply by two because we have extra points at the beginning and end
                        EEGOUT_PE.times=times_pad;
                        EEGOUT_PE.trials    = size(dataset,3);
                        EEGOUT_PE.chanlocs = EEG.chanlocs;

                        % Cut reflected data after decomposition 
                        cut_idx(1)=extra_points+1;
                        cut_idx(2)=find(times_pad==1848);
        
                        % 1.6 Run Morlet Wavelet analysis
                        baseline_window = [-100 0];
                        [~, ~,as_complex] = f_tf_decomp_ERP_ST_reflection(EEGOUT_PE,clusters,baseline_window,cut_idx);

                
                        % save(strcat(outputpath,sublists{isub},'\',subnames{isub},'_',condition,'_TF_phase.mat'),'as_complex')
                        % as_complex=[];
                    end
                
        end
end