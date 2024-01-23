% The purpose of this script is to perform z score normalisation, to rectify cortical maps and to compare the difference ERPs
%between Pain and Control group

clear;close all;clc;

binlabels = {'B1(condition1)','B2(condition2)','B3(condition3)','B4(condition4)'};

COMP=[3,1;3,4;2,1];
PE={'sPE','cPE','pPE'};

brainstorm

%Select which group to analyse
prompt='Select group (Pain/Controls):';
group=input(prompt,'s');

if strcmp(group,'Pain')
    inputpath  = 'D:\Jorge\Data\1.Raw_data\';
    resultpath='D:\Jorge\Data\Source_localisation_data\Pain\';

     sublists = {'EEG1','EEG2','EEG3','EEG4','EEG5','EEG6','EEG7','EEG8','EEG9','EEG10','EEG11','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20','EEG21','EEG22','EEG23','EEG24','EEG25'};
     subnames = lower(sublists);


    cd(inputpath)
    load("GA_input_pain.mat")
else
    inputpath  = 'D:\Jorge\Data\Fan_healthy_data\All_paradigm_data\1 Preprocessing_eeglab\';
    resultpath='D:\Jorge\Data\Source_localisation_data\Controls\';

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

    cd('D:\Jorge\Data\Fan_healthy_data\All_paradigm_data\1 Preprocessing_eeglab\');
    load("GA_input_controls.mat")
end

%% Perform grand-averaging on z score normalised difference ERPs

% Create link to z score file
GA_all=[GA_input_pain GA_input_controls];
temp=erase(GA_all,'.mat');
GA_all_z=strcat(temp,'_zscore.mat');

group_idx=[1,25;26,64];

% Grand - average the absolute values of the z scores

% for igroup=2 % for both pain and control groups
%     for icomp=1:3
%         % Process: Average: By trial group (grand average)
%         sFiles = bst_process('CallProcess', 'process_average', {GA_all_z{icomp,group_idx(igroup,1):group_idx(igroup,2)}}, [], ...
%         'avgtype',         7, ...  % By trial group (grand average)
%         'avg_func',        2, ...  % Average absolute values:  mean(abs(x))
%         'weighted',        0, ...
%         'scalenormalized', 0);
%     end
% end

%% Calculate difference of difference ERPs between groups (possibly not necessary?)

% Input files
sFiles = {...
    'Group_analysis/Within_group_differences_z_scores/results_average_230720_0957.mat', ...
    'Group_analysis/Within_group_differences_z_scores/results_average_230720_0958.mat', ...
    'Group_analysis/Within_group_differences_z_scores/results_average_230720_0961.mat'};
sFiles2 = {...
    'Group_analysis/Within_group_differences_z_scores/results_average_230720_0959.mat', ...
    'Group_analysis/Within_group_differences_z_scores/results_average_230720_0960.mat', ...
    'Group_analysis/Within_group_differences_z_scores/results_average_230720_1000.mat'};

% Process: Difference: A-B, abs
sFiles = bst_process('CallProcess', 'process_diff_ab', sFiles, sFiles2, ...
    'source_abs', 1);

% Process: Low-pass:40Hz
sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
    'highpass',    0, ...
    'lowpass',     40, ...
    'tranband',    0, ...
    'attenuation', 'strict', ...  % 60dB
    'ver',         '2019', ...  % 2019
    'mirror',      0, ...
    'overwrite',   1);


%% Perform non-parametric permutation t-test between the two groups

%This permutation test is at a vertex level
sFiles=[];
sFiles2=[];
TOI=[1.432 1.488; 1.730 1.840; 1.432 1.510; 1.250 1.330];% the last comparison is MMN in pain group
for i=1:64
    GA_all_z{4,i}=GA_all_z{3,i}; 
end

% for icomp=4
% 
%     % Input files
%     sFiles = {GA_all_z{icomp,1:25}};
%     sFiles2 = {GA_all_z{icomp,26:64}};
% 
%     % Process: Perm absmean test [1432ms,1488ms]          H0:(|mean(A)|=|mean(B)|), H1:(|mean(A)|<>|mean(B)|)
%     sFiles = bst_process('CallProcess', 'process_test_permutation2', sFiles, sFiles2, ...
%         'timewindow',     [TOI(icomp,1) TOI(icomp,2)], ... % if we use 4 then it's the second pPE comparioson for MMN
%         'scoutsel',       [], ...
%         'scoutfunc',      1, ...  % Mean
%         'isnorm',         0, ...
%         'avgtime',        1, ...
%         'iszerobad',      1, ...
%         'Comment',        '', ...
%         'test_type',      'absmean', ...  % Absolute mean test:    (works with unconstrained sources)T = (|mean(A)|-|mean(B)|) / sqrt(|var(A)|/nA + |var(B)|/nB)
%         'randomizations', 1000, ...
%         'tail',           'two');  % Two-tailed
% end
%% Extract one value from each participant averaged over TOI for each comparison
% To be used in MANCOVA

for icomp=1:4
    sFiles=[];
    sFiles={GA_all_z{icomp,:}};
    % Process: Extract values: [1432ms,1488ms] 68 scouts abs
    sFiles = bst_process('CallProcess', 'process_extract_values', sFiles, [], ...
    'timewindow', [TOI(icomp,1) TOI(icomp,2)],...
    'scoutsel',{'Destrieux', {'G_Ins_lg_and_S_cent_ins L', 'G_Ins_lg_and_S_cent_ins R', 'G_and_S_cingul-Ant L', 'G_and_S_cingul-Ant R', 'G_and_S_cingul-Mid-Ant L', 'G_and_S_cingul-Mid-Ant R', 'G_and_S_cingul-Mid-Post L', 'G_and_S_cingul-Mid-Post R', 'G_and_S_frontomargin L', 'G_and_S_frontomargin R', 'G_and_S_occipital_inf L', 'G_and_S_occipital_inf R', 'G_and_S_paracentral L', 'G_and_S_paracentral R', 'G_and_S_subcentral L', 'G_and_S_subcentral R', 'G_and_S_transv_frontopol L', 'G_and_S_transv_frontopol R', 'G_cingul-Post-dorsal L', 'G_cingul-Post-dorsal R', 'G_cingul-Post-ventral L', 'G_cingul-Post-ventral R', 'G_cuneus L', 'G_cuneus R', 'G_front_inf-Opercular L', 'G_front_inf-Opercular R', 'G_front_inf-Orbital L', 'G_front_inf-Orbital R', 'G_front_inf-Triangul L', 'G_front_inf-Triangul R', 'G_front_middle L', 'G_front_middle R', 'G_front_sup L', 'G_front_sup R', 'G_insular_short L', 'G_insular_short R', 'G_oc-temp_lat-fusifor L', 'G_oc-temp_lat-fusifor R', 'G_oc-temp_med-Lingual L', 'G_oc-temp_med-Lingual R', 'G_oc-temp_med-Parahip L', 'G_oc-temp_med-Parahip R', 'G_occipital_middle L', 'G_occipital_middle R', 'G_occipital_sup L', 'G_occipital_sup R', 'G_orbital L', 'G_orbital R', 'G_pariet_inf-Angular L', 'G_pariet_inf-Angular R', 'G_pariet_inf-Supramar L', 'G_pariet_inf-Supramar R', 'G_parietal_sup L', 'G_parietal_sup R', 'G_postcentral L', 'G_postcentral R', 'G_precentral L', 'G_precentral R', 'G_precuneus L', 'G_precuneus R', 'G_rectus L', 'G_rectus R', 'G_subcallosal L', 'G_subcallosal R', 'G_temp_sup-G_T_transv L', 'G_temp_sup-G_T_transv R', 'G_temp_sup-Lateral L', 'G_temp_sup-Lateral R', 'G_temp_sup-Plan_polar L', 'G_temp_sup-Plan_polar R', 'G_temp_sup-Plan_tempo L', 'G_temp_sup-Plan_tempo R', 'G_temporal_inf L', 'G_temporal_inf R', 'G_temporal_middle L', 'G_temporal_middle R', 'Lat_Fis-ant-Horizont L', 'Lat_Fis-ant-Horizont R', 'Lat_Fis-ant-Vertical L', 'Lat_Fis-ant-Vertical R', 'Lat_Fis-post L', 'Lat_Fis-post R', 'Pole_occipital L', 'Pole_occipital R', 'Pole_temporal L', 'Pole_temporal R', 'S_calcarine L', 'S_calcarine R', 'S_central L', 'S_central R', 'S_cingul-Marginalis L', 'S_cingul-Marginalis R', 'S_circular_insula_ant L', 'S_circular_insula_ant R', 'S_circular_insula_inf L', 'S_circular_insula_inf R', 'S_circular_insula_sup L', 'S_circular_insula_sup R', 'S_collat_transv_ant L', 'S_collat_transv_ant R', 'S_collat_transv_post L', 'S_collat_transv_post R', 'S_front_inf L', 'S_front_inf R', 'S_front_middle L', 'S_front_middle R', 'S_front_sup L', 'S_front_sup R', 'S_interm_prim-Jensen L', 'S_interm_prim-Jensen R', 'S_intrapariet_and_P_trans L', 'S_intrapariet_and_P_trans R', 'S_oc-temp_lat L', 'S_oc-temp_lat R', 'S_oc-temp_med_and_Lingual L', 'S_oc-temp_med_and_Lingual R', 'S_oc_middle_and_Lunatus L', 'S_oc_middle_and_Lunatus R', 'S_oc_sup_and_transversal L', 'S_oc_sup_and_transversal R', 'S_occipital_ant L', 'S_occipital_ant R', 'S_orbital-H_Shaped L', 'S_orbital-H_Shaped R', 'S_orbital_lateral L', 'S_orbital_lateral R', 'S_orbital_med-olfact L', 'S_orbital_med-olfact R', 'S_parieto_occipital L', 'S_parieto_occipital R', 'S_pericallosal L', 'S_pericallosal R', 'S_postcentral L', 'S_postcentral R', 'S_precentral-inf-part L', 'S_precentral-inf-part R', 'S_precentral-sup-part L', 'S_precentral-sup-part R', 'S_suborbital L', 'S_suborbital R', 'S_subparietal L', 'S_subparietal R', 'S_temporal_inf L', 'S_temporal_inf R', 'S_temporal_sup L', 'S_temporal_sup R', 'S_temporal_transverse L', 'S_temporal_transverse R'}}, ...
    'scoutfunc',  1, ...  % Mean
    'isnorm',     1, ...
    'avgtime',    1, ...
    'dim',        2, ...  % Concatenate time (dimension 2)
    'Comment',    '');

end
%%
% Manually import the TOI values from the brainstorm UI and export to excel and SPSS for MANCOVA with age as a covariate

fields= cellstr(sPE_TOI.Description);

for i=1:68
    fields{i}=strrep(fields{i},' ','');
end
Values=sPE_TOI.Value';
fields=fields';

%% Import the MANCOVA results and perform FDR correction
cd('D:\Jorge\Data\Source_localisation_data')
ROIS=table2array(readtable("sPE z score output.xlsx","Range",'B206:B274'));% for some reason we need to start on the previous column

% Determine which comparison is of interest
prompt='Select comparison (sPE/cPE/pPE/pPE_MMN/sPE_dest/cPE_dest/pPE_dest/pPE_MMN_dest):';
Comparison=input(prompt,'s');
if strcmp(Comparison,'sPE')
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","sPE","Range",'G207:G274'));
elseif strcmp(Comparison,'pPE')
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","pPE","Range",'G207:G274'));
elseif strcmp(Comparison,'cPE')
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","cPE","Range",'G207:G274'));
elseif strcmp(Comparison,'pPE_MMN')
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","pPE_MMN","Range",'G207:G274'));
elseif strcmp(Comparison,'sPE_dest')
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","sPE_dest","Range",'G447:G594'));
    cd('D:\Jorge\Data\Source_localisation_data')
    load("Dest_locs.mat")
    ROIS=[];
    ROIS=sPE_TOI.Description;
elseif strcmp(Comparison,'cPE_dest')
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","cPE_dest","Range",'G447:G594'));
    cd('D:\Jorge\Data\Source_localisation_data')
    load("Dest_locs.mat")
    ROIS=[];
    ROIS=sPE_TOI.Description;
elseif strcmp(Comparison,'pPE_dest')
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","pPE_dest","Range",'G447:G594'));
    fvals=table2array(readtable("sPE z score output.xlsx","Sheet","pPE_dest","Range",'F447:F594'));
    cd('D:\Jorge\Data\Source_localisation_data')
    load("Dest_locs.mat")
    ROIS=[];
    ROIS=sPE_TOI.Description;
else
    pvals=table2array(readtable("sPE z score output.xlsx","Sheet","pPE_MMN_dest","Range",'G447:G594'));
    fvals=table2array(readtable("sPE z score output.xlsx","Sheet","pPE_MMN_dest","Range",'F447:F594'));
    cd('D:\Jorge\Data\Source_localisation_data')
    load("Dest_locs.mat")
    ROIS=[];
    ROIS=sPE_TOI.Description;
    
end

addpath('D:\Jorge\Scripts\Processing\functions')
FDR_corrected=f_fdr_correct(ROIS,pvals,fvals);

for irow=1:148
    temp{irow,1}=erase(FDR_corrected{irow,1},' ');
end

%%
cd('D:\Jorge\Data')

VAS_pain=table2array(readtable("AGES chronic pain eeg.xlsx",'Range','k2:k26'));
VAS_pain=VAS_pain';
temp=z.Value;
[RHO,pval]=corr(VAS_pain,temp,'type','Spearman');