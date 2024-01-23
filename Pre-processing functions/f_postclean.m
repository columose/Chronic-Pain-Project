function EEGOUT = f_postclean(session,sublist,subname,interpath,outputpath)
    %% Determine which paradigm
    if strcmp(session,'PE')
        winlength = 5;
    elseif strcmp(session,'Pr')
        winlength = 5;
    elseif strcmp(session,'ZT')
        winlength = 2;
    else
        dis('No session matched. Please check the input session and restart again.');
        return
    end

    %% Manual epoch cleaning

    % 1. Load EEG file
    EEG = pop_loadset(strcat(interpath,char(sublist),'\',char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean_ICA_clean.set'));

    % 2. Manual cleaning
    eegplot(EEG.data, 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000, 'eloc_file', EEG.chanlocs, 'events',EEG.event,'command','string','winlength',winlength);
    
    uiwait;

    TMPREJ = evalin("base",'TMPREJ');
    [trialrej, elecrej] = eegplot2trial(TMPREJ,EEG.pnts,EEG.trials);%elecrej
    EEGOUT_artrej = pop_rejepoch(EEG,trialrej,0);

    % 3. Save data
    EEG = pop_saveset (EEGOUT_artrej, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean_ICA_clean_artRej.set'), 'filepath', strcat(interpath,char(sublist)));
    
    EEGOUT = pop_saveset (EEGOUT_artrej, 'filename', ...
        strcat(char(subname),'_',session,'_final.set'), 'filepath', strcat(interpath,char(sublist)));%outputpath
    
end