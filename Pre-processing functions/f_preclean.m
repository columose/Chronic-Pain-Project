function [EEGOUT,badchan_index] = f_preclean(session,sublist,subname,outputpath)
    %% Determine which paradigm
    if strcmp(session,'PE')
        Binlist = 'C:\Users\SC10\OneDrive - Trinity College Dublin\Documents\Trinity Lab Members\Anusha_Fan\1 Auditory Illusion Project\4 Scripts\1 ProcessingScripts\Binlister_PE.txt';
        epoch = [-100 1850];
        winlength = 5;
    elseif strcmp(session,'Pr')
        Binlist = 'D:\1 Projects (offline)\1 Auditory Illusion Project\4 Scripts\1 ProcessingScripts\Binlister_Pr.txt';
        epoch = [-298 1300];
        winlength = 5;
    elseif strcmp(session,'ZT')
        Binlist = 'D:\1 Projects (offline)\1 Auditory Illusion Project\4 Scripts\1 ProcessingScripts\Binlister_ZT_single_trls.txt';
        epoch = [-600 4500];
        winlength = 2;
    else
        dis('No session matched. Please check the input session and restart again.');
        return
    end

    %% Preprocessing 

    % 1. Upload one EDF file
    clear ALLEEG;
    EEGOUT_load = pop_biosig();

    % 1.1 Merge data frame (selective)
    prompt_merge = 'Do you need to merge any other file (0 or any more files to add):';
    mg_num = input(prompt_merge);

    if mg_num == 0 % No additional files
        EEG = pop_select(EEGOUT_load,'nochannel',[65:79]);
    else  % Has additional files
        % Append the first file
        ALLEEG(1) = EEGOUT_load;
        for img = 2:(mg_num+1)
            ALLEEG(img) = pop_biosig(); % Add all file
        end
        % Merge together
        EEGOUT_merge = pop_mergeset(ALLEEG,[1:mg_num + 1]);
        EEG = pop_select(EEGOUT_merge,'nochannel',[65:79]);
    end
    
    
    % 2. Channel Location
    EEGOUT_chans = pop_chanedit(EEG,'load',{'Chanloc_BioSemi64.loc','filetype','loc'});
    
    % 3. Down sampling rate
    EEGOUT_resamp = pop_resample(EEGOUT_chans, 500);
    EEG = pop_saveset (EEGOUT_resamp, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc.set'), 'filepath', strcat(outputpath,char(sublist)));
    
    % 4. Re-reference
    EEGOUT_reRef = pop_reref(EEG,[]);
    EEG = pop_saveset (EEGOUT_reRef, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_reRef.set'), 'filepath', strcat(outputpath,char(sublist)));
    
    % 5. Filtering
    EEGOUT_filt = pop_eegfiltnew(EEG, 0.55, 44);
    EEG = pop_saveset (EEGOUT_filt, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_reRef_filter.set'), 'filepath', strcat(outputpath,char(sublist)));
    
    % 6. Plot EEG and check bad channels
    display({'Please check bad channels.'})

    eegplot(EEG.data, 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000, 'eloc_file', EEG.chanlocs, 'events',EEG.event,'command','string','winlength',winlength);
    
    uiwait;
    
    % 7. Input the bad channel
    prompt_badcha = 'Please input the bad channel you want to remove(e.g. F1 F4 Cz):';
    badcha = input(prompt_badcha, 's');
    badchan = strsplit(badcha);
    badchan_index = [];
    
    % 8. Remove bad channels
    if strcmp(badcha,'0')
        display(strcat('No bad channels are selected for subject:',sublist));
    else
        for c = 1:length(badchan)
            badchan_index(c) = find(strcmp(char(EEG.chanlocs.labels), badchan(c)));
        end
        % Redo all above procedures
        EEGOUT_resamp = EEG;
    
        % Remove bad channels
        EEGOUT_badchan = pop_select(EEG,'nochannel', badchan);
        EEG = pop_saveset (EEGOUT_badchan, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_badchan.set'), 'filepath', strcat(outputpath,char(sublist)));
    
        % Reference
        EEGOUT_reRef = pop_reref(EEG,[]);
        EEG = pop_saveset (EEGOUT_reRef, 'filename', ...
            strcat(char(subname),'_',session,'_chanloc_reRef.set'), 'filepath', strcat(outputpath,char(sublist)));
        
        % Filtering
        EEGOUT_filt = pop_eegfiltnew(EEG, 0.55, 44);
        EEG = pop_saveset (EEGOUT_filt, 'filename', ...
            strcat(char(subname),'_',session,'_chanloc_reRef_filter.set'), 'filepath', strcat(outputpath,char(sublist)));            
    end
    
    % 9. Create an event list
    EEGOUT_elist = pop_creabasiceventlist(EEG,'Eventlist',strcat(outputpath,char(sublist),'\',char(subname),'_',session,'_elist.txt'), ...
        'AlphanumericCleaning', 'on', 'BoundaryNumeric', {-99}, 'BoundaryString', {'boundary'});
    EEG = pop_saveset (EEGOUT_elist, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist.set'), 'filepath', strcat(outputpath,char(sublist)));
    
    % 10. Assign bins to events
    EEGOUT_bin = pop_binlister(EEG,'BDF', ...
                    Binlist, ...
                    'ImportEL',strcat(outputpath,char(sublist),'\',char(subname),'_',session,'_elist.txt'), ...
                    'ExportEL',strcat(outputpath,char(sublist),'\',char(subname),'_',session,'_ba.txt'), 'SendEL2', 'EEG');
    EEG = pop_saveset (EEGOUT_bin, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist_bin.set'), 'filepath', strcat(outputpath,char(sublist)));
    
    % 11. Extract bin-based epochs
    EEGOUT_epoch = pop_epochbin(EEG,epoch,'none');
    EEG = pop_saveset (EEGOUT_epoch, 'filename', ...
        strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch.set'), 'filepath', strcat(outputpath,char(sublist)));
    
    % 12. Plot EEG and reject artifact epoches

    display({'Please reject artifact epoches.'})

    eegplot(EEG.data, 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000, 'eloc_file', EEG.chanlocs, 'events',EEG.event,'command','string','winlength',winlength);
    
    uiwait;

    TMPREJ = evalin("base",'TMPREJ');
    [trialrej elecrej] = eegplot2trial(TMPREJ,EEG.pnts,EEG.trials);
    EEGOUT_clean = pop_rejepoch(EEG,trialrej,0);
    EEGOUT = pop_saveset (EEGOUT_clean, 'filename', ...
         strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean.set'), 'filepath', strcat(outputpath,char(sublist)));
    
%     % 12.1 (IF BAD CHANNELS) Interpolate bad channels
%     if isempty(badchan_index) == 0 
%          
%         EEGOUT_interp = interpol(EEG, EEGOUT_chans.chanlocs, 'spherical');
%         EEGOUT = pop_saveset (EEGOUT_interp, 'filename', ...
%             strcat(char(subname),'_',session,'_chanloc_reRef_filter_elist_bin_epoch_clean.set'), 'filepath', strcat(outputpath,char(sublist)));
%     end
end