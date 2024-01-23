function f_binoperateERP(session,sublist,subname,inputpath,outputpath)
    % Determine the bin operation list by sessions
    if strcmp(session,'PE')
        bin_operation = {'b5 = b2 - b1 label = dev_std','b6 = b3 - b4 label = omi_dev_std'};
    elseif strcmp(session,'Pr')
        bin_operation = {'b13 = b11 - b12 label = YES_0_NO_0','b14 = b11 - b5 label = YES_0_YES75'};
    end

    % Load final EEG data
    EEG = pop_loadset(strcat(inputpath,char(sublist),'\',char(subname),'_',char(session),'_final.set'));
    
    % Bin operation
    ERPOUT = pop_averager(EEG, 'DSindex', 1, 'Criterion', 'good', 'ExcludeBoundary', 'on','SEM', 'on');
    if strcmp(session,'PE') || strcmp(session,'Pr')
        ERPOUT_bin = pop_binoperator(ERPOUT,bin_operation);
        ERP = pop_savemyerp(ERPOUT_bin,'erpname', char(subname), 'filename', strcat(char(subname),'_',char(session),'_ERP.erp'), 'filepath',strcat(outputpath,char(sublist)));
    else
        ERP = pop_savemyerp(ERPOUT,'erpname', char(subname), 'filename', strcat(char(subname),'_',char(session),'_ERP.erp'), 'filepath',strcat(outputpath,char(sublist)));
    end
    
end
