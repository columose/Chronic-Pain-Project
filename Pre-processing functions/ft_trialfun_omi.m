function [trl, event] = ft_trialfun_omi(cfg)
    
    % read the header information and the events from the data
    hdr   = ft_read_header(cfg.dataset);
    event = ft_read_event(cfg.dataset);
    
    % search for "trigger" events
    binlabel  = {event(find(strcmp('trigger', {event.type}))).binlabel}';
    sample = [event(find(strcmp('trigger', {event.type}))).sample]';
    
    % determine the number of samples before and after the trigger
    pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
    posttrig =  round(cfg.trialdef.post * hdr.Fs);
    
    
    % look for triggers
    trl = [];
        
    for j = 1:(length(binlabel)-1)
      trg = binlabel(j);
      if strcmp(trg, 'B3(condition3)') == 1 
        trlbegin = sample(j) + pretrig;       
        trlend   = sample(j) + posttrig;       
        offset   = pretrig;
        newtrl   = [trlbegin trlend offset];
        trl      = [trl; newtrl];
      end
    end
end