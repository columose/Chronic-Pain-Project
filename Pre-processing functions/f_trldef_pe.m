function cfg = f_trldef_pe(cfg,condition)
    cfg.trialdef.pre =0.1;
    cfg.trialdef.post = 1.848;
    
    if strcmp(condition,'std')
        cfg.trialfun = 'ft_trialfun_std';
    elseif strcmp(condition,'dev')
        cfg.trialfun = 'ft_trialfun_dev';
    elseif strcmp(condition,'omi')
        cfg.trialfun = 'ft_trialfun_omi';
    elseif strcmp(condition,'omistd')
        cfg.trialfun = 'ft_trialfun_omistd';
    end

    cfg = ft_definetrial(cfg);
    cfg.continuous = 'yes';
end
