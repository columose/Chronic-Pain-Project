function sample_mean = f_timelock_sampleMean(target,sample)
    % This function is used to create a distribution of sample mean (by the concept of Central Limit Theorem) 
    %
    % Target = Timelocked fieldtrip data with the less number of trials to target the number
    %          of trails per sampling.
    % Sample = Timelocked fieldtrip data with the more number of trials to
    %          generate the sample mean
    
    % 1. distribution of sampling mean.
    sample_mean_dataset = [];
    for n = 1:1000 %1000 times sampling
        display(['Distribution of sampling mean: loop ',num2str(n)])
        % Extract same number of trails as target dataset from sample dataset
        randnum_trial = randperm(size(sample.trial,1),size(target.trial,1));
        
        % Ealculate mean ERP
        cfg = [];
        cfg.trials     = randnum_trial;
        cfg.channel    = 'all';
        cfg.latency    = 'all';
        cfg.keeptrials = 'no';
        data = ft_timelockanalysis(cfg, sample);  
        sample_mean_dataset{n} = data;
    end

    % 2. extract mean of the distribution of sample mean per channel per timepoint
    display(['Start extract mean of the sample mean'])
    sample_mean = [];
    cfg = [];
    cfg.channel = 'all';
    cfg.latency = 'all';
    cfg.parameter = 'avg';
    sample_mean   = ft_timelockgrandaverage(cfg, sample_mean_dataset{:});
end