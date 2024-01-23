function [tf_chan, as_chan, as_complex_final] = f_tf_decomp_ERP_ST(EEG,Channels,baseline_window,cut_idx)

% tf chan = total power
% as chan = power per trial
% as_complex_final = phase locking values for each participants. Values transformed from complex to real numbers

% Note that this function is adapated from Mike Cohen's ANTS series. More information can be found in his youtube videos

% convert baseline time into indices
[~,baseidx(1)] = min(abs(EEG.times-baseline_window(1)));
[~,baseidx(2)] = min(abs(EEG.times-baseline_window(2)));

% wavelet parameters
min_freq =  1;
max_freq = 44;
num_frex = 44;

% set frequency range
frex = logspace(log10(min_freq),log10(max_freq),num_frex); % (min, max, number of frequency)
time = -2:1/500:2;
half_wave = (length(time)-1)/2;

% set a few different wavelet widths (number of wavelet cycles)
range_cycles = [3 13];
nCycles =  logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex); 

% FFT parameters
nKern = length(time);
nData = EEG.pnts*EEG.trials;
nConv = nKern+nData-1;

% Create empty variables
tf_chan = [];
as_chan = [];
as_complex_chan = [];

% initialize output time-frequency data
tf = zeros(length(frex),EEG.pnts);

% FFT of data (doesn't change on frequency iteration)
for elec = 1 : length(Channels)
    disp(['Working in channel ' num2str(elec)]);
    
    [~,chan_index] = find(strcmp({EEG.chanlocs.labels},Channels{elec}));
    dataX = fft(reshape(EEG.data(chan_index,:,:), 1,[]),nConv); %for eeglab matrix as opposed to FT

    as_freq = [];
    for fi=1:length(frex)

        % create wavelet and get its FFT
        s = nCycles(fi)/(2*pi*frex(fi));
        cmw = exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*s^2));
        cmwX = fft(cmw,nConv);
        cmwX = cmwX./max(cmwX);

        % run convolution
        as = ifft(cmwX.*dataX,nConv);
        as = as(half_wave+1:end-half_wave); % trim edges, cut half the wavelet at the beginning and half at the end
        as = reshape(as,EEG.pnts,EEG.trials);
        as_freq(fi,:,:) = as;
        
        % put power data into big matrix
        tf(fi,:) = mean(abs(as).^2,2);
        as_final(fi,:,:) = abs(as).^2;
        
        % Baseline normalization (decibel)
        as_final(fi,:,:) = 10*log10( bsxfun(@rdivide, as_final(fi,:,:), mean(as_final(fi,baseidx(1):baseidx(2),:),2)));
        
    end
    
    % decibel conversion
    tf_db = 10*log10( bsxfun(@rdivide, tf, mean(tf(:,baseidx(1):baseidx(2)),2)));
    
    tf_chan(elec,:,:) = tf_db(:,cut_idx(1):cut_idx(2));% Remove reflected data before saving
    as_chan(elec,:,:,:) = as_final;
    as_complex_chan(elec,:,:,:) = as_freq(:,cut_idx(1):cut_idx(2),:);
    as_complex_final=as_complex_chan;



% This final loop is only needed for the comparison of conditions with an unequal trial count

for iter=1:1000
     x=size(as_complex_chan,4);
     randx=randperm(x,35); % 35 in this case because every participant had a minimum of 35 trials per condition
     x=[];


     % Perform the conversion from complex to real numbers in the function to decrease computation time
     for ifreq=1:length(frex)
         temp.as(ifreq).data=squeeze(as_complex_chan(:,ifreq,:,randx));% chan x freq x time x trial

         temp.final(ifreq,:)=abs(mean(exp(1i*angle(temp.as(ifreq).data)),2)); % time x trial matrix
     end
     temp.reps(iter,:,:)=temp.final;
end

    temp.final=squeeze(mean(temp.reps,1));
    as_complex_final=temp.final;
    temp=[];

end
end

