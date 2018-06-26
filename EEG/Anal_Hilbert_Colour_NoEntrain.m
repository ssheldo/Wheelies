

% Select electrodes
elec = exp.singletrialselecs; %O1, Pz, Cz, O2, FCz
n_chan = length(elec); % number of channels

% Alpha
hipass = 7; % Lower bound of filter
lopass = 14; % Upper bound of filter
    
    
for xx = 1:length(exp.participants)
    
    % Data matrix, time points by channels (excluding the EOGs and refs)
    Y = permute(ALLEEG(xx).data(exp.brainelecs,:,:),[2 1 3]); %permute so times x chans x trials
    Fs = ALLEEG(xx).srate; % Sampling Frequency in Hz
    times = ALLEEG(xx).times; % Time in ms of data relative to event 

    Y_slow = Y; % Y would be your data file with time points by channels
    n_chan = size(Y,2); % number of channels
    n_trial = size(Y,3); % number of trials
    
    %//////////////////////////////////////////////////////////////////////
    % Filter to get the frequency band data
    order = 3;
    bandpass = [(hipass*2)/Fs, (lopass*2)/Fs]; % Bandwidth of bandpass filter
    [Bbp,Abp] = butter(order,bandpass); % Generation of Xth order Butterworth highpass filter

    for t = 1:n_trial %loop through trials
        for c = 1:n_chan %loop through channels
            Y_slow(:,c,t) = filtfilt(Bbp,Abp,double(Y(:,c,t))); % Butterworth bandpass filtering of Y
        end
    end
    
    clear t c bandpass order Bbp Abp Fs
    
    %//////////////////////////////////////////////////////////////////////
    %//////////////////////////////////////////////////////////////////////
    % Get power and phase information

    % Take the instantaneous angle of the filtered time series in radians using 
    % a hilbert transform
    Y_phase{xx} = angle(hilbert(Y_slow)); %phase
    Y_power{xx} = abs(hilbert(Y)).^2;  %power
end
clear xx

% reorder so electrode x time x trials
hilb_pwr = cell(1,length(exp.participants)); %pre-allocate
hilb_phi = cell(1,length(exp.participants)); %pre-allocate
for xx = 1:length(exp.participants)
    hilb_pwr{1,xx} = permute(Y_power{xx},[2 1 3]);
    hilb_phi{1,xx} = permute(Y_phase{xx},[2 1 3]);
end


