function Analysis_Colour_NoEntrain(exp)


anal.tf = 'on'; % if loading TF data
anal.singletrials = 'on'; % if loading single trial data
anal.segments ='on'; % if loading epochs
anal.tfelecs = exp.brainelecs; %electrodes
anal.singletrialselecs = exp.singletrialselecs; %single trial electrodes


nparts = length(exp.participants);
nsets = length(exp.setname);
nevents = length(exp.events(1,:));




[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% 
% targlatency = cell([nsets,nevents,nparts]);
% event_trials = cell([nsets,nevents,nparts]);
% all_event_phase = cell([nsets,nevents,nparts,anal.tfelecs(end)]);
% 

% #########################################################################
%% Load the data
%The main loop loops over events, then participants, then sets.
for i_set = 1:nsets
    exp.setname(i_set)
    tic
    
    for i_part = 1:nparts
        sprintf(['Loading Participant ' num2str(exp.participants{i_part}) '...' ])
        
        nevents = length(exp.events(i_set,:));
        
        for i_event = 1:nevents
            
            filename = [exp.participants{i_part} '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}];
            
            % Load the Time frequency data, if needed.
            if strcmp('on',anal.tf) == 1 % only load these variables if we are loading time-frequency data
                
                %The variable ersp will be a 6D variable: (participants x sets x events x electrodes x frequencies x timepoints).
                ersp(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.pathname '\' exp.settings '\TimeFrequency\' 'TF_' filename '.mat'],'ersp'));
                itc(i_part,i_set,i_event,:,:,:) = struct2array(load([exp.pathname '\' exp.settings '\TimeFrequency\' 'TF_' filename '.mat'],'itc'));
                
                if i_part == 1 && i_set == 1 && i_event == 1 %load time and freq data
                    times = struct2array(load([exp.pathname '\' exp.settings '\TimeFrequency\' 'TF_' filename '.mat'],'times'));
                    freqs = struct2array(load([exp.pathname '\' exp.settings '\TimeFrequency\' 'TF_' filename '.mat'],'freqs'));
                end
                
            end

            % Load the EEGLAB datasets, if needed.
            if strcmp('on',anal.segments) == 1 || strcmp('on',anal.singletrials) == 1 % only load these variables if we are loading either ERP or single trial data
                try
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                catch
                    WaitSecs(.5)
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                end
            elseif strcmp('on',anal.tf) == 1 % if we are loading time-frequency data only, then we just need one of these.
                if i_part == 1 && i_set == 1 && i_event == 1
                    EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                end
            end
            
            % Load the Single Trial complex values, if needed
            if strcmp('on',anal.singletrials) == 1 % only load these variables if we are loading single trial data
                
%                 ntrigs = length(exp.events{i_set});
                
                if i_part == 1 && i_set == 1 && i_event == 1
                    times = struct2array(load([exp.pathname '\' exp.settings '\TimeFrequency\' 'TF_' filename '.mat'],'times'));
                    freqs = struct2array(load([exp.pathname '\' exp.settings '\TimeFrequency\' 'TF_' filename '.mat'],'freqs'));
                end
                
                %This block finds the latency of each event listed in exp.events, and the trials it appeared on.
%                 for nperevent = 1:ntrigs
%                     for i_trial = 1:EEG.trials
%                         if any(strcmp(num2str(exp.events{i_set,i_event}(nperevent)),EEG.epoch(i_trial).eventtype)) == 1
%                             targlatency{i_set,i_event,i_part} = [targlatency{i_set,i_event,i_part} EEG.epoch(i_trial).eventlatency(find(strcmp(num2str(exp.events{i_set,i_event}(nperevent)),EEG.epoch(i_trial).eventtype)))];
%                             event_trials{i_set,i_event,i_part} = [event_trials{i_set,i_event,i_part} i_trial];
%                         end
%                     end
%                 end

                for i_chan = anal.singletrialselecs
                    % all_ersp is (participant x electrode).trials(freq x time x trial)
                    try %Unfortunately, this load procedure can break sometimes in a non-reproducible way. So if an error happens here, we wait half a second and try again.
                        channeldata = load([exp.pathname '\' exp.settings '\SingleTrials\' exp.participants{i_part} '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
%                         all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        all_ersp(i_part) = struct2cell(channeldata);
                    catch
                        WaitSecs(.5)
                        channeldata = load([exp.pathname '\' exp.settings '\SingleTrials\' exp.participants{i_part} '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
%                         all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        all_ersp(i_part) = struct2cell(channeldata);
                    end
                    clear channeldata
                    
%                    
                end
            end
        end
    end
    toc
end

% #########################################################################

    
    
freqlim = [7 14];
freqband = find(freqs>=freqlim(1) & freqs<=freqlim(2)); %this code finds the frequencies you want from the freqs variable

timelim = [-100 0];
timewin = find(times>=timelim(1) & times<=timelim(2)); %this code finds the times you want from the timess variable


i_part = 2;

% all_ersp is (participant x electrode).trials(freq x time x trial)
current_complex = [squeeze(all_ersp{i_part}(freqband,timewin,:))];
current_power = squeeze(mean(mean(abs(current_complex),1),2));

current_phase = squeeze(angle(current_complex));
for pha = 1:length(current_phase(1,1,:))
%     [circ_mean_v,range_v,X,Y,cos_a,sin_a] = circle_mean(circ_rad2ang(current_phase(:,:,1))); 
    circ_mean_v(pha) = circle_mean(circ_rad2ang(current_phase(:,:,pha))); 
end



[n,m] = size(ALLEEG(i_part).rejtrial);
pip = 1;
for ni = 1:n
    for mi = 1:m
        if ~isempty(ALLEEG(i_part).rejtrial(ni,mi).ids)
            rejlist{pip} = ALLEEG(i_part).rejtrial(ni,mi).ids;
            pip = 1 + pip;
        end
    end
end


tmplist = [rejlist{1:end}];

err_deg_tmp = ALLEEG(i_part).error_deg;
err_deg_tmp(tmplist) = [];
color_grp_tmp = ALLEEG(i_part).color_grp;
color_grp_tmp(tmplist) = [];


%standardize
% err_norm = normcdf(err_deg_tmp,0,62.73);
% err_norm_g = normcdf(error_deg_g,0,48.36);
% err_norm_b = normcdf(error_deg_b,0,88);
% err_norm_r = normcdf(error_deg_r,0,42.59);

% Split trials by target color
error_deg_g = err_deg_tmp(color_grp_tmp == 1);
error_deg_b = err_deg_tmp(color_grp_tmp == 2);
error_deg_r = err_deg_tmp(color_grp_tmp == 3);

current_power_g = current_power(color_grp_tmp == 1);
current_power_b = current_power(color_grp_tmp == 2);
current_power_r = current_power(color_grp_tmp == 3);

circ_mean_v_g = circ_mean_v(color_grp_tmp == 1);
circ_mean_v_b = circ_mean_v(color_grp_tmp == 2);
circ_mean_v_r = circ_mean_v(color_grp_tmp == 3);



figure; scatter(current_power, abs(err_deg_tmp))
figure; scatter(current_power_g, abs(error_deg_g))
figure; scatter(current_power_b, abs(error_deg_b))
figure; scatter(current_power_r, abs(error_deg_r))

figure; scatter(current_power, err_deg_tmp)
figure; scatter(current_power_g, error_deg_g)
figure; scatter(current_power_b, error_deg_b)
figure; scatter(current_power_r, error_deg_r)

[rho,pval] = corr(current_power, abs(err_deg_tmp)')
[rho,pval] = corr(current_power_g, abs(error_deg_g)')
[rho,pval] = corr(current_power_b, abs(error_deg_b)')
[rho,pval] = corr(current_power_r, abs(error_deg_r)')

[rho,pval] = corr(current_power, err_deg_tmp')
[rho,pval] = corr(current_power_g, error_deg_g')
[rho,pval] = corr(current_power_b, error_deg_b')
[rho,pval] = corr(current_power_r, error_deg_r')


figure; scatter(circ_mean_v, abs(err_deg_tmp))

[rho,pval] = corr(circ_mean_v', abs(err_deg_tmp)')
[rho,pval] = corr(circ_mean_v_g', abs(error_deg_g)')
[rho,pval] = corr(circ_mean_v_b', abs(error_deg_b)')
[rho,pval] = corr(circ_mean_v_r', abs(error_deg_r)')

[rho,pval] = corr(circ_mean_v', err_deg_tmp')
[rho,pval] = corr(circ_mean_v_g', error_deg_g')
[rho,pval] = corr(circ_mean_v_b', error_deg_b')
[rho,pval] = corr(circ_mean_v_r', error_deg_r')

clear circ_mean_v rejlist color_grp_tmp err_deg_tmp current_power


