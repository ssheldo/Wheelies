function Load_Colour_NoEntrain()


% load processing settings
load('byTargets_v4_Settings.mat');
% load('byWheel_v2_Settings.mat');


anal.tf = 'on'; % if loading TF data
anal.singletrials = 'on'; % if loading single trial data
anal.segments ='on'; % if loading epochs
anal.tfelecs = exp.brainelecs; %electrodes
anal.singletrialselecs = exp.singletrialselecs; %single trial electrodes


nparts = length(exp.participants);
nsets = length(exp.setname);


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
        
%         nevents = length(exp.events(i_set,:));
        nevents = 1;
        
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

                for ii = 1:length(exp.singletrialselecs)
                    i_chan = exp.singletrialselecs(ii);
                    % all_ersp is (participant x electrode).trials(freq x time x trial)
                    try %Unfortunately, this load procedure can break sometimes in a non-reproducible way. So if an error happens here, we wait half a second and try again.
                        channeldata = load([exp.pathname '\' exp.settings '\SingleTrials\' exp.participants{i_part} '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
%                         all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        all_ersp(i_part,ii) = struct2cell(channeldata);
                    catch
                        WaitSecs(.5)
                        channeldata = load([exp.pathname '\' exp.settings '\SingleTrials\' exp.participants{i_part} '\' EEG.chanlocs(i_chan).labels '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'elec_all_ersp');
%                         all_ersp(i_part,i_chan) = struct2cell(channeldata);
                        all_ersp(i_part,ii) = struct2cell(channeldata);
                    end
                    clear channeldata

                end
            end
        end
    end
    toc
end

eeglab redraw
% #########################################################################






