
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% Variables working with:
% ersp(i_sub,i_cond,i_perm,i_chan,:,:)
% itc(i_sub,i_cond,i_perm,i_chan,:,:)
% powbase,times,freqs

% The variables ersp and itc will be a 6D variable: 
% (participants x conditions x events x electrodes x frequencies x timepoints)
% (participants x sets x events x electrodes x frequencies x timepoints)

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% period = 1/EEG.srate; 
% time (in s) = [EEG.event.latency]*period
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%%               BEH Data Corrected For Rejected Trials
% /////////////////////////////////////////////////////////////////////////
% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
targ_colorgrp = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % ---------------------
    [n,m] = size(ALLEEG(i_part).rejtrial);
    pip = 1;
    for ni = 1:n
        for mi = 1:m
            if ~isempty(ALLEEG(i_part).rejtrial(ni,mi).ids)
                rejlist{pip} = ALLEEG(i_part).rejtrial(ni,mi).ids;
                pip = 1 + pip;
            end
        end
        clear mi
    end
    if pip > 1 %if trials were rejected
        tmplist = [rejlist{1:end}];
        err_deg_tmp = ALLEEG(i_part).error_deg;
        err_deg_tmp(tmplist) = [];
        color_grp_tmp = ALLEEG(i_part).color_grp;
        color_grp_tmp(tmplist) = [];
    elseif pip == 1 %if no trials were rejected, rejlist variable not created
        err_deg_tmp = ALLEEG(i_part).error_deg;
        color_grp_tmp = ALLEEG(i_part).color_grp;
    end
    % create variable with selected BEH 
    resp_errdeg{i_part} = err_deg_tmp;
    targ_colorgrp{i_part} = color_grp_tmp;
    
    clear tmplist rejlist n m err_deg_tmp color_grp_tmp pip tmplist ni
end
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////


%finds the frequencies you want from the freqs variable
freqband = [7 14]; %alpha
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-260 -250];
timephi = find(times>=timewin(1) & times<=timewin(2));

current_phase = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
x_phase = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
n_phase = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
   
    % Get error +/-1 IQRs of the mean
    mean_errdeg = nanmedian(resp_errdeg{i_part});
    sd_errdeg = iqr(resp_errdeg{i_part})/2;
    uplim = mean_errdeg + (sd_errdeg);
    lowlim = mean_errdeg - (sd_errdeg);
    
    for i_elect = 1:length(exp.singletrialselecs)
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_phase{i_part,i_elect}(:,:,i_trial) = squeeze(angle(part_ersp(freqlim,timephi,i_trial)));
        end
        clear part_ersp i_trial

        for i_freq = 1:length(freqlim)
            % get phase at each frequency by high/low error trials
            x_phase{i_part,i_elect}(i_freq,:,:) = current_phase{i_part,i_elect}(i_freq,1,[find(resp_errdeg{i_part}>=uplim) find(resp_errdeg{i_part}<=lowlim)]);
            n_phase{i_part,i_elect}(i_freq,:,:) = current_phase{i_part,i_elect}(i_freq,1,resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim);
        end
        clear i_freq  
        
    end
end
clear i_elect i_part

        % plot trial power on a histogram
        i_freq = 3;
        i_elect = 2;
        i_part = 8;
        
        phase_tmp = x_phase{i_part,i_elect};
        figure; rose(squeeze(circ_mean(phase_tmp(:,:),[],1)))
        ylabel('Count');
        title(['X Subj ' num2str(exp.participants{i_part}) ': ' exp.singtrlelec_name{i_elect} '  ' num2str(freqs(freqlim(i_freq))) ' Hz'])
        clear phase_tmp
        phase_tmp = n_phase{i_part,i_elect};
        figure; rose(squeeze(circ_mean(phase_tmp(:,:),[],1)))
        ylabel('Count');
        title(['N Subj ' num2str(exp.participants{i_part}) ': ' exp.singtrlelec_name{i_elect} '  ' num2str(freqs(freqlim(i_freq))) ' Hz'])
        clear phase_tmp






for i_part = 1:length(exp.participants) % --------------
    for i_elect = 1:length(exp.singletrialselecs)
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),this_pwr{i_part}(i_elect,:));
        % correlation betwen power and errors and plot
        figure; polarscatter(circ_ang2rad(resp_errdeg{i_part}),this_pwr{i_part}(i_elect,:))
        title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
     end
end
% ---------------------


i_part = 1; %participant
i_elect = 1;
% all_ersp is {participant x electrode}(freq x time x trial)
figure; imagesc(times,freqs,angle(all_ersp{i_part,i_elect}(:,:,1))); set(gca,'ydir','normal'); colorbar;

figure; rose(squeeze(angle(all_ersp{i_part,i_elect}(12,56,:))))











