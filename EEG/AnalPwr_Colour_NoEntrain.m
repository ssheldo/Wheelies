
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%                             INFORMATION
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

% SPECTOGRAM
% A spectogram is a 3d figure that plots time on the x-axis, frequency on the 
% y-axis, and shows you the power or phase-locking value for each point. 
% We compute spectograms if we have power and phase information, averaged 
% across trials, for at least one electrode. 
% This can help us understand the changes of power and phase throughout the 
% trial.

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
%%     ERSP plots averaged over subjects and channels
% /////////////////////////////////////////////////////////////////////////


% (participants x conditions x events x electrodes x frequencies x timepoints)
out_grand_ersp = squeeze(mean(mean(ersp(:,1,1,:,:,:),4),1));


figure; 
CLim = [0 20];
colormap('jet')

imagesc(times,freqs,(out_grand_ersp));
title('ERS: Grand Mean'); 
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
ylim([2 30])
ylabel('Freq (Hz)');
xlabel('Time (ms)');
colorbar



% /////////////////////////////////////////////////////////////////////////
%%         ERSP plots averaged over subjects for each channel
% /////////////////////////////////////////////////////////////////////////

CLim = [0 20];

for ii = 1:length(exp.electrode)

    i_elect = exp.electrode(ii); %get electrode number
    
    % (participants x conditions x events x electrodes x frequencies x timepoints)
    out_ers = squeeze(mean(ersp(:,1,1,i_elect,:,:),1));

    figure; 
    colormap('jet')

    imagesc(times,freqs,(out_ers));
    title(['ERS: ' exp.elec_names{ii}]); 
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
    ylim([2 30])
    ylabel('Freq (Hz)');
    xlabel('Time (ms)');
    colorbar

    clear out_ers
end
clear ii i_elect







% /////////////////////////////////////////////////////////////////////////
%%         ITC plots averaged over subjects for each channel
% /////////////////////////////////////////////////////////////////////////

for ii = 1:length(exp.electrode)
    i_elect = exp.electrode(ii); %get electrode number
    % (participants x conditions x events x electrodes x frequencies x timepoints)
    out_ers = squeeze(mean(abs(itc(:,1,1,i_elect,:,:)),1));

    figure; 
    colormap('jet')

    imagesc(times,freqs,(out_ers));
    title(['ITC: ' exp.elec_names{ii}]); 
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
    ylim([2 30])
    ylabel('Freq (Hz)');
    xlabel('Time (ms)');
    colorbar

    clear out_ers
end

clear ii i_elect










% /////////////////////////////////////////////////////////////////////////
%%   Compute alpha power on each trial for each participant by electrode
% /////////////////////////////////////////////////////////////////////////

%finds the frequencies you want from the freqs variable
freqband = [7 14]; %alpha
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-500 -250];
timelim = find(times>=timewin(1) & times<=timewin(2));

current_power = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants)
    for i_elect = 1:length(exp.singletrialselecs)
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial

        % plot trial power on a histogram
        figure; hist(current_power{i_part,i_elect},30)
        ylabel('Count');
        title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}])
    end
end
clear i_elect


% /////////////////////////////////////////////////////////////////////////
%%        Correlate alpha power with degrees error
% /////////////////////////////////////////////////////////////////////////

% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
targ_colorgrp = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % ------------------------------------
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

%convert the -deg errors into 360 angle coordinates
% err_angle = cell(length(exp.participants),1); %pre-allocate
% for i_part = 1:length(exp.participants) % ------------------------------------
%     tmp_err = resp_errdeg{i_part};
%     for err = 1:length(tmp_err)
%         if tmp_err(err) < 0
%             err_angle{i_part}(err) = 360 + tmp_err(err);
%         else
%             err_angle{i_part}(err) = tmp_err(err);
%         end
%     end
%     clear err tmp_err
% end
% % -------------------------------------------------------------------------


for i_part = 1:length(exp.participants) % --------------
    for i_elect = 1:length(exp.singletrialselecs)
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect});
        % correlation betwen power and errors and plot
        figure; polarscatter(circ_ang2rad(resp_errdeg{i_part}),current_power{i_part,i_elect})
        title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
     end
end
% ---------------------


% Separate by color
% for i_part = 1:4 % --------
for i_part = 5:length(exp.participants) % --------    
    % Split trials by target color
    error_deg_g{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 1);
    error_deg_b{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 2);
    error_deg_r{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 3);

    for i_elect = 1:length(exp.singletrialselecs)
        
        pwr_g{i_part,i_elect} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 1);
        pwr_b{i_part,i_elect} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 2);
        pwr_r{i_part,i_elect} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 3);
         
        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(error_deg_g{i_part}),pwr_g{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(error_deg_g{i_part}),pwr_g{i_part,i_elect})
        title(['Green Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear rho pval
        
        [rho,pval] = circ_corrcl(circ_ang2rad(error_deg_b{i_part}),pwr_b{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(error_deg_b{i_part}),pwr_b{i_part,i_elect})
        title(['Blue Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear rho pval
        
        [rho,pval] = circ_corrcl(circ_ang2rad(error_deg_r{i_part}),pwr_r{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(error_deg_r{i_part}),pwr_r{i_part,i_elect})
        title(['Red Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear rho pval
        
     end
end

% ---------------------------


% /////////////////////////////////////////////////////////////////////////
%% Separate trials by +/-2SD errors
% /////////////////////////////////////////////////////////////////////////

% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
targ_colorgrp = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % ------------------------------------
    [n,m] = size(ALLEEG(i_part).rejtrial);
    % Get list of rejected trials
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

%finds the frequencies you want from the freqs variable
freqband = [7 14]; %alpha
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-500 -250];
timelim = find(times>=timewin(1) & times<=timewin(2));
for i_part = 5:length(exp.participants) % ------------------------------------
% for i_part = 1:4 % ------------------------------------
    % Get error +/-2 SDs of the mean
    mean_errdeg = nanmedian(resp_errdeg{i_part});
    sd_errdeg = iqr(resp_errdeg{i_part})/2;
    uplim = mean_errdeg + (sd_errdeg);
    lowlim = mean_errdeg - (sd_errdeg);

    x_errdeg = [resp_errdeg{i_part}(resp_errdeg{i_part}>=uplim) resp_errdeg{i_part}(resp_errdeg{i_part}<=lowlim)];
    n_errdeg = [resp_errdeg{i_part}(resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim)];

    % Calculate power
    for i_elect = 1:length(exp.singletrialselecs)
%     for eek = 1:length(exp.singletrialselecs)    
%         i_elect = exp.singletrialselecs(eek);
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial

        % get power trials
        x_pwr{i_part,i_elect} = [current_power{i_part,i_elect}(resp_errdeg{i_part}>=uplim) current_power{i_part,i_elect}(resp_errdeg{i_part}<=lowlim)];
        n_pwr{i_part,i_elect} = current_power{i_part,i_elect}(resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim);

        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(x_errdeg),x_pwr{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(x_errdeg),x_pwr{i_part,i_elect})
        title(['X Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(n_errdeg),n_pwr{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(n_errdeg),n_pwr{i_part,i_elect})
        title(['N Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
    end
    clear x_errdeg n_errdeg mean_errdeg sd_errdeg uplim lowlim
end

% Get means for statistics
for i_part = 1:length(exp.participants)
    for i_elect = 1:length(exp.singletrialselecs)
        mean_n_pwr(i_part,i_elect) = median([n_pwr{i_part,i_elect}]);
        mean_x_pwr(i_part,i_elect) = median([x_pwr{i_part,i_elect}]);
    end
end
%t-test
[h p] = ttest(mean_n_pwr(:,5),mean_x_pwr(:,5))


clear n_pwr x_pwr mean_n_pwr mean_x_pwr



% /////////////////////////////////////////////////////////////////////////
%% Separate trials by high/low alpha power
% /////////////////////////////////////////////////////////////////////////

% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
targ_colorgrp = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % ------------------------------------
    [n,m] = size(ALLEEG(i_part).rejtrial);
    % Get list of rejected trials
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

%finds the frequencies you want from the freqs variable
% freqband = [4 7]; %theta
freqband = [10 14]; %alpha
% freqband = [19 24]; %beta1
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-500 -250];
timelim = find(times>=timewin(1) & times<=timewin(2));

% for i_part = 5:length(exp.participants) % ------------------------------------
for i_part = 1:4 % ------------------------------------
    % Calculate power
    for i_elect = 1:length(exp.singletrialselecs)
%     for eek = 1:length(exp.singletrialselecs)    
%         i_elect = exp.singletrialselecs(eek);
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial
        
        % Get power +/-2 IQRs of the median
        mdn_pwr = nanmedian(current_power{i_part,i_elect});

        h_errdeg{i_part,i_elect} = resp_errdeg{i_part}(current_power{i_part,i_elect}>=mdn_pwr);
        l_errdeg{i_part,i_elect} = resp_errdeg{i_part}(current_power{i_part,i_elect}<mdn_pwr);

        % get power trials
        h_pwr{i_elect} = current_power{i_part,i_elect}(current_power{i_part,i_elect}>=mdn_pwr);
        l_pwr{i_elect} = current_power{i_part,i_elect}(current_power{i_part,i_elect}<mdn_pwr);

        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(h_errdeg{i_part,i_elect}),h_pwr{i_elect});
        figure; polarscatter(circ_ang2rad(h_errdeg{i_part,i_elect}),h_pwr{i_elect})
        title(['Theta High Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(l_errdeg{i_part,i_elect}),l_pwr{i_elect});
        figure; polarscatter(circ_ang2rad(l_errdeg{i_part,i_elect}),l_pwr{i_elect})
        title(['Theta Low Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
    end
    clear h_pwr l_pwr mean_errdeg sd_errdeg uplim lowlim
end


% /////////////////////////////////////////////////////////////////////////
%% Model fit error by high/low alpha power
% /////////////////////////////////////////////////////////////////////////

model = AllGuessingModel; %for when the fitting doesn't work
for i_part = 1:length(exp.participants) % ------------------------------------
  
    for i_elect = 2
        model_out_h{i_part,i_elect} = MemFit(h_errdeg{i_part,i_elect});
%         model_out_h{i_part,i_elect} = MemFit(h_errdeg{i_part,i_elect},model);
        model_out_l{i_part,i_elect} = MemFit(l_errdeg{i_part,i_elect});
%         model_out_l{i_part,i_elect} = MemFit(l_errdeg{i_part,i_elect},model);
        

        % Save output from model fit
        g_out_h(i_part,i_elect) = model_out_h{i_part,i_elect}.maxPosterior(1);
        sd_out_h(i_part,i_elect) = model_out_h{i_part,i_elect}.maxPosterior(2);
        g_out_l(i_part,i_elect) = model_out_l{i_part,i_elect}.maxPosterior(1);
        sd_out_l(i_part,i_elect) = model_out_l{i_part,i_elect}.maxPosterior(2);    
    end
end


nanmean(g_out_h(:,i_elect))
nanmean(g_out_l(:,i_elect))

[h p] = ttest(g_out_h(:,i_elect),g_out_l(:,i_elect))

nanmean(sd_out_h(:,i_elect))
nanmean(sd_out_l(:,i_elect))

[h p] = ttest(sd_out_h(:,i_elect),sd_out_l(:,i_elect))


figure; boxplot([g_out_h(:,i_elect),g_out_l(:,i_elect)])


clear h_errdeg l_errdeg sd_out_l g_out_l sd_out_h g_out_h
clear sd_out_l g_out_l sd_out_h g_out_h


% #########################################################################
% #########################################################################
%%                               BY COLOR
% #########################################################################
% #########################################################################


% /////////////////////////////////////////////////////////////////////////
%%        Correlate alpha power with degrees error
% /////////////////////////////////////////////////////////////////////////

% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
targ_colorgrp = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % ------------------------------------
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

% Separate by color
% for i_part = 1:4 % --------
for i_part = 5:length(exp.participants) % --------    
    % Split trials by target color
    error_deg_g{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 1);
    error_deg_b{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 2);
    error_deg_r{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 3);

    for i_elect = 1:length(exp.singletrialselecs)
        
        pwr_g{i_part,i_elect} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 1);
        pwr_b{i_part,i_elect} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 2);
        pwr_r{i_part,i_elect} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 3);
         
        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(error_deg_g{i_part}),pwr_g{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(error_deg_g{i_part}),pwr_g{i_part,i_elect})
        title(['Green Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear rho pval
        
        [rho,pval] = circ_corrcl(circ_ang2rad(error_deg_b{i_part}),pwr_b{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(error_deg_b{i_part}),pwr_b{i_part,i_elect})
        title(['Blue Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear rho pval
        
        [rho,pval] = circ_corrcl(circ_ang2rad(error_deg_r{i_part}),pwr_r{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(error_deg_r{i_part}),pwr_r{i_part,i_elect})
        title(['Red Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear rho pval
        
     end
end

% ---------------------------


% /////////////////////////////////////////////////////////////////////////
%% Separate trials by +/-2SD errors
% /////////////////////////////////////////////////////////////////////////

% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
targ_colorgrp = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % ------------------------------------
    [n,m] = size(ALLEEG(i_part).rejtrial);
    % Get list of rejected trials
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

%finds the frequencies you want from the freqs variable
freqband = [7 14]; %alpha
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-500 -250];
timelim = find(times>=timewin(1) & times<=timewin(2));

% Separate errors and power trials by color   
for i_part = 1:length(exp.participants) % ------------------------------------    
    
    % Split trials by target color
    errdeg_c{i_part,1} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 1);
    errdeg_c{i_part,2} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 2);
    errdeg_c{i_part,3} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 3);
    
    for i_elect = 1:length(exp.singletrialselecs) % Calculate power
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial

        % Split trials by color
        pwr_c{i_part,i_elect,1} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 1);
        pwr_c{i_part,i_elect,2} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 2);
        pwr_c{i_part,i_elect,3} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 3);
    end
end    

cname = {'Green';'Blue';'Red'};
for i_color = 1:3 %green, blue, red
%     for i_part = 1:length(exp.participants) % ------------------------------------
    for i_part = 8 % ------------------------------------

        % Get error +/-2 SDs of the mean
        mean_errdeg = nanmedian(errdeg_c{i_part,i_color});
        sd_errdeg = iqr(errdeg_c{i_part,i_color})/2;
        uplim = mean_errdeg + sd_errdeg;
        lowlim = mean_errdeg - sd_errdeg;

        x_errdeg = [errdeg_c{i_part,i_color}(errdeg_c{i_part,i_color}>=uplim) errdeg_c{i_part,i_color}(errdeg_c{i_part,i_color}<=lowlim)];
        n_errdeg = [errdeg_c{i_part,i_color}(errdeg_c{i_part,i_color}<uplim & errdeg_c{i_part,i_color}>lowlim)];

        % Calculate power
        for i_elect = 1:length(exp.singletrialselecs)

            % get power trials
            x_pwr{i_part,i_elect,i_color} = [pwr_c{i_part,i_elect,i_color}(errdeg_c{i_part,i_color}>=uplim) pwr_c{i_part,i_elect,i_color}(errdeg_c{i_part,i_color}<=lowlim)];
            n_pwr{i_part,i_elect,i_color} = pwr_c{i_part,i_elect,i_color}(errdeg_c{i_part,i_color}<uplim & errdeg_c{i_part,i_color}>lowlim);

            % correlation betwen power and errors and plot
            [rho,pval] = circ_corrcl(circ_ang2rad(x_errdeg),x_pwr{i_part,i_elect,i_color});
            figure; polarscatter(circ_ang2rad(x_errdeg),x_pwr{i_part,i_elect,i_color})
            title([cname{i_color} ' X Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
                ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
            clear x y rho pval
            % correlation betwen power and errors and plot
            [rho,pval] = circ_corrcl(circ_ang2rad(n_errdeg),n_pwr{i_part,i_elect,i_color});
            figure; polarscatter(circ_ang2rad(n_errdeg),n_pwr{i_part,i_elect,i_color})
            title([cname{i_color} ' N Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
                ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
            clear x y rho pval
        end
        clear x_errdeg n_errdeg mean_errdeg sd_errdeg uplim lowlim
    end
end



% Get means for statistics
for i_part = 1:length(exp.participants)
    for i_elect = 1:length(exp.singletrialselecs)
        mean_n_pwr(i_part,i_elect,i_color) = median([n_pwr{i_part,i_elect}]);
        mean_x_pwr(i_part,i_elect,i_color) = median([x_pwr{i_part,i_elect}]);
    end
end
%t-test
[h p] = ttest(mean_n_pwr(:,5),mean_x_pwr(:,5))


clear n_pwr x_pwr mean_n_pwr mean_x_pwr



% /////////////////////////////////////////////////////////////////////////
%% Separate trials by high/low alpha power
% /////////////////////////////////////////////////////////////////////////

% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
targ_colorgrp = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants) % ------------------------------------
    [n,m] = size(ALLEEG(i_part).rejtrial);
    % Get list of rejected trials
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

%finds the frequencies you want from the freqs variable
freqband = [7 14]; %alpha
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)+0.5));
%finds the times you want from the timess variable
timewin = [-500 -250];
timelim = find(times>=timewin(1) & times<=timewin(2));
% Separate errors and power trials by color % ------------------------------------   
for i_part = 1:length(exp.participants)    
    % Split trials by target color
    errdeg_c{i_part,1} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 1);
    errdeg_c{i_part,2} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 2);
    errdeg_c{i_part,3} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 3);
    
    for i_elect = 1:length(exp.singletrialselecs) % Calculate power
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial

        % Split trials by color
        pwr_c{i_part,i_elect,1} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 1);
        pwr_c{i_part,i_elect,2} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 2);
        pwr_c{i_part,i_elect,3} = current_power{i_part,i_elect}(targ_colorgrp{i_part} == 3);
    end
end  
% -------------------------------------------------------------------------

cname = {'Green';'Blue';'Red'};
% for i_part = 1:length(exp.participants) % ------------------------------------
for i_part = 5:8 % ------------------------------------
    for i_color = 1:3 %green, blue, red
        % Calculate power
        for i_elect = 3

            % Get power +/-2 IQRs of the median
            mdn_pwr = nanmedian(pwr_c{i_part,i_elect,i_color});

            h_errdeg{i_part,i_elect,i_color} = errdeg_c{i_part,i_color}(pwr_c{i_part,i_elect,i_color}>=mdn_pwr);
            l_errdeg{i_part,i_elect,i_color} = errdeg_c{i_part,i_color}(pwr_c{i_part,i_elect,i_color}<mdn_pwr);

            % get power trials
            h_pwr{i_part,i_elect,i_color} = pwr_c{i_part,i_elect,i_color}(pwr_c{i_part,i_elect,i_color}>=mdn_pwr);
            l_pwr{i_part,i_elect,i_color} = pwr_c{i_part,i_elect,i_color}(pwr_c{i_part,i_elect,i_color}<mdn_pwr);

            
            % correlation betwen power and errors and plot
            [rho,pval] = circ_corrcl(circ_ang2rad(h_errdeg{i_part,i_elect,i_color}),h_pwr{i_part,i_elect,i_color});
            figure; polarscatter(circ_ang2rad(h_errdeg{i_part,i_elect,i_color}),h_pwr{i_part,i_elect,i_color})
            title(['H ' cname{i_color} ' Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
                ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
            clear x y rho pval
            % correlation betwen power and errors and plot
            [rho,pval] = circ_corrcl(circ_ang2rad(l_errdeg{i_part,i_elect,i_color}),l_pwr{i_part,i_elect,i_color});
            figure; polarscatter(circ_ang2rad(l_errdeg{i_part,i_elect,i_color}),l_pwr{i_part,i_elect,i_color})
            title(['L ' cname{i_color} ' Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
                ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
            clear x y rho pval
            
        end
        clear mdn_pwr mean_errdeg sd_errdeg uplim lowlim
    end
end

% /////////////////////////////////////////////////////////////////////////
%% Model fit error by high/low alpha power
% /////////////////////////////////////////////////////////////////////////

model = AllGuessingModel; %for when the fitting doesn't work
for i_part = 4:length(exp.participants) % ------------------------------------
  
    for i_elect = 2
        model_out_h{i_part,i_elect} = MemFit(h_errdeg{i_part,i_elect});
%         model_out_h{i_part,i_elect} = MemFit(h_errdeg{i_part,i_elect},model);
        model_out_l{i_part,i_elect} = MemFit(l_errdeg{i_part,i_elect});
%         model_out_l{i_part,i_elect} = MemFit(l_errdeg{i_part,i_elect},model);
        

        % Save output from model fit
        g_out_h(i_part,i_elect) = model_out_h{i_part,i_elect}.maxPosterior(1);
        sd_out_h(i_part,i_elect) = model_out_h{i_part,i_elect}.maxPosterior(2);
        g_out_l(i_part,i_elect) = model_out_l{i_part,i_elect}.maxPosterior(1);
        sd_out_l(i_part,i_elect) = model_out_l{i_part,i_elect}.maxPosterior(2);    
    end
end


nanmean(g_out_h(:,i_elect))
nanmean(g_out_l(:,i_elect))

[h p] = ttest(g_out_h(:,i_elect),g_out_l(:,i_elect))

nanmean(sd_out_h(:,i_elect))
nanmean(sd_out_l(:,i_elect))

[h p] = ttest(sd_out_h(:,i_elect),sd_out_l(:,i_elect))





clear h_errdeg l_errdeg





% #########################################################################
% #########################################################################


% /////////////////////////////////////////////////////////////////////////
%% TF Spectogram: alpha power by errors
% /////////////////////////////////////////////////////////////////////////

x_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
n_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
x_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
n_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % ----------------------
% for i_part = 1:4 % ------------------------------------
    % Get error +/-1 SDs of the mean
%     mean_errdeg = nanmedian(resp_errdeg{i_part});
%     sd_errdeg = iqr(resp_errdeg{i_part})/2;
%     uplim = mean_errdeg + (sd_errdeg);
%     lowlim = mean_errdeg - (sd_errdeg);
    mean_errdeg = nanmean(resp_errdeg{i_part});
    sd_errdeg = nanstd(resp_errdeg{i_part});
    uplim = mean_errdeg + (sd_errdeg);
    lowlim = mean_errdeg - (sd_errdeg);

    x_errdeg_m{i_part} = [resp_errdeg{i_part}(resp_errdeg{i_part}>=uplim) resp_errdeg{i_part}(resp_errdeg{i_part}<=lowlim)];
    n_errdeg_m{i_part} = [resp_errdeg{i_part}(resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim)];

    % Calculate power
    for i_elect = 1:length(exp.singletrialselecs)
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = abs(all_ersp{i_part,i_elect}); %get single subject's power
        % get power on small/large error trials
        x_pwr{1,i_elect}(i_part,:,:) = squeeze(mean(part_ersp(:,:,[find(resp_errdeg{i_part}>=uplim) find(resp_errdeg{i_part}<=lowlim)]),3));
        n_pwr{1,i_elect}(i_part,:,:) = squeeze(mean(part_ersp(:,:,[find((resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim))]),3));
    end
    clear mean_errdeg sd_errdeg uplim lowlim part_ersp
end
clear ii i_elect


% Plot Spectogram
for i_elect = 1:length(exp.singletrialselecs)
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %mean across subjs
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %mean across subjs
    
    CLim = [-0.5 0.5];
    % Plot Large Errors
%     figure; colormap('jet')
%     imagesc(times,freqs,plot_ers_x,CLim);
%     title(['Large Errors: ' exp.singtrlelec_name{i_elect}]); set(gca,'Ydir','Normal')
%     line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) %vertical line
%     ylim([2 30])
%     ylabel('Freq (Hz)'); xlabel('Time (ms)');
%     colorbar
%     % Plot Small Errors
%     figure; colormap('jet')
%     imagesc(times,freqs,plot_ers_n,CLim);
%     title(['Small Errors: ' exp.singtrlelec_name{i_elect}]); set(gca,'Ydir','Normal')
%     line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) %vertical line
%     ylim([2 30])
%     ylabel('Freq (Hz)'); xlabel('Time (ms)');
%     colorbar
    % Plot Large-Small Errors
    figure; colormap('jet')
    imagesc(times,freqs,plot_ers_x-plot_ers_n,CLim);
    title(['Large-Small Errors: ' exp.singtrlelec_name{i_elect}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) %vertical line
    ylim([2 30])
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    
    clear plot_ers_x plot_ers_n
end
clear ii i_elect









