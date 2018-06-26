
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
%     ylim([2 30])
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
%%   Get BEH Data
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
%%   Compute alpha power on each trial for each participant by electrode
% /////////////////////////////////////////////////////////////////////////

%finds the frequencies you want from the freqs variable
% freqband = [7 14]; %alpha
freqband = [5 7]; %theta
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)));
%finds the times you want from the timess variable
% timewin = [50 250];
% timewin = [-600 -500];
timewin = [550 650];
timelim = find(times>=timewin(1) & times<=timewin(2));

current_pwr_lim = cell(length(exp.participants),length(exp.singletrialselecs)+1); %pre-allocate
tmp_pwr_lim = cell(1,length(exp.singletrialselecs));%pre-allocate
for i_part = 1:length(exp.participants)
    for i_elect = 1:length(exp.singletrialselecs)
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_pwr_lim{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial
        
        % plot trial power on a histogram
%         figure; hist(current_power{i_part,i_elect},30)
%         ylabel('Count');
%         title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}])
    end
    % Average power for O1 & O2
    current_pwr_lim{i_part,6}(:)=(current_pwr_lim{i_part,1}(:)+current_pwr_lim{i_part,4}(:))./2;
end
% clear i_elect
exp.singtrlelec_name{6} = 'Avg(O1,O2)';


% /////////////////////////////////////////////////////////////////////////
%%        Correlate alpha power with degrees error
% /////////////////////////////////////////////////////////////////////////

% Plot Correlations
for i_part = 1:length(exp.participants) % --------------
    for i_elect = 2    %length(exp.singletrialselecs)+1
        [rho,pval] = circ_corrcl(circ_ang2rad(resp_errdeg{i_part}),current_pwr_lim{i_part,i_elect});
        % correlation betwen power and errors and plot
        figure; polarscatter(circ_ang2rad(resp_errdeg{i_part}),current_pwr_lim{i_part,i_elect})
        title(['Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ' [' num2str(freqband(1)) '-'  num2str(freqband(2)) ' Hz]: rho=' num2str(round(rho,2))...
            ' pval=' num2str(round(pval,2))])
        clear x y rho pval
     end
end
% ---------------------


% /////////////////////////////////////////////////////////////////////////
% Separate by power by color trials
% for i_part = 1:4 % --------
for i_part = 5:length(exp.participants) % --------    
    % Split trials by target color
    error_deg_g{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 1);
    error_deg_b{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 2);
    error_deg_r{i_part} = resp_errdeg{i_part}(targ_colorgrp{i_part} == 3);

    for i_elect = 1:length(exp.singletrialselecs)
        
        pwr_g{i_part,i_elect} = current_pwr_lim{i_part,i_elect}(targ_colorgrp{i_part} == 1);
        pwr_b{i_part,i_elect} = current_pwr_lim{i_part,i_elect}(targ_colorgrp{i_part} == 2);
        pwr_r{i_part,i_elect} = current_pwr_lim{i_part,i_elect}(targ_colorgrp{i_part} == 3);
         
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

clear current_pwr_lim

% /////////////////////////////////////////////////////////////////////////
%% Separate trials by +/-1SD errors
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
freqlim = find(freqs>=(freqband(1)-0.5) & freqs<=(freqband(2)));
%finds the times you want from the timess variable
timewin = [-500 -250];
timelim = find(times>=timewin(1) & times<=timewin(2));

% only_ths_elect = [1,7,8,15]; %from exp.singletrialselecs

x_pwr_meanerr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
n_pwr_meanerr = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
% Calculate and scatter polar plot
for i_part = 1:3    %length(exp.participants) % ------------
    % Get error +/- SDs of the mean
    mean_errdeg = nanmean(resp_errdeg{i_part});
    sd_errdeg = nanstd(resp_errdeg{i_part});
    uplim = mean_errdeg + (sd_errdeg);
    lowlim = mean_errdeg - (sd_errdeg);

    x_errdeg = [resp_errdeg{i_part}(resp_errdeg{i_part}>=uplim) resp_errdeg{i_part}(resp_errdeg{i_part}<=lowlim)];
    n_errdeg = [resp_errdeg{i_part}(resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim)];

    % Calculate power
    for i_elect = 1:length(exp.singletrialselecs)-1%not FCz
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_power{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial

        % get power trials
        x_pwr_meanerr{i_part,i_elect} = [current_power{i_part,i_elect}(resp_errdeg{i_part}>=uplim) current_power{i_part,i_elect}(resp_errdeg{i_part}<=lowlim)];
        n_pwr_meanerr{i_part,i_elect} = current_power{i_part,i_elect}(resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim);

        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(x_errdeg),x_pwr_meanerr{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(x_errdeg),x_pwr_meanerr{i_part,i_elect})
        title(['X Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(n_errdeg),n_pwr_meanerr{i_part,i_elect});
        figure; polarscatter(circ_ang2rad(n_errdeg),n_pwr_meanerr{i_part,i_elect})
        title(['N Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ': rho=' num2str(round(rho,2)) ' pval=' num2str(round(pval,2))])
        clear x y rho pval
    end
    clear x_errdeg n_errdeg mean_errdeg sd_errdeg uplim lowlim
end
clear current_power


% ====== Get means for statistics ============
for i_part = 1:length(exp.participants)
    for i_elect = 1:length(exp.singletrialselecs)
        mean_n_pwr_byerr(i_part,i_elect) = mean([n_pwr_meanerr{i_part,i_elect}]);
        mean_x_pwr_byerr(i_part,i_elect) = mean([x_pwr_meanerr{i_part,i_elect}]);
    end
end
%t-test
[h p] = ttest(mean_n_pwr_byerr(:,5),mean_x_pwr_byerr(:,5))


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


% /////////////////////////////////////////////////////////////////////////

%finds the frequencies you want from the freqs variable
% freqband = [7 14]; %alpha
freqband = [5 7]; %theta
freqlim = find(freqs>=(freqband(1)-0.25) & freqs<=(freqband(2)+0.25));
%finds the times you want from the timess variable
% timewin = [-400 -150];
timewin = [100 650];
timelim = find(times>=timewin(1) & times<=timewin(2));

current_splitpower = cell(length(exp.participants),length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % -------------
    % Calculate power
    for i_elect = 1:length(exp.singletrialselecs)
        % all_ersp is (participant x electrode).trials(freq x time x trial)
        part_ersp = all_ersp{i_part,i_elect}; %get single subject's ersp
        for i_trial = 1:(size(part_ersp,3))
            current_splitpower{i_part,i_elect}(i_trial) = squeeze(mean(mean(abs(part_ersp(freqlim,timelim,i_trial)),1),2));
        end
        clear part_ersp i_trial
    end
    % Average power for O1 & O2
    current_splitpower{i_part,6}(:)=(current_splitpower{i_part,1}(:)+current_splitpower{i_part,4}(:))./2;
    exp.singtrlelec_name{6} = 'Avg(O1,O2)';
end


% correlation betwen power and errors and plot        
for i_part = 1:length(exp.participants) % -------------
    for i_elect = 2    %length(exp.singletrialselecs)  
        % Get trials by median power split
        mdn_pwr = nanmedian(current_splitpower{i_part,i_elect});
        h_errdeg{i_part,i_elect} = resp_errdeg{i_part}(current_splitpower{i_part,i_elect}>=mdn_pwr);
        l_errdeg{i_part,i_elect} = resp_errdeg{i_part}(current_splitpower{i_part,i_elect}<mdn_pwr);
        % get power trials
        h_pwr{i_elect} = current_splitpower{i_part,i_elect}(current_splitpower{i_part,i_elect}>=mdn_pwr);
        l_pwr{i_elect} = current_splitpower{i_part,i_elect}(current_splitpower{i_part,i_elect}<mdn_pwr);
        
        [rho,pval] = circ_corrcl(circ_ang2rad(h_errdeg{i_part,i_elect}),h_pwr{i_elect});
        figure; polarscatter(circ_ang2rad(h_errdeg{i_part,i_elect}),h_pwr{i_elect})
        title(['High Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ' [' num2str(freqband(1)) '-'  num2str(freqband(2)) ' Hz]: rho=' num2str(round(rho,2))...
            ' pval=' num2str(round(pval,2))])
        clear x y rho pval
        % correlation betwen power and errors and plot
        [rho,pval] = circ_corrcl(circ_ang2rad(l_errdeg{i_part,i_elect}),l_pwr{i_elect});
        figure; polarscatter(circ_ang2rad(l_errdeg{i_part,i_elect}),l_pwr{i_elect})
        title(['Low Subj ' num2str(exp.participants{i_part}) '-' exp.singtrlelec_name{i_elect}...
            ' [' num2str(freqband(1)) '-'  num2str(freqband(2)) ' Hz]: rho=' num2str(round(rho,2))...
            ' pval=' num2str(round(pval,2))])
        clear x y rho pval
    end
    clear h_pwr l_pwr mean_errdeg sd_errdeg uplim lowlim
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

% When fit fails
model_out_h{i_part,i_elect} = MemFit(h_errdeg{i_part,i_elect},model); %just plots
g_out_h(i_part,i_elect) = NaN;
sd_out_h(i_part,i_elect) = NaN;
% When fit fails
model_out_l{i_part,i_elect} = MemFit(l_errdeg{i_part,i_elect},model); %just plots
g_out_l(i_part,i_elect) = NaN;
sd_out_l(i_part,i_elect) = NaN;



 
nanmean(g_out_h(:,i_elect))
nanmean(g_out_l(:,i_elect))

[h p] = ttest(g_out_h(:,i_elect),g_out_l(:,i_elect))

nanmean(sd_out_h(:,i_elect))
nanmean(sd_out_l(:,i_elect))

[h p] = ttest(sd_out_h(:,i_elect),sd_out_l(:,i_elect))


figure; boxplot([g_out_h(:,i_elect),g_out_l(:,i_elect)])
figure; boxplot([sd_out_h(:,i_elect),sd_out_l(:,i_elect)])



figHand = PlotModelParametersAndData(StandardMixtureModel,...
    model_out_l{i_part,i_elect}.posteriorSamples,l_errdeg{i_part,i_elect}) 



clear h_errdeg l_errdeg sd_out_l g_out_l sd_out_h g_out_h
clear sd_out_l g_out_l sd_out_h g_out_h




% /////////////////////////////////////////////////////////////////////////
%% TF Spectogram: alpha power by errors - mean
% /////////////////////////////////////////////////////////////////////////

x_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
n_errdeg_m = cell(1,length(exp.participants)); %pre-allocate
x_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
n_pwr = cell(1,length(exp.singletrialselecs)); %pre-allocate
for i_part = 1:length(exp.participants) % ----------------------
    
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

% Make average O1 and O2
x_pwr{1,6} = (x_pwr{1,1}+x_pwr{1,4})./2;
n_pwr{1,6} = (n_pwr{1,1}+n_pwr{1,4})./2;
exp.singtrlelec_name{6} = 'Avg(O1,O2)';



% Plot Spectogram
for i_elect = 1:(length(exp.singletrialselecs)+1) %add one for Avg O1 & O2
    plot_ers_x = squeeze(mean(x_pwr{1,i_elect}(:,:,:),1)); %mean across subjs
    plot_ers_n = squeeze(mean(n_pwr{1,i_elect}(:,:,:),1)); %mean across subjs
    
    CLim = [-0.7 0.7];
    xmin = -600; xmax = 1200;
    % Plot Large-Small Errors
    figure; colormap('jet')
    imagesc(times,freqs,plot_ers_x-plot_ers_n,CLim);
    title(['ERS Difference: Large-Small Errors: ' exp.singtrlelec_name{i_elect}]); set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineWidth',1.5) %vertical line
    line([50 50],[min(freqs) max(freqs)],'color','b','LineStyle',':','LineWidth',1.5) %vertical line for mask onset
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5) %vertical line for color wheel
    ylim([1 30]); xlim([xmin xmax]); xticks(xmin:100:xmax)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    colorbar
    
    clear plot_ers_x plot_ers_n
end
clear ii i_elect

% (EEG.urevent(39).latency-EEG.urevent(37).latency)*(1/EEG.srate)

% select electrodes (bases on exp.singletrialselecs list)
% i_elect = 2; %Pz
% i_elect = 3; %Cz
% i_elect = 5; %FCz
% i_elect = 6; %O1 & O2
% i_elect = 1; %O1



%finds the frequencies you want from the freqs variable
% freqband = [7 14]; %alpha
freqband = [5 7]; %theta
% freqband = [1 4]; %delta
freqlim = find(freqs>=(freqband(1)-0.25) & freqs<=(freqband(2)+0.25));

%finds the times you want from the timess variable
% timewin = [50 250]; i_elect = 2; %Pz %alpha (p=0.0385)
% timewin = [0 250]; i_elect = 2; %Pz %alpha (p=0.0412)
% timewin = [550 650]; i_elect = 2; %Pz %theta (p=0.0370)
% timewin = [100 250]; i_elect = 2; %Pz %theta (p=0.0062)

% timewin = [-350 -100]; i_elect = 1; %O1 %alpha (p=0.0153)
% timewin = [50 250]; i_elect = 1; %O1 %alpha (p=0.0483)
% timewin = [50 200]; i_elect = 1; %O1 %alpha (p=0.0466)
% timewin = [50 200]; i_elect = 1; %O1 %theta (p=0.0073)
% timewin = [550 650]; i_elect = 1; %O1 %delta (p=0.0493)

% timewin = [-350 150]; i_elect = 4; %O2 %alpha (p=0.0140)
% timewin = [-200 150]; i_elect = 4; %O2 %alpha (p=0.0275)
% timewin = [950 1050]; i_elect = 4; %O2 %theta (p=0.0056)
% timewin = [900 1050]; i_elect = 4; %O2 %theta (p=0.0078)

% timewin = [-400 200]; i_elect = 6; %O1 & O2 %alpha (p=0.0183)
% timewin = [-400 -50]; i_elect = 6; %O1 & O2 %alpha (p=0.0173)
% timewin = [950 1050]; i_elect = 6; %O1 & O2 %theta (p=0.0092)
% timewin = [-400 -100]; i_elect = 6; %O1 & O2 %alpha (p=0.0163)
timelim = find(times>=timewin(1) & times<=timewin(2));

plot_ers_x = squeeze(mean(mean(x_pwr{1,i_elect}(:,freqlim,timelim),2),3)); 
plot_ers_n = squeeze(mean(mean(n_pwr{1,i_elect}(:,freqlim,timelim),2),3)); 

[h p] = ttest(plot_ers_x,plot_ers_n)


















% #########################################################################
%%     Correlate alpha power with degrees error BY COLOR
% #########################################################################

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

% #########################################################################
% #########################################################################






