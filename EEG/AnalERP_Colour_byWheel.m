

%% BEH data
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
%%                  ERPs by Color - Mean
% /////////////////////////////////////////////////////////////////////////

for i_part = 1:length(exp.participants)   
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %because FCz is 17
        erp_out_c(i_part,ii,:,1) = mean(ALLEEG(i_part).data(i_elect,:,find(targ_colorgrp{i_part}==1)),3);
        erp_out_c(i_part,ii,:,2) = mean(ALLEEG(i_part).data(i_elect,:,find(targ_colorgrp{i_part}==2)),3);
        erp_out_c(i_part,ii,:,3) = mean(ALLEEG(i_part).data(i_elect,:,find(targ_colorgrp{i_part}==3)),3);
    end
end

erp_out_bycolor(:,:,1) = squeeze(mean(erp_out_c(:,:,:,1),1));
erp_out_bycolor(:,:,2) = squeeze(mean(erp_out_c(:,:,:,2),1));
erp_out_bycolor(:,:,3) = squeeze(mean(erp_out_c(:,:,:,3),1));
% cname = {'Green';'Blue';'Red'};

% Make average O1 and O2
erp_out_bycolor(17,:,1) = (erp_out_bycolor(1,:,1)+erp_out_bycolor(15,:,1))./2;
erp_out_bycolor(17,:,2) = (erp_out_bycolor(1,:,2)+erp_out_bycolor(15,:,2))./2;
erp_out_bycolor(17,:,3) = (erp_out_bycolor(1,:,3)+erp_out_bycolor(15,:,3))./2;
exp.elec_names{17} = 'Avg(O1,O2)';


% i_elect = 7; %Pz
% i_elect = 8; %Cz
% i_elect = 15; %O2
% i_elect = 1; %O1
i_elect = 16; %FCz
% i_elect = 17; %O1 & O2

figure
plot(EEG.times,erp_out_bycolor(i_elect,:,1),'g','LineWidth',1.5)
hold on
plot(EEG.times,erp_out_bycolor(i_elect,:,3),'-.r','LineWidth',1.5)
hold on
plot(EEG.times,erp_out_bycolor(i_elect,:,2),'b','LineWidth',1)
hold on
line([-500 1000],[0 0],'color','k') %horizontal line
line([0 0],[-4 8],'color','k') %vertical line
line([50 50],[-8 10],'color','b') %vertical line for mask onset
line([567 567],[-8 10],'color','b') %vertical line for color wheel onset
set(gca,'ydir','reverse'); xlim([-500 1000]); ylim([-4 8])
title([exp.elec_names{i_elect} ': ERPs by Color Trials']); 
xlabel('Time (ms)'); ylabel('Voltage (uV)'); 
legend('Green','Red','Blue'); xticks(-500:250:1000);

 
i_elect = 7; %Pz
% i_elect = 8; %Cz
% i_elect = 9; %Fz
% i_elect = 1; %O1
% i_elect = 15; %O2
% i_elect = 16; %FCz
ths_elect = [1;7;8;15;16];

for epp = 1:length(ths_elect)
    i_elect = ths_elect(epp);
    figure;
    boundedline(EEG.times,squeeze(mean(erp_out_c(:,i_elect,:,1),1)),...
        squeeze(std(erp_out_c(:,i_elect,:,1),[],1))./sqrt(nparts),'g',...
        EEG.times,squeeze(mean(erp_out_c(:,i_elect,:,2),1)),...
        squeeze(std(erp_out_c(:,i_elect,:,2),[],1))./sqrt(nparts),'r',...
        EEG.times,squeeze(mean(erp_out_c(:,i_elect,:,3),1)),...
        squeeze(std(erp_out_c(:,i_elect,:,3),[],1))./sqrt(nparts),'b');
    hold on
    set(gca,'ydir','reverse'); xlim([-500 1000]); ylim([-8 10])
    line([-500 1000],[0 0],'color','k') %horizontal line
    line([0 0],[-8 10],'color','k') %vertical line
    title([exp.elec_names{i_elect} ': ERPs by Response Error']); 
    xlabel('Time (ms)'); ylabel('Voltage (uV)'); xticks(-500:250:1000);
    legend('Green','Red','Blue');
    hold off
end
clear epp


% Select electrode
i_elect = 7; %Pz
% i_elect = 8; %Cz
% i_elect = 9; %Fz
% i_elect = 1; %O1
% i_elect = 15; %O2
% i_elect = 16; %FCz

% Set the range of time to consider
itlims = [-350 -200];
% itlims = [-350 -250];
% itlims = [0 200];
% itlims = [250 550];
% itlims = [500 600];
% itlims = [500 550];
% itlims = [800 999];

%this code finds the times you want from the timess variable
time_lims = find(EEG.times>= itlims(1),1):find(EEG.times>= itlims(2),1)-1;

erp_out_g = squeeze(mean(erp_out_c(:,:,time_lims,1),3));
erp_out_b = squeeze(mean(erp_out_c(:,:,time_lims,2),3));
erp_out_r = squeeze(mean(erp_out_c(:,:,time_lims,3),3));

[h p] = ttest(erp_out_g(:,i_elect), erp_out_r(:,i_elect))
[h p] = ttest(erp_out_b(:,i_elect), erp_out_r(:,i_elect))





% /////////////////////////////////////////////////////////////////////////
%%                  ERPs by Errors - Mean
% /////////////////////////////////////////////////////////////////////////

for i_part = 1:length(exp.participants) % ------------------------------------
    % Get error +/-1 SD of the mean
    mean_errdeg = mean(resp_errdeg{i_part});
    sd_errdeg = std(resp_errdeg{i_part});
    uplim = mean_errdeg + 1*(sd_errdeg);
    lowlim = mean_errdeg - 1*(sd_errdeg);
    % Calculate ERP
    for ii = 1:length(exp.electrode)
        i_elect = exp.electrode(ii); %because FCz is 17
        % Get trials with large errors
        erp_out_x(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,[find((resp_errdeg{i_part}<uplim & resp_errdeg{i_part}>lowlim))]),3));
        % Get trials with small errors
        erp_out_n(i_part,ii,:) = squeeze(mean(ALLEEG(i_part).data(i_elect,:,[find(resp_errdeg{i_part}>=uplim) find(resp_errdeg{i_part}<=lowlim)]),3));
    end
    clear mean_errdeg sd_errdeg uplim lowlim
end

% average across subjects
erp_out_byerr(:,:,1) = squeeze(mean(erp_out_x(:,:,:),1));
erp_out_byerr(:,:,2) = squeeze(mean(erp_out_n(:,:,:),1));

% Make average O1 and O2 across subjects
erp_out_byerr(17,:,1) = (erp_out_byerr(1,:,1)+erp_out_byerr(15,:,1))./2;
erp_out_byerr(17,:,2) = (erp_out_byerr(1,:,2)+erp_out_byerr(15,:,2))./2;
exp.elec_names{17} = 'Avg(O1,O2)';
% Make average O1 and O2 for each subjects
erp_out_x(:,17,:) = (erp_out_x(:,1,:)+erp_out_x(:,15,:))./2;
erp_out_n(:,17,:) = (erp_out_n(:,1,:)+erp_out_n(:,15,:))./2;


% (EEG.urevent(35).latency-EEG.urevent(32).latency)*(1/EEG.srate)

% Select electrode to plot
i_elect = 7; %Pz
% i_elect = 8; %Cz
% i_elect = 15; %O2
% i_elect = 1; %O1
% i_elect = 16; %FCz
% i_elect = 17; %O1 & O2

plot_ths = [7;8;16;17];

for ii = 1:length(plot_ths)
    i_elect = plot_ths(ii);
    % get axes limits
%     ymin = floor(min([erp_out_byerr(i_elect,:,1) erp_out_byerr(i_elect,:,2)]));
%     ymax = ceil(max([erp_out_byerr(i_elect,:,1) erp_out_byerr(i_elect,:,2)]));
    ymin = -5; ymax = 9;    
    xmin = -1400; xmax = 600;
    % Plotting ERPs
    figure
    plot(EEG.times,erp_out_byerr(i_elect,:,1),'c','LineWidth',1.5)
    hold on
    plot(EEG.times,erp_out_byerr(i_elect,:,2),'-.m','LineWidth',1.5)
    hold on
    line([xmin xmax],[0 0],'color','k') %horizontal line
    line([0 0],[ymin ymax],'color','k') %vertical line
%     line([-517 -517],[ymin ymax],'color','b') %vertical line for mask onset
    line([-567 -567],[ymin ymax],'color','b') %vertical line for target onset
    set(gca,'ydir','reverse'); xlim([xmin xmax]); ylim([ymin ymax])
    title([exp.elec_names{i_elect} ': ERPs by Response Error']); 
    xlabel('Time (ms)'); ylabel('Voltage (uV)'); xticks(xmin:200:xmax);
    legend('Large','Small');
    clear ymin ymax
end


% +++++++ Statistics ++++++++
% Select electrode
% i_elect = 7; %Pz
% i_elect = 8; %Cz
% i_elect = 9; %Fz
% i_elect = 1; %O1
% i_elect = 15; %O2
% i_elect = 16; %FCz
% i_elect = 17; %O1 & O2

% Set the range of time to consider
% itlims = [-750 -450]; i_elect = 7; %Pz (p=0.0214)
% itlims = [-800 -500]; i_elect = 7; %Pz (p=0.0116)
% itlims = [-800 -400]; i_elect = 8; %Cz (p=0.0294)
% itlims = [-300 -200]; i_elect = 16; %FCz (p=0.0424)
% itlims = [-300 0]; i_elect = 16; %FCz (p=0.0437)
% itlims = [-800 -500]; i_elect = 17; %O1 & O2 (p=0.0444)

% itlims = [500 600];
% itlims = [450 650];
% itlims = [-750 -550];

%this code finds the times you want from the timess variable
time_lims = find(EEG.times>= itlims(1),1):find(EEG.times>= itlims(2),1)-1;

erp_err_x = squeeze(mean(erp_out_x(:,:,time_lims),3));
erp_err_n = squeeze(mean(erp_out_n(:,:,time_lims),3));


[h p] = ttest(erp_err_x(:,i_elect), erp_err_n(:,i_elect))



% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% '''''''''''''''''''''''''    Topographys     '''''''''''''''''''''''''''
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


% Set the range of time to consider
tWin{1} = [100 200];
tWin{2} = [200 300];
tWin{3} = [300 400];
tWin{4} = [400 500];
tWin{5} = [500 600];
% tWin{6} = [450 600];

% tWin{1} = [175 275];

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% erp_out(time_window, standard/target, electrode, in/out, subjects)
erp_diff_out = squeeze(erp_out(:,2,:,:,:)-erp_out(:,1,:,:,:)); 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 % Topographys (Targets - Standards)
CLims = [-7 7];
 
for tw_i = 1:length(tWin) 
 
    itWin = tWin{tw_i}; %select each time range if looping
    %this code finds the times you want from the timess variable
    time_window = find(EEG.times>= itWin(1),1):find(EEG.times>= itWin(2),1)-1;
    
    figure('Color',[1 1 1]);

    for i_cond = 1:nconds
        
        subtightplot(1,2,i_cond,[0.02,0.02],[0.05,0.07],[0.05,0.05]);
        set(gca,'Color',[1 1 1]);
        temp = mean(mean(erp_diff_out(time_window,:,i_cond,:),4),1)';
        temp(16:18) = NaN;
        topoplot(temp,electrode_loc, 'whitebk','on','plotrad',.6,'maplimits',CLims,...
            'plotchans',electrode,'emarker',{'.','k',11,1})
        title(conds{i_cond});
        t = colorbar('peer',gca);
        set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
%         text( 2, 1, 'P3 (300-430ms) ')

        clear temp
    end
    
    % Overall subplot title
    supertitle(['Targets - Standards: ' num2str(itWin(1)) ' to ' num2str(itWin(2)) ' ms'],...
        'FontSize',10.5)
    
    clear itWin time_window

end


clear tWin erp_diff_out tw_i




