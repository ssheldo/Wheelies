



% =========================================================================
% /////////////////////////////////////////////////////////////////////////
% =========================================================================

% List of subject IDs
% parts = {'100';'101';'102';'103';'104';'105';'106';'107'};
parts = {'100';'102';'103';'104';'105';'106';'107'};
% parts = {'101';'102';'103';'104';'105';'106';'107';'108';'109';'110';'111'};

% Lags between last entrainer and target
lag = [5,10,15,20];
lag_times = [83.33/2; 83.33; 83.33+(83.33/2); 83.33*2];

% Get values for hue color
y_color = rgb2hsv([1,1,0]); %0.167
m_color = rgb2hsv([1,0,1]); %0.83
c_color = rgb2hsv([0,1,1]); %0.5
y_hue = y_color(1); %0.167
m_hue = m_color(1); %0.83
c_hue = c_color(1); %0.5

% =========================================================================
% /////////////////////////////////////////////////////////////////////////
% =========================================================================
for i_part = 6:length(parts) %loop through subjects
% for i_part = 3:5 %loop through subjects

    % Load each subject's data
    load(['M:\Experiments\Colour\Entrain_Colour_Wheel\Trials\Colour_NoEntrain_Data\' parts{i_part} '_Colour_NoEntrain.mat'])
%     load(['M:\Experiments\Colour\Entrain_Colour_Wheel\Trials\Colour_Entrain_Data\' parts{i_part} '_Colour_Entrain.mat'])

    % Get color of target
    for trialIndex = 1:prefs.nTrials
        data.target_colors(:,trialIndex) = prefs.colorwheel_trial{trialIndex}(data.presentedColorDegrees(trialIndex), :)';
    end
    
    % Remove practice trials (first 20 trials)
    error_deg = data.errorDegrees(1,21:end);
    error_rads = data.errorRads(1,21:end);
    lag_errors = data.lags(1,21:end);
    targets = data.target(1,21:end);
    color_report_deg = data.reportedColorDegrees(1,21:end);
    color_report_rads = data.reportedColorRads(1,21:end);
    color_present_deg = data.presentedColorDegrees(1,21:end);
    color_present_rads = data.presentedColorRads(1,21:end);
    color_diff_deg = data.error_differenceDegrees(1,21:end);
    color_diff_rads = data.error_differenceRads(1,21:end);
    target_colors = data.target_colors(:,21:end);
    
    % Remove trials where no target appeared
    error_deg = error_deg(targets == 1);
    error_rads = error_rads(targets == 1);
    lag_errors = lag_errors(targets == 1);
    color_report_deg = color_report_deg(targets == 1);
    color_report_rads = color_report_rads(targets == 1);
    color_present_deg = color_present_deg(targets == 1);
    color_present_rads = color_present_rads(targets == 1);
    color_diff_deg = color_diff_deg(targets == 1);
    color_diff_rads = color_diff_rads(targets == 1);
    target_colors = target_colors(:,targets == 1);
    
    %--------------------------------------------------
    % Split trials by color of target
    target_colors_hsv = rgb2hsv((target_colors./255)'); %change color map
    color_grp(1:length(target_colors_hsv)) = NaN; %pre-allocate

    for trl = 1:length(target_colors_hsv)

        hue = target_colors_hsv(trl,1); %get hue on trial

        if hue >= y_hue && hue < c_hue %yellow 2 cyan
            color_grp(1,trl) = 1; %green
        elseif hue >= c_hue && hue < m_hue %cyan 2 magenta
            color_grp(1,trl) = 2; %blue
        elseif hue >= m_hue || hue < y_hue %magenta 2 yellow
            color_grp(1,trl) = 3; %red
        end

    end
    %--------------------------------------------------
    
    % Split trials by target color
    error_deg_g = error_deg(color_grp == 1);
    error_deg_b = error_deg(color_grp == 2);
    error_deg_r = error_deg(color_grp == 3);
    
    % Save errors
    out_error_deg_g{i_part,1} = error_deg_g;
    out_error_deg_b{i_part,1} = error_deg_b;
    out_error_deg_r{i_part,1} = error_deg_r;
    
    
    % Get parameters from standard mixture model
    model = WithBias(StandardMixtureModel); %allows for central tendency to not be fixed at 0
    model_out{i_part,1} = MemFit(error_deg,model);
    model_out_g{i_part,1} = MemFit(error_deg_g,model);
    model_out_b{i_part,1} = MemFit(error_deg_b,model);
    model_out_r{i_part,1} = MemFit(error_deg_r,model);
%     model_out{i_part,1} = MLE(error_deg); %returns just parameters
   

    % Save output from model fit
    mu_out(i_part,1) = model_out{i_part,1}.maxPosterior(1);
    g_out(i_part,1) = model_out{i_part,1}.maxPosterior(2);
    sd_out(i_part,1) = model_out{i_part,1}.maxPosterior(3);
    
    mu_out_g(i_part,1) = model_out_g{i_part,1}.maxPosterior(1);
    g_out_g(i_part,1) = model_out_g{i_part,1}.maxPosterior(2);
    sd_out_g(i_part,1) = model_out_g{i_part,1}.maxPosterior(3);
    
    mu_out_b(i_part,1) = model_out_b{i_part,1}.maxPosterior(1);
    g_out_b(i_part,1) = model_out_b{i_part,1}.maxPosterior(2);
    sd_out_b(i_part,1) = model_out_b{i_part,1}.maxPosterior(3);
    
    mu_out_r(i_part,1) = model_out_r{i_part,1}.maxPosterior(1);
    g_out_r(i_part,1) = model_out_r{i_part,1}.maxPosterior(2);
    sd_out_r(i_part,1) = model_out_r{i_part,1}.maxPosterior(3);
%     g_out_phase(i_part,1) = model_out_phase{i_part,1}(1);
%     sd_out_phase(i_part,1) = model_out_phase{i_part,1}(2);

    clear error_deg error_rads lag_errors color_report_deg color_report_rads...
        color_present_deg color_present_rads color_diff_deg color_diff_rads...
        target_colors targets color_grp target_colors_hsv error_deg_g...
        error_deg_b error_deg_r

end




mean(g_in_phase)
mean(g_out_phase)

mean(sd_in_phase)
mean(sd_out_phase)



