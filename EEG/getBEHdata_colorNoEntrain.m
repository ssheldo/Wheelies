function [error_deg,color_grp] = getBEHdata_colorNoEntrain(subj)




% List of subject IDs
% parts = {'108';'109';'110';'111';'112';'113'};

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


% Load each subject's data
load(['M:\Experiments\Colour\Entrain_Colour_Wheel\Trials\Colour_NoEntrain_Data\' subj '_Colour_NoEntrain_v4.mat'])
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
% out_error_deg_g{i_part,1} = error_deg_g;
% out_error_deg_b{i_part,1} = error_deg_b;
% out_error_deg_r{i_part,1} = error_deg_r;






