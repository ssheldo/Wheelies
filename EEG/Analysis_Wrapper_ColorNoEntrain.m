%% ANALYSIS WRAPPER

% ==== TRIGGERS ====
%
% -- Targets Present --
%Trial start: lag (1,2,3,4)
%Entrainers: 61-68
%Target: 20+lag (21,22,23,24)
%Colour Wheel: 40+lag (41,42,43,44)
%Mask: 90+lag (91,92,93,94)
%Response: 80+lag (81,82,83,84)
%
% -- Targets Not Present --
%Trial start: 10+lag (11,12,13,14)
%Entrainers: 61-68
%Target: 30+lag (31,32,33,34)
%Colour Wheel: 50+lag (51,52,53,54)%
%Mask: 95+lag (96,97,98,99)
%Response: 70+lag (71,72,73,74)



%clear and close everything
ccc

%%%%
% description of the dataset settings
% 1) byTargets: winsize is 256, no ERSP baseline. Epoched to targets
% 2) byTargets_v2: winsize is 128, no ERSP baseline. Epoched to targets
% 2) byTargets_v3: winsize is 256, no ERSP baseline. Epoched to targets;
%    postocularthresh = [-500 500]
% 3) byWheel: winsize is 256, no ERSP baseline. Epoched to color wheel
% 2) byTargets_v4: winsize is 256, no ERSP baseline. Epoched to targets;
%    postocularthresh = [-500 500]. ERP baseline [-600 -400] 
% 3) byWheel_v2: winsize is 256, no ERSP baseline. Epoched to color wheel;
%    postocularthresh = [-500 500]. ERP baseline [-1200 -800]
%%%%

%% Load data
exp.name = 'Colour_NoEntrain';
exp.conds = '';
% exp.participants = {'108';'109';'110';'111';'112';'113';'114'};
exp.participants = {'108';'109';'110';'112';'113';'114';'115';'116'};
% exp.participants = {'115';'116'};
exp.pathname = 'M:\Data\Colour\Colour_NoEntrain\';
exp.setname = {'byWheel_v2'}; % name each epoched set

%% Blink Correction
% the Blink Correction wants dissimilar events (different erps) seperated by 
% commas and similar events (similar erps) seperated with spaces. See 'help gratton_emcp'
% exp.selection_cards = {'11 21','13 23'};
%%%indicates where you want to center your data (where time zero is)
% exp.selection_cards = {'21 22 23 24'}; %must be list == length(exp.setname) 
exp.selection_cards = {'41 42 43 44'}; %must be list == length(exp.setname) 

%% Artifact rejection. 
% Choose the threshold to reject trials. More lenient threshold followed by an (optional) stricter threshold 
exp.preocularthresh = [-1000 1000]; %First happens before the ocular correction.
% exp.postocularthresh = [ ]; %Second happens after. Leave blank [] to skip
exp.postocularthresh = [-500 500]; %Second happens after. Leave blank [] to skip

%% Events and event labels
%%%for each condition ( lag 1-4 in this case), numbers correspond to
%%%triggers that will be kept for each condition. All other triggers will
%%%be removed
% exp.events = {[21,22,23,24]};%can be list or matrix (sets x events) 
exp.events = {[41,42,43,44]};%can be list or matrix (sets x events)   
exp.event_names = {'Wheel','Wheel','Wheel','Wheel'}; %must be list or matrix (sets x events)
exp.suffix = {'byWheel'};

%% Electrode location
%Where are your electrodes? (.ced file)
exp.electrode_locs = 'M:\Analysis\Colour\Colour_NoEntrain\EEG\electrode_map_RL_FCz.ced';
%electrode information
exp.electrode = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17];
exp.elec_names = {'O1';'P7';'T7';'P3';'C3';'F3';'Pz';'Cz';'Fz';'P4';'C4';'F4';'P8';'T8';'O2';'FCz'};

%% Re-referencing the data
exp.refelec = 16; %which electrode do you want to re-reference to?
exp.brainelecs = [1:15, 17]; %list of every electrode collecting brain data (exclude mastoid reference, EOGs, HR, EMG, etc.

%% Filter the data?
exp.filter = 'off'; %filter all files
exp.hicutoff = 40; %higher edge of the frequency pass band (Hz)
exp.locutoff = 0.1; %lower edge of the frequency pass band (Hz)

%% FFT/Wavelet settings
%how long is your window going to be? (Longer window == BETTER frequency 
% resolution & WORSE time resolution)
exp.winsize = 256; %use numbers that are 2^x, e.g. 2^10 == 1024ms

%baseline will be subtracted from the power variable. It is relative to 
% your window size. 
% exp.erspbaseline = [-200 0]; 
exp.erspbaseline = NaN;
%e.g., [-200 0] will use [-200-exp.winsize/2 0-exp.winsize/2]; 
% Can use just NaN for no baseline

%Instead of choosing a windowsize, you can choose a number of cycles per 
% frequency. See "help popnewtimef"
exp.cycles = [0]; %leave it at 0 to use a consistent time window

exp.freqrange = [1 40]; % what frequencies to consider? default is [1 50]

%% Epoching the data
exp.epoch = 'on'; %on to epoch data
%%%indicates where you want to center your data (where time zero is)
exp.epochs = {}; %must be list == length(exp.setname)
exp.epochs_name = {};
exp.epochslims = [-1.6 0.6]; %in seconds; epoched trigger is 0 e.g. [-1 2]
exp.epochbaseline = [-1200 -800]; %remove the baseline for each epoched set, in ms. e.g. [-200 0] 


%% Time-Frequency settings
%Do you want to run time-frequency analyses? (on/off)
exp.tf = 'on';
%Do you want to save the single-trial data? (on/off) (Memory intensive!!!)
exp.singletrials = 'on';
%Do you want to use all the electrodes or just a few? Leave blank [] for 
% all (will use same as exp.brainelecs)
exp.tfelecs = [];
%Saving the single trial data is memory intensive. Just use the electrodes
% you need. 
exp.singletrialselecs = [1 7 8 15 17]; %O1, Pz, Cz, O2, FCz
exp.singtrlelec_name = {'O1';'Pz';'Cz';'O2';'FCz'}; %because FCz is 17 

%//////////////////////////////////////////////////////////////////////////
%% Save your pipeline settings
% The settings will be saved as a new folder. It lets you save multiple datasets with different preprocessing parameters.
exp.settings = char(exp.setname); %name settings
% `````````````````````````````````````````````````````````````````````````
% Saving will help you remember what settings were used in each dataset
save([exp.settings '_Settings'],'exp') %save these settings as a .mat file. 
%//////////////////////////////////////////////////////////////////////////


%% Run preprocessing code
Preprocessing_ColourNoEntrain(exp)
%

%run analysis
% exp.electrode = 3
% Analysis_Att_Ent(exp)
%
