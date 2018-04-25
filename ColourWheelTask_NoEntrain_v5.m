% ColourWheelTask_NoEntrain_v4() 
% Runs a color working memory task a la Zhang & Luck (2008). The task 
% requires memory for the color of briefly presented squares. Participants 
% then report the color of a single probed square using a contnuous report 
% task.
%
% Preferences can be found down at the bottom, beginning on line 197.
% 
% Fixation appears for one entrainer length instead of a blank screen.
%
% /////////////////////////////////////////////////////////////////////////
%
% ==== TRIGGERS ====
%
% -- Targets Present --
% Trial start: lag (1,2,3,4)
% Entrainers: 61-68
% Target: 20+lag (21,22,23,24)
% Colour Wheel: 40+lag (41,42,43,44)  
% Mask: 90+lag (91,92,93,94)
% Response: 80+lag (81,82,83,84)
%
% -- Targets Not Present --
% Trial start: 10+lag (11,12,13,14)
% Entrainers: 61-68
% Target: 30+lag (31,32,33,34)
% Colour Wheel: 50+lag (51,52,53,54)%
% Mask: 95+lag (96,97,98,99)
% Response: 70+lag (71,72,73,74)
% 
% /////////////////////////////////////////////////////////////////////////


function ColourWheelTask_NoEntrain_v5()

ccc

try
    prepareEnvironment;
    
    part_num = input('Participant Number:','s');
    Filename = ['M:\Experiments\Colour\Entrain_Colour_Wheel\Trials\Colour_NoEntrain_Data\' part_num '_Colour_NoEntrain_v4.mat'];
    
    
    window = openWindow();
    prefs = getPreferences();
    
    % counter of trials per block
    nblock = prefs.trials_per_block + 20; %including the first 20 practice trials
    
    
    % Get presentation timing information
    refresh = Screen('GetFlipInterval', window.onScreen); % Get flip refresh rate
    slack = refresh/2; % Divide by 2 to get slack
    
   
    % Get rects for each item.
    rects = cell(1, max(prefs.setSizes));
    for i = 1:max(prefs.setSizes)
        rects{i} = circularArrayRects([0, 0, prefs.squareSize, prefs.squareSize], ...
            i, prefs.radius, window.centerX, window.centerY)';
    end
    
    % Put up instructions and wait for keypress.
    instruct(window);
    
    % Location of color wheel on screen
    colorWheelLocations = colorwheelLocations(window,prefs);
    
    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %% Set up the random mix of lags and entrainers for each block
    %Now to find a lag you pick the next random index between 1:n_lags from i_lags
    %Then you find that index in all_lags, it tells you where to look in lags
    
    all_lags = [1:prefs.n_lags];
    if prefs.lags_per_block > 1
        for i_lag = 2:prefs.lags_per_block
            all_lags = [all_lags 1:prefs.n_lags];
        end
    end
    i_lags = randperm(prefs.lags_per_block * (prefs.n_lags));
    all_lags = all_lags(i_lags);
    
    
    %set up the catch trials on every n-lagsth trial
    p = 1/prefs.p_catchtrials;
    q = 1/prefs.p_catchtrials;
    present = [1];
    for i_pres = 2:prefs.lags_per_block * (prefs.n_lags)
        if i_pres == p
            present = [present 0];
            p = p + q;
        else
            present = [present 1];
        end
    end
    
    rand_pres = randperm(prefs.lags_per_block * (prefs.n_lags));
    present = present(rand_pres);
    
    trials = prefs.lags_per_block * (prefs.n_lags);
    
    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %% Set up the random rotation of the color wheel per trial
    for ii = 1:prefs.nTrials
        % Randomize which color to draw first = the color wheel appears to rotate
        irow = randi([1 360],1); %randomly select first row to start matrix
        trial_color_wheel = NaN(360,3); %pre-allocate variable
        % Fill new color matrix from irow to 360
        idx_row = irow;
        for i_color = 1:length(irow:360)
            trial_color_wheel(i_color,:) = prefs.colorwheel(idx_row,:);
            idx_row  = idx_row + 1;
        end
        clear i_color idx_row
        % Fill new color matrix from 1 to irow-1
        idx_row = length(irow:360) + 1;
        for i_color = 1:(irow - 1)
            trial_color_wheel(idx_row,:) = prefs.colorwheel(i_color,:);
            idx_row  = idx_row + 1;
        end
        clear i_color idx_row
        % Save new color matrix
        prefs.colorwheel_trial{ii} = trial_color_wheel;
        clear trial_color_wheel
    end
    clear ii
    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% START TASK
    for trialIndex = 1:length(prefs.fullFactorialDesign)
        
        lag = prefs.lags(all_lags(trialIndex))/5;
        
        % Determine how many items there are on this trial and the duration.
        nItems = prefs.setSizes(prefs.fullFactorialDesign(prefs.order(trialIndex), 1));
        retentionInterval = prefs.retentionIntervals; %(prefs.fullFactorialDesign(prefs.order(trialIndex), 2));
        
        % Pick an item to test.
%         itemToTest(trialIndex) = randsample(1:nItems); %legacy code that no longer works
        itemToTest(trialIndex) = nItems; %nItems should always be 1 because there will only be 1 target

        % Pick the colors for this trial.
        colorsInDegrees{trialIndex} = ceil(rand(1, nItems)*360);
        
        
        if present(trialIndex) == 1
            pause(window);
            
            %Present Fixation
            Screen('FillRect', window.onScreen, window.gray);
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FillRect',window.onScreen, Vpixx2Vamp(lag), prefs.trigger_size);
            t_fixate_onset = Screen('Flip', window.onScreen);
            
        elseif present(trialIndex) == 0
            pause(window);
            
            %Present Fixation
            Screen('FillRect', window.onScreen, window.gray);
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FillRect',window.onScreen, Vpixx2Vamp(10 + lag), prefs.trigger_size);
            t_fixate_onset = Screen('Flip', window.onScreen);
            
        end
        
        % Interval
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        tblank_onset = Screen('Flip', window.onScreen,t_fixate_onset + prefs.fixation_length*refresh - slack);
        

% =========================================================================       
        %% Entrainers
        if prefs.n_entrs > 0
%             Screen('FillOval', window.onScreen, (prefs.entr_grey+window.gray),...
%                 [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth,...
%                 rects{nItems}(4)+prefs.maskwidth]);
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
            Screen('FillRect',window.onScreen,Vpixx2Vamp(61),prefs.trigger_size);
            tentr_onset = Screen(window.onScreen, 'Flip', tblank_onset + prefs.preblank_length*refresh - slack);
            if prefs.n_entrs > 1
                for i_entr = 2:prefs.n_entrs
                    Screen('FillOval', window.onScreen, window.gray, rects{nItems});
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
                    tblank_onset = Screen(window.onScreen, 'Flip', tentr_onset + prefs.entr_length*refresh - slack);
                    Screen('FillOval', window.onScreen,  (prefs.entr_grey+window.gray),...
                        [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth,...
                        rects{nItems}(4)+prefs.maskwidth]);
                    Screen('FillOval', window.onScreen, window.gray, rects{nItems});
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(60 + i_entr),prefs.trigger_size);
                    tentr_onset = Screen(window.onScreen, 'Flip', tblank_onset + prefs.entr_gap_length*refresh - slack);
                end
            end
        end
        
        if present(trialIndex) == 1 %two options depending on whether the target is present or absent
% =========================================================================             
            %% Lag
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
%             Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            tlag_onset = Screen(window.onScreen, 'Flip', tentr_onset + prefs.entr_length*refresh - slack);
            
            % Draw the stimulus.
%             colorsToDisplay = prefs.colorwheel(colorsInDegrees{trialIndex}, :)';
            colorsToDisplay = prefs.colorwheel_trial{trialIndex}(colorsInDegrees{trialIndex}, :)';
            prefs.colorsToDisplay{trialIndex} = colorsToDisplay; %save color
            Screen('FillOval', window.onScreen, colorsToDisplay, rects{nItems});
            Screen('FillRect',window.onScreen,Vpixx2Vamp(20 + lag),prefs.trigger_size);
            
            % Post the stimulus
            ttarget_onset = Screen('Flip', window.onScreen, tlag_onset + prefs.lagISI(all_lags(trialIndex))*refresh - slack);
            
            %Interval
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            tISI_onset = Screen('Flip', window.onScreen, ttarget_onset + prefs.targ_length*refresh - slack);
            
        else
            %% Lag
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
%             Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            tlag_onset = Screen(window.onScreen, 'Flip', tentr_onset + prefs.entr_length*refresh - slack);
            
            %% present the Missing Target
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('FillRect',window.onScreen,Vpixx2Vamp(30 + lag),prefs.trigger_size);
            ttarget_onset = Screen(window.onScreen, 'Flip', tlag_onset + prefs.lagISI(all_lags(trialIndex))*refresh - slack);
            
            WaitSecs(.01);
            
            %% blank Inter stimulus interval
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            tISI_onset = Screen(window.onScreen, 'Flip', ttarget_onset + prefs.targ_length*refresh - slack);
            
        end
% =========================================================================         
        %% Mask
        
        %Create uniform RGB planes to make the texture       
%         noiseimg = reshape(prefs.colorwheel(ceil(360*rand(prefs.maskwidth/prefs.tilesize, prefs.maskwidth/prefs.tilesize)),:),[floor(prefs.maskwidth/prefs.tilesize), floor(prefs.maskwidth/prefs.tilesize),3]);
        noiseimg = reshape(prefs.colorwheel_trial{trialIndex}(ceil(360*rand(prefs.maskwidth/prefs.tilesize, prefs.maskwidth/prefs.tilesize)),:),[floor(prefs.maskwidth/prefs.tilesize), floor(prefs.maskwidth/prefs.tilesize),3]);

        % Convert it to a texture 'tex':
        tex=Screen('MakeTexture', window.onScreen, noiseimg);
 
        % Draw a mask-sized aperture to view the texture
        aperture =Screen('OpenOffscreenwindow', window.onScreen, 128);     
        Screen('FillOval', aperture, [255 255 255 0], [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth, rects{nItems}(4)+prefs.maskwidth]);
        Screen('FillOval', aperture, [128 128 128 255], rects{nItems});
    
        Screen('BlendFunction', window.onScreen, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        % Draw the texture 
        if present(trialIndex) == 1
            Screen('DrawTexture', window.onScreen, tex, [], [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth, rects{nItems}(4)+prefs.maskwidth], [], 0);
            Screen('DrawTexture', window.onScreen, aperture, [], [], [], 0); 
            Screen('FillRect',window.onScreen,Vpixx2Vamp(90 + lag),prefs.trigger_size);
        elseif present(trialIndex) == 0
            Screen('DrawTexture', window.onScreen, tex, [], [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth, rects{nItems}(4)+prefs.maskwidth], [], 0);
            Screen('DrawTexture', window.onScreen, aperture, [], [], [], 0); 
            Screen('FillRect',window.onScreen,Vpixx2Vamp(95 + lag),prefs.trigger_size);
        end
        
        t_maskonset = Screen('Flip', window.onScreen, tISI_onset + prefs.maskISI*refresh - slack);
        
        % Interval
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        Screen('Flip', window.onScreen, t_maskonset + prefs.mask_length*refresh - slack);
        
        % Retention interval
        WaitSecs(retentionInterval);
        
        % Choose a circle to test, then display the response screen.
        data.presentedColorRads(trialIndex) = deg2rad(colorsInDegrees{trialIndex}(itemToTest(trialIndex)));
        data.presentedColorDegrees(trialIndex) = colorsInDegrees{trialIndex}(itemToTest(trialIndex));
        colorsOfTest = repmat([120 120 120], nItems, 1);
        colorsOfTest(itemToTest(trialIndex), :) = [145 145 145];
        
        if present(trialIndex) == 1
            Screen('FillRect',window.onScreen, Vpixx2Vamp(40 + lag), prefs.trigger_size);
            Screen('Flip', window.onScreen);
        else
            Screen('FillRect',window.onScreen, Vpixx2Vamp(50 + lag), prefs.trigger_size);
            Screen('Flip', window.onScreen);
        end

% ========================================================================= 
        %% Response
        
        % Draw color wheel
        drawColorWheel(window, prefs, prefs.colorwheel_trial{trialIndex});
        
        % Wait for click.
        SetMouse(window.centerX, window.centerY);
        ShowCursor('Arrow');
        
        % If mouse button is already down, wait for release.
        GetMouse(window.onScreen);
        buttons = 0;
        while any(buttons)
            [x, y, buttons] = GetMouse(window.onScreen);
        end
        
        everMovedFromCenter = false;
        while ~any(buttons)
            
            drawColorWheel(window, prefs, prefs.colorwheel_trial{trialIndex});
            
            [x,y,buttons] = GetMouse(window.onScreen);
            [minDistance, minDistanceIndex] = min(sqrt((colorWheelLocations(1, :) - x).^2 + (colorWheelLocations(2, :) - y).^2));
            
            if(minDistance < 500)
                everMovedFromCenter = true;
            end
            
            if(everMovedFromCenter)
%                 colorsOfTest(itemToTest(trialIndex), :) = prefs.colorwheel(minDistanceIndex,:);
                colorsOfTest(itemToTest(trialIndex), :) = prefs.colorwheel_trial{trialIndex}(minDistanceIndex,:);
            else
                colorsOfTest(itemToTest(trialIndex), :) = [145 145 145];
            end
            
            Screen('FillOval', window.onScreen, colorsOfTest', rects{nItems});
            
                        drawColorWheel(window, prefs, prefs.colorwheel_trial{trialIndex});
                        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
                        Screen('Flip', window.onScreen);
            
        end
        
        
            if present(trialIndex) == 1
                Screen('FillRect',window.onScreen,Vpixx2Vamp(80 + lag),prefs.trigger_size);
                Screen('Flip', window.onScreen);
                
            elseif present(trialIndex) == 0
                Screen('FillRect',window.onScreen,Vpixx2Vamp(70 + lag),prefs.trigger_size);
                Screen('Flip', window.onScreen);
            end
        
        data.reportedColorRads(trialIndex) = deg2rad(minDistanceIndex);
        data.reportedColorDegrees(trialIndex) = minDistanceIndex;
        
        HideCursor
        
        if trialIndex == length(prefs.fullFactorialDesign)
            finish(window);
        elseif trialIndex == nblock
            rest(window);
            nblock = nblock + prefs.trials_per_block;
        elseif trialIndex == 20 % first 20 trials are practice trials
            practice(window);
        end
        
    end
    
% =========================================================================     
    %% Preliminary analysis of results.
    %errors
    data.errorDegrees = (180/pi) .* (angle(exp(1i*data.reportedColorRads)./exp(1i*data.presentedColorRads)));
    data.errorRads = deg2rad(data.errorDegrees);
    data.error_differenceRads = data.reportedColorRads - data.presentedColorRads;
    data.error_differenceDegrees = data.reportedColorDegrees - data.presentedColorDegrees;
    %lags
    data.lags = all_lags;
    %target present or absent
    data.target = present;
    
    save(Filename, 'data', 'prefs');
    postpareEnvironment;

% =========================================================================     
% =========================================================================     
catch
    postpareEnvironment;
    psychrethrow(psychlasterror);
    
end % end try/catch

% ========================================================================= 
% ========================================================================= 
end % end whole colorworkingmemoryscript

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ========================================================================= 
% -------------------------------------------------------------------------
% #########################################################################
% #########################################################################
% -------------------------------------------------------------------------
% ========================================================================= 


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function prepareEnvironment

clear all;

% HideCursor; % Comment out when debugging

commandwindow; % Select the command window to avoid typing in open scripts

% Seed the random number generator.
% RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100*clock)));
RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', sum(100*clock)));

%   ListenChar(2); % Don't print to MATLAB command window
% Screen('Preference', 'SkipSyncTests', 1);

% /////////////////////////////////////////////////////////////////////////
%% Set up parallel port
%initialize the inpoutx64 low-level I/O driver
% config_io;
% %optional step: verify that the inpoutx64 driver was successfully installed
% global cogent;
% if( cogent.io.status ~= 0 )
%     error('inp/outp installation failed');
% end
% %write a value to the default LPT1 printer output port (at 0x378)
% address_eeg = hex2dec('B010');
% 
% outp(address_eeg,0);  %set pins to zero
% /////////////////////////////////////////////////////////////////////////
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function postpareEnvironment
ShowCursor;
%   ListenChar(0);
Screen('CloseAll');
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function instruct(window)

%Screen 1
prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
% Screen('DrawText', window.onScreen, 'The next few screens will briefly explain the task.', 100, 560, 255);
Screen('DrawText', window.onScreen, 'The next few screens will briefly explain the task.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'After the instructions, you will perform several practice trials followed by the experiment.',(window.centerX-400),(window.centerY+50), 255); %line every +30
Screen('DrawText', window.onScreen, 'On each screen click the mouse to continue.',(window.centerX-400),(window.centerY+80), 255); %line every +25
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);

%Screen 2
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'A fixation dot will appear.  Please remain focused on the white central dot during the entire task.',(window.centerX-600),(window.centerY+30), 255);
Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %255 is color (white)
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 4
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'After a few moments the dot will disappear and a target consisting of a coloured circle will then quickly appear.',...
    (window.centerX-600),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'It might be difficult to see the target, but you should try to detect what colour it is.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
Screen('FillOval', window.onScreen, [0,177.5,103.5417], [950,530,970,550]);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 5
Screen('TextSize', window.onScreen, window.fontsize);
% Screen('DrawText', window.onScreen, 'A multicoloured donut will then quickly appear and then disappear.',...
%     (window.centerX-600),(window.centerY+50), 255);
% Screen('DrawText', window.onScreen, 'Remember that you want to detect the colour of the circle, not the donut.',...
%     (window.centerX-400),(window.centerY+80), 255); %line every +30
noiseimg = reshape(prefs.colorwheel(ceil(360*rand(prefs.maskwidth/prefs.tilesize, prefs.maskwidth/prefs.tilesize)),:),[floor(prefs.maskwidth/prefs.tilesize), floor(prefs.maskwidth/prefs.tilesize),3]);
% Convert it to a texture 'tex':
tex=Screen('MakeTexture', window.onScreen, noiseimg);
% Draw a mask-sized aperture to view the texture
aperture =Screen('OpenOffscreenwindow', window.onScreen, 128);   
rects{1} = [950;530;970;550];
Screen('FillOval', aperture, [255 255 255 0], [rects{1}(1)-prefs.maskwidth, rects{1}(2)-prefs.maskwidth,rects{1}(3)+prefs.maskwidth, rects{1}(4)+prefs.maskwidth]);
% Screen('FillOval', aperture, [128 128 128 255], rects{1});
Screen('BlendFunction', window.onScreen, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% Draw the texture 
Screen('DrawTexture', window.onScreen, tex, [], [rects{1}(1)-prefs.maskwidth, rects{1}(2)-prefs.maskwidth,rects{1}(3)+prefs.maskwidth, rects{1}(4)+prefs.maskwidth], [], 0);
Screen('DrawTexture', window.onScreen, aperture, [], [], [], 0); 
Screen('DrawText', window.onScreen, 'A multicoloured donut will then quickly appear and then disappear.',...
    (window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Remember that you want to detect the colour of the circle, not the donut.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
Screen('FillOval', window.onScreen, window.gray, rects{1});
Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 6
colorWheelLocations = [cosd(1:360).*prefs.colorWheelRadius + window.centerX; ...
    sind(1:360).*prefs.colorWheelRadius + window.centerY];
colorWheelSizes = 20;
Screen('DrawDots', window.onScreen, colorWheelLocations, colorWheelSizes, prefs.colorwheel', [], 1);
Screen('FillOval', window.onScreen, [157.5,157.5,157.5], [950,530,970,550]);
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'A colour wheel will then appear.',(window.centerX-880),(window.centerY-50), 255);
Screen('DrawText', window.onScreen, 'Using the mouse, click on the part of the',(window.centerX-880),(window.centerY), 255);
Screen('DrawText', window.onScreen, 'colour wheel you believe represents the',(window.centerX-880),(window.centerY+30), 255);
Screen('DrawText', window.onScreen, 'colour of the target that appeared.',(window.centerX-880),(window.centerY+60), 255);
Screen('DrawText', window.onScreen, 'If you are unsure of the colour of the target,',(window.centerX-880),(window.centerY+110), 255);
Screen('DrawText', window.onScreen, 'or if you did not see a target, please provide',(window.centerX-880),(window.centerY+140), 255);
Screen('DrawText', window.onScreen, 'your best guess.',(window.centerX-880),(window.centerY+170), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);

%Screen 7
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You will now complete several practice trials so that you become more familiar with the task.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Remember to remain focused and fixated on the white central dot that will appear.',(window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Click the left mouse button when you are ready to begin the practice trials.',(window.centerX-400),(window.centerY+80), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function pause(window)
prefs = getPreferences();
Screen('FillRect', window.onScreen, window.gray);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(0.5);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function practice(window)
prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You have now finished the set of practice trials.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Please let the experimenter know by using the call box.',(window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Do NOT begin the experiment.',(window.centerX-400),(window.centerY+80), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function rest(window)
prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'Feel free to take a break at this time.  Click the mouse when you are ready to continue.',(window.centerX-400),(window.centerY+20), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function finish(window)
prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You are now finished the experiment. Thank you for your time.',(window.centerX-400),(window.centerY+20), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

% function drawFixation(window, fixationX, fixationY, fixationSize)
% prefs = getPreferences();
% Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
% Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
% end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function offsets = circularArrayOffsets(n, centerX, centerY, radius, rotation)
degreeStep = 360/n;
offsets = [sind(0:degreeStep:(360-degreeStep) + rotation)'.* radius, ...
    cosd(0:degreeStep:(360-degreeStep) + rotation)'.* radius];
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function rects = circularArrayRects(rect, nItems, radius, centerX, centerY)
coor = circularArrayOffsets(nItems, centerX, centerY, radius, 0) + repmat([centerX centerY], nItems, 1);
rects = [coor(:, 1)-rect(3)/2 , coor(:, 2)-rect(3)/2, coor(:, 1)+rect(3)/2, coor(:, 2)+rect(3)/2];
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

% function returnToFixation(window, fixationX, fixationY, fixationSize)
% prefs = getPreferences();
% Screen('FillRect', window.onScreen, window.gray);
% Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
% Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
% Screen('Flip', window.onScreen);
% end
%
% function returnToFixation_pres(window, fixationX, fixationY, fixationSize)
% prefs = getPreferences();
% Screen('FillRect', window.onScreen, window.gray);
% Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
% Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
% t_fixate_onset = Screen('Flip', window.onScreen);
% end
%
% function returnToFixation_abs(window, fixationX, fixationY, fixationSize)
% prefs = getPreferences();
% Screen('FillRect', window.onScreen, window.gray);
% Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
% Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
% t_fixate_onset = Screen('Flip', window.onScreen);
% end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% ~~~ Open the main window and get dimensions ~~~
function window = openWindow()

window.screenNumber = max(Screen('Screens'));
window.onScreen = Screen('OpenWindow', window.screenNumber, [127.5 127.5 127.5]);
Screen('BlendFunction', window.onScreen, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[window.screenX, window.screenY] = Screen('WindowSize', window.onScreen); % check resolution
window.screenRect  = [0, 0, window.screenX, window.screenY]; % screen rect
window.centerX = window.screenX * 0.5; % center of screen in X direction
window.centerY = window.screenY * 0.5; % center of screen in Y direction
window.centerXL = floor(mean([0, window.centerX])); % center of left half of screen in X direction
window.centerXR = floor(mean([window.centerX, window.screenX])); % center of right half of screen in X direction

% Basic drawing and screen variables.
window.black    = BlackIndex(window.onScreen);
window.white    = WhiteIndex(window.onScreen);
window.gray     = mean([window.black window.white]);
window.fontsize = 26; % size of instruction text
window.bcolor   = window.gray;

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% Draw the color wheel
function drawColorWheel(window, prefs, trial_color_wheel)
prefs = getPreferences();
colorWheelLocations = [cosd(1:360).* prefs.colorWheelRadius + window.centerX; ...
    sind(1:360).* prefs.colorWheelRadius + window.centerY];
colorWheelSizes = 20;
% colorWheelSizes = 60;

% Now draws "rotated" color wheel
Screen('DrawDots', window.onScreen, colorWheelLocations, colorWheelSizes, trial_color_wheel', [], 1);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
function L = colorwheelLocations(window,prefs)
L = [cosd(1:360).*prefs.colorWheelRadius + window.centerX; ...
    sind(1:360).*prefs.colorWheelRadius + window.centerY];
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% ~~~ Set task variables for timing and stimuli ~~~
function prefs = getPreferences

prefs.retentionIntervals = [0.50]; %time between target and color wheel (I think)

% ---------------
% Target --------
% ---------------
prefs.targ_length = 1; %in refresh cycles; each refresh is 1000msec / 120 Hz = 8.333 msec
prefs.setSizes = [1]; %number of targets presented (should always be 1)
prefs.squareSize = 20; % size of each stimulus object, in pixels
prefs.radius = 0; %how far apart the targets are from each other (if more than 1 target)


% ---------------
% Entrainers ----
% ---------------
prefs.entr_grey = 0; %points darker than background (do not want to see entrainers)
prefs.entr_length = 1 ; %refresh cycles of entrainer = refresh cycles of Draw.entrainer
prefs.n_entrs = [1]; %number of entrainers on a trial (can be list: entrs = [6:1:8]);
                     %number of entrainers set to 1 because 0 will break the code


% ---------------
% Lags ----------
% ---------------
prefs.lags_per_block = 105; %how many of each lag in each block 105

prefs.rate = 10; %refresh cycles before next entrainer (1000msec /120 Hz) = 8.333 msec; 1000msec / 10 Hz  = 100 msec; 100 msec / 8.333 msec = 12 cycles
                          %1000/8.33*6 == 20Hz, 1000/8.33*8 == 15Hz, 1000/8.33*10 == 12Hz, 1000/8.33*14 == 8.5Hz, 1000/8.33*30 == 4.0Hz
                          %number of refreshes, so formula = 1000/(8.33*desired frequency) desired frequency = 12, 15, 20, etc. 

%number of refreshes after last entrainer before target (can be list: lags = [4:2:8])
prefs.lags = [prefs.rate/2:prefs.rate/2:prefs.rate*2]; %lag rate*8.333 = time ms
prefs.n_lags = length(prefs.lags); %number of unique lags
prefs.lagISI = prefs.lags - 1;


% ---------------
% Mask ----------
% ---------------
prefs.SOA = 6; %refresh target onset to mask onset (50 ms optimal/8.3333 msec = 6 cycles)
prefs.maskISI = prefs.SOA - 1; %target OFFSET to mask onset
prefs.maskwidth = 20;
prefs.mask_length = 1; %refresh cycles of mask
prefs.tilesize = 5; % how big the coloured squares are


% Variables
prefs.trigger_size = [0 0 1 1];
prefs.p_catchtrials = 0.2; %what proportion of trials will be catch trials
prefs.entr_gap_length = prefs.rate - prefs.entr_length;

% ---------------
% Fixation ------
% ---------------
prefs.fixation_length = 60; %500ms
prefs.preblank_length = 24; %200ms
prefs.fixationSize = 4; %size of dot


% ---------------
% Color Wheel ---
% ---------------
prefs.colorWheelRadius = 360; %size of color wheel

prefs.colorwheel_path = 'M:\Experiments\Colour\Entrain_Colour_Wheel\Colour Wheels\';
% prefs.colorwheel = load('colorwheel360.mat', 'fullcolormatrix');
% prefs.colorwheel = prefs.colorwheel.fullcolormatrix;
% prefs.colorwheel = load([prefs.colorwheel_path 'newest_colorwheel.mat'], 'newest_colorwheel');
% prefs.colorwheel = prefs.colorwheel.newest_colorwheel;
% prefs.colorwheel = load([prefs.colorwheel_path 'lum_colorwheel.mat'], 'lum_colorwheel');
% prefs.colorwheel = prefs.colorwheel.lum_colorwheel;
prefs.colorwheel = load([prefs.colorwheel_path 'bright_colour_wheel.mat'], 'bright_colour_wheel');
prefs.colorwheel = prefs.colorwheel.bright_colour_wheel;


% -------------------
% Trials & Blocks ---
% -------------------
% Randomize trial order of full factorial design order.
prefs.fullFactorialDesign = fullfact([length(prefs.setSizes), ...
    prefs.lags_per_block, ...
    prefs.n_lags]);
prefs.order = Shuffle(1:length(prefs.fullFactorialDesign));

% Total number of trials
prefs.nTrials = length(prefs.fullFactorialDesign);
% Number of blocks - must me a multiple of prefs.nTrials-20 (8 blocks when 400 trials)
prefs.nBlocks = 8;
% Trials per block
prefs.trials_per_block = (prefs.nTrials - 20)/prefs.nBlocks;

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
