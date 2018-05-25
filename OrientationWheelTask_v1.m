% OrientationWheelTask_v1() 
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


function OrientationWheelTask_v1()

ccc

try
    
    prepareEnvironment;
% -------------------------------------------------------------------------    
    part_num = input('Participant Number:','s');
    Filename = ['M:\Experiments\OrientationWheel\Orient_Data\' part_num '_Orient_v1.mat'];
% -------------------------------------------------------------------------    
    window = openWindow(); %get info for psychtoolbox
    prefs = getPreferences(); %get task variable info
% -------------------------------------------------------------------------    
    maskwin = createMasktex(window,prefs); %create offscreen window with mask
% ------------------------------------------------------------------------- 

    % counter of trials per block
    nblock = prefs.trials_per_block + 20; %including the first 20 practice trials
    
    % total number of trials including practice trials
    ntotaltrial = prefs.nTrials + 20; %add 20 practice trials
    
    % Get presentation timing information
    refresh = Screen('GetFlipInterval', window.onScreen); % Get flip refresh rate
    slack = refresh/2; % Divide by 2 to get slack
    
   
    % Get rects for each item.
    rects = cell(1, max(prefs.setSizes));
    for i = 1:max(prefs.setSizes)
        rects{i} = circularArrayRects([0, 0, prefs.squareSize, prefs.squareSize], ...
            i, prefs.radius, window.centerX, window.centerY)';
    end
    
    % Center the target oval on the centre of the screen (for drawing
    % target stimuli)
    centeredRectresp = CenterRectOnPointd(prefs.baseCircleresp,window.centerX,window.centerY); 
            
% -------------------------------------------------------------------------    
    % Put up instructions and wait for keypress.
    instruct(window);
% -------------------------------------------------------------------------    
    % Location of color wheel on screen
    orientWheelLocations = orientwheelLocations(window,prefs);
% -------------------------------------------------------------------------    
    
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
    
    
    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% START TASK
    for trialIndex = 1:ntotaltrial
        
        lag = prefs.lags(all_lags(trialIndex))/5;
        
        % Determine how many items there are on this trial and the duration (this is left over code from original task)
        nItems = 1; %should always be 1 because there will only be 1 target
        retentionInterval = prefs.retentionIntervals; %(prefs.fullFactorialDesign(prefs.order(trialIndex), 2));
        
        % Pick an item to test.
%         itemToTest(trialIndex) = randsample(1:nItems); %legacy code that no longer works
        itemToTest(trialIndex) = nItems; %nItems should always be 1 because there will only be 1 target

        % Pick the orientation for this trial.
        trialOrient{trialIndex} = prefs.degslocs(prefs.selectorientTrial(trialIndex));
        
        
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
% =========================================================================            
            %% Draw the target stimulus
            
            yadj = prefs.orientwheel(2,trialOrient{trialIndex}); %adjust lines y position 
            xadj = prefs.orientwheel(1,trialOrient{trialIndex}); %adjust lines x position 
            yadj2 = prefs.orientwheel2(2,trialOrient{trialIndex}); %adjust lines y position 
            xadj2 = prefs.orientwheel2(1,trialOrient{trialIndex}); %adjust lines x position 
            % Set the coordinates (these are all relative to zero we will let
            % the drawing routine center it for us)
            xCoords = [xadj -xadj2 xadj xadj2 xadj 0];
            yCoords = [yadj yadj2 yadj -yadj2 yadj 0];
            allCoords_targ{trialIndex} = [xCoords; yCoords];
            % Draw the orientation lines, set it to the center of our screen and
            % set good quality antialiasing
            Screen('DrawLines', window.onScreen, allCoords_targ{trialIndex}, prefs.lineWidthPix, window.black, [window.centerX window.centerY], 2);
            % Draw the central circle to the screen
            Screen('FillOval', window.onScreen, window.black, centeredRectresp');

            Screen('FillRect',window.onScreen,Vpixx2Vamp(20 + lag),prefs.trigger_size);
           
            % Post the stimulus
            ttarget_onset = Screen('Flip', window.onScreen, tlag_onset + prefs.lagISI(all_lags(trialIndex))*refresh - slack);
            
            % Interval
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            tISI_onset = Screen('Flip', window.onScreen, ttarget_onset + prefs.targ_length*refresh - slack);
            
% =========================================================================            
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
        % Draw the texture 
        if present(trialIndex) == 1
            Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
            Screen('FillRect',window.onScreen,Vpixx2Vamp(90 + lag),prefs.trigger_size);
        elseif present(trialIndex) == 0
            Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
            Screen('FillRect',window.onScreen,Vpixx2Vamp(95 + lag),prefs.trigger_size);
        end
        
        t_maskonset = Screen('Flip', window.onScreen, tISI_onset + prefs.maskISI*refresh - slack);
        
        % Interval
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        Screen('Flip', window.onScreen, t_maskonset + prefs.mask_length*refresh - slack);
        
        % Retention interval
        WaitSecs(retentionInterval);
        
% ========================================================================= 
        %% Response
        
        % Choose a circle to test, then display the response screen.
        data.presentedOrientRads(trialIndex) = deg2rad(trialOrient{trialIndex});
        data.presentedOrientDegrees(trialIndex) = trialOrient{trialIndex};
%         colorsOfTest = repmat([120 120 120], nItems, 1);
%         colorsOfTest(itemToTest(trialIndex), :) = [145 145 145];
        
        if present(trialIndex) == 1
            Screen('FillRect',window.onScreen, Vpixx2Vamp(40 + lag), prefs.trigger_size);
            Screen('Flip', window.onScreen);
        else
            Screen('FillRect',window.onScreen, Vpixx2Vamp(50 + lag), prefs.trigger_size);
            Screen('Flip', window.onScreen);
        end
        
        %------------------------------------------------------------------
        % Wait for click.
        SetMouse(window.centerX, window.centerY); %mouse at center
        ShowCursor('Arrow');
        %------------------------------------------------------------------
        % If mouse button is already down, wait for release.
        GetMouse(window.onScreen);
        buttons = 0;
        while any(buttons)
            [x, y, buttons] = GetMouse(window.onScreen);
        end
        %------------------------------------------------------------------
        everMovedFromCenter = false;
        while ~any(buttons) % keep track of mouse location if moved
            
            [x,y,buttons] = GetMouse(window.onScreen);
            [minDistance, minDistanceIndex] = min(sqrt((orientWheelLocations(1, :) - x).^2 + (orientWheelLocations(2, :) - y).^2));
            
            if(minDistance < 3) %make sure subject moved mouse
                everMovedFromCenter = true;
            end
            
            if(everMovedFromCenter) %as the mouse is moved, update coords of response stimulus to match
                % new coordinates of start and end of lines based on mouse position
                yloc = prefs.orientwheel(2,minDistanceIndex); %adjust lines y position 
                xloc = prefs.orientwheel(1,minDistanceIndex); %adjust lines x position 
                yloc2 = prefs.orientwheel2(2,minDistanceIndex); %for edge lines 
                xloc2 = prefs.orientwheel2(1,minDistanceIndex); %for edge lines 
            else
                %starting location of lines if mouse is not moved
                yloc = -30;
                xloc = 0;
                yloc2 = 0;
                xloc2 = 2;
            end
            
            % Set the new coordinates
            xCoords = [xloc -xloc2 xloc xloc2 xloc 0];
            yCoords = [yloc yloc2 yloc -yloc2 yloc 0];
            allCoords = [xCoords; yCoords];

            % Draw the orient lines, set it to the center of our screen and
            % set good quality antialiasing
            Screen('DrawLines', window.onScreen, allCoords, prefs.lineWidthPix, window.black, [window.centerX window.centerY], 2);
            % Draw the rect to the screen
            Screen('FillOval', window.onScreen, window.black, centeredRectresp');
            
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            Screen('Flip', window.onScreen);
            
            clear xCoords yCoords allCoords
            
        end
        %------------------------------------------------------------------
        
        if present(trialIndex) == 1
            Screen('FillRect',window.onScreen,Vpixx2Vamp(80 + lag),prefs.trigger_size);
            Screen('Flip', window.onScreen);

        elseif present(trialIndex) == 0
            Screen('FillRect',window.onScreen,Vpixx2Vamp(70 + lag),prefs.trigger_size);
            Screen('Flip', window.onScreen);
        end
        
        [respDistance, respDistanceIndex] = min(sqrt((orientWheelLocations(1,:) - x).^2 + (orientWheelLocations(2,:) - y).^2));
        data.reportedOrientRads(trialIndex) = deg2rad(respDistanceIndex);
        data.reportedOrientDegrees(trialIndex) = respDistanceIndex;
        data.reportedRespDistance(trialIndex) = respDistance;
        
        HideCursor
        
% =========================================================================        
        if trialIndex == ntotaltrial
            finish(window);
        elseif trialIndex == nblock
            rest(window);
            nblock = nblock + prefs.trials_per_block;
        elseif trialIndex == 20 % first 20 trials are practice trials
            practice(window);
        end
% =========================================================================   
    end %end trials
    
% =========================================================================     
    %% Preliminary analysis of results.
    %errors
    data.errorDegrees = (180/pi) .* (angle(exp(1i*data.reportedOrientRads)./exp(1i*data.presentedOrientRads)));
    data.errorRads = deg2rad(data.errorDegrees);
    data.error_differenceRads = data.reportedOrientRads - data.presentedOrientRads;
    data.error_differenceDegrees = data.reportedOrientDegrees - data.presentedOrientDegrees;
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
end % end whole script

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

ccc

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


%Screen 3
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'After a few moments the dot will disappear and a target pointing in a certain direction will quickly appear.',...
    (window.centerX-600),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'It might be difficult to see the target, but you should try to detect which way it is pointing.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
%draw target stimulus
xCoords = [15 1.73 15 -1.73 15 0];
yCoords = [-25.98 1 -25.98 -1 -25.98 0];
allCoords_targ = [xCoords; yCoords];
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, window.black, [window.centerX window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerX, window.centerY); 
Screen('FillOval', window.onScreen, window.black, centeredRectresp');
% Screen('FillOval', window.onScreen, [0,177.5,103.5417], [950,530,970,550]);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 4
Screen('TextSize', window.onScreen, window.fontsize);
maskwin = createMasktex(window,prefs); %create offscreen window with mask
Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
Screen('DrawText', window.onScreen, 'A star-like shape will then quickly appear and then disappear.',...
    (window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Remember that you want to detect the orientation of the target, not the star.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
% Screen('FillOval', window.onScreen, window.gray, rects{1});
Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 5
%draw target stimulus
xCoords = [-21.21 -1.41 -21.21 1.41 -21.21 0];
yCoords = [21.21 -1.41 21.21 1.41 21.21 0];
allCoords_targ = [xCoords; yCoords];
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, window.black, [window.centerX window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerX, window.centerY); 
Screen('FillOval', window.onScreen, window.black, centeredRectresp');
% Screen('FillOval', window.onScreen, [157.5,157.5,157.5], [950,530,970,550]);
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'Another target will then appear.',(window.centerX-880),(window.centerY-50), 255);
Screen('DrawText', window.onScreen, 'Using the mouse, move the new target until',(window.centerX-880),(window.centerY), 255);
Screen('DrawText', window.onScreen, 'it is pointing in the exact same direction',(window.centerX-880),(window.centerY+30), 255);
Screen('DrawText', window.onScreen, 'of the original target that appeared.',(window.centerX-880),(window.centerY+60), 255);
Screen('DrawText', window.onScreen, 'If you are unsure of the direction of the target,',(window.centerX-880),(window.centerY+110), 255);
Screen('DrawText', window.onScreen, 'or if you did not see a target, please provide',(window.centerX-880),(window.centerY+140), 255);
Screen('DrawText', window.onScreen, 'your best guess.',(window.centerX-880),(window.centerY+170), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 6
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
function drawOrientWheel(window)
prefs = getPreferences();
orientWheelLocations = [sind(1:360).*prefs.orientWheelRadius + window.centerX;...
    cosd(1:360).*prefs.orientWheelRadius + window.centerY];
% colorWheelSizes = 60;

% Now draws "rotated" color wheel
% Screen('DrawDots', window.onScreen, orientWheelLocations, colorWheelSizes, trial_color_wheel', [], 1);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
function L = orientwheelLocations(window,prefs)
L = [sind(1:360).*prefs.orientWheelRadius + window.centerX;...
    cosd(1:360).*prefs.orientWheelRadius + window.centerY];
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
function maskwin = createMasktex(window,prefs)
% Set-up off screen window with mask
maskwin = Screen('OpenOffscreenwindow', window.onScreen, window.gray); 

% Set up alpha-blending for smooth (anti-aliased) lines in mask window
Screen('BlendFunction', maskwin, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% degslocs = prefs.degslocs; %deg of orientation in orientwheel for mask
degslocs = [15:15:360,360:-15:15]; %deg of orientation for drawing mask
% Draws each target stimulus to create mask
for ii = 1:length(degslocs)

    yadj = prefs.orientwheel(2,degslocs(ii)); %adjust lines y position 
    xadj = prefs.orientwheel(1,degslocs(ii)); %adjust lines x position 
    yadj2 = prefs.orientwheel2(2,degslocs(ii)); %for edge lines 
    xadj2 = prefs.orientwheel2(1,degslocs(ii)); %for edge lines

    % Set the coordinates (these are all relative to zero we will let
    % the drawing routine center to our monitor for us)
    % Set the new coordinates
    xCoords = [xadj -xadj2 xadj xadj2 xadj 0];
    yCoords = [yadj yadj2 yadj -yadj2 yadj 0];
    allCoords = [xCoords; yCoords];

    % Draw the orientation lines, set it to the center of our screen and
    % set good quality antialiasing
    Screen('DrawLines', maskwin, allCoords, prefs.lineWidthPix, window.black,...
        [window.centerX window.centerY], 2);
end

% Draw the circle to the screen
centeredRect = CenterRectOnPointd([0,0,14,14], window.centerX, window.centerY); 
Screen('FillOval', maskwin, window.black, centeredRect');

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
prefs.baseCircleresp = [0,0,10,10]; %Size of target's central circle in pixels
% Used to draw the targets and response wheel
prefs.orientwheel = [sind(1:360).*30; cosd(1:360).*30]; %all possible x & y positions for center line
prefs.orientwheel2 = [cosd(1:360).*2; sind(1:360).*2]; %all possible x & y positions edge lines
prefs.lineWidthPix = 2; %Set the line width for our response stimulus


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


% ---------------------
% Orientation Wheel ---
% ---------------------
% To estimate orientation based on location of mouse (not actually drawn)
prefs.orientWheelRadius = 35; %size of orientation wheel (not drawn)


% -------------------
% Trials & Blocks ---
% -------------------
% Randomize trial order of full factorial design order.
% prefs.fullFactorialDesign = fullfact([length(prefs.setSizes), ...
%     prefs.lags_per_block, ...
%     prefs.n_lags]);
% prefs.order = Shuffle(1:length(prefs.fullFactorialDesign));

% All possible orientations for targets
prefs.degslocs = [15:15:360]; 
% Number of blocks
prefs.nBlocks = 8;
% Trials per block (multiple of # of target orientations)
prefs.trials_per_block = length(prefs.degslocs)*2;
% Total number of trials (384)
prefs.nTrials = prefs.nBlocks*prefs.trials_per_block;
% prefs.nTrials = 12;

% Random order of orientations on each trial
prefs.selectorientTrial = randi(length(prefs.degslocs),[prefs.nTrials,1]);


end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
