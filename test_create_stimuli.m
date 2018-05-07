

% Clear the workspace and the screen
sca;
close all;
clearvars;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Draw to the external screen if avaliable
screenNumber = max(screens);

% Define black and white
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = white/2; %create grey for background

% Open an on screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get center of screen
centerX = screenXpixels * 0.5; % center of screen in X direction
centerY = screenYpixels * 0.5; % center of screen in Y direction

% Query the frame duration
refresh = Screen('GetFlipInterval', window);
slack = refresh/2; % Divide by 2 to get slack

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Make a base Rect in pixels
baseRect = [0,10,5,40;0,0,12,12];

% Center the oval on the centre of the screen
% centeredRect = CenterRectOnPointd(baseRect(1,:), xCenter, yCenter-12);

% Set the color of the rect to red
% rectColor = [1 0 0];
rectColor = black; %set color to black

% Draw the rect to the screen
% Screen('FillOval', window, rectColor, centeredRect');

% Draw line
% Screen('DrawLine', window, color, fromH, fromV, toH, toV, linewidth);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter-1, yCenter, 2);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter+1, yCenter, 2);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter-2, yCenter, 2);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter+2, yCenter, 2);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter-3, yCenter, 2);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter+3, yCenter, 2);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter-4, yCenter, 2);
Screen('DrawLine', window, rectColor, xCenter, yCenter+30, xCenter+4, yCenter, 2);


% Center the oval on the centre of the screen
centeredRect = CenterRectOnPointd(baseRect(2,:), xCenter, yCenter);

% Draw the rect to the screen
Screen('FillOval', window, rectColor, centeredRect');

% Flip to the screen
Screen('Flip', window);

% =========================================================================
% Wait for click.
SetMouse(xCenter, yCenter); %place the mouse in the center of screen
ShowCursor('Arrow');

% If mouse button is already down, wait for release.
GetMouse(window); %location of mouse relative to screen dimensions
buttons = 0;
while any(buttons)
    [x, y, buttons] = GetMouse(window);
end

% everMovedFromCenter = false;
% while ~any(buttons)
% 
% %     drawColorWheel(window, prefs, prefs.colorwheel_trial{trialIndex});
% 
%     [x,y,buttons] = GetMouse(window);
%     [minDistance, minDistanceIndex] = min(sqrt((colorWheelLocations(1, :) - x).^2 + (colorWheelLocations(2, :) - y).^2));
% 
%     if(minDistance < 500)
%         everMovedFromCenter = true;
%     end
% 
%     if(everMovedFromCenter)
% %                 colorsOfTest(itemToTest(trialIndex), :) = prefs.colorwheel(minDistanceIndex,:);
%         colorsOfTest(itemToTest(trialIndex), :) = prefs.colorwheel_trial{trialIndex}(minDistanceIndex,:);
%     else
%         colorsOfTest(itemToTest(trialIndex), :) = [145 145 145];
%     end
% 
%     Screen('FillOval', window.onScreen, colorsOfTest', rects{nItems});
% 
%     drawColorWheel(window, prefs, prefs.colorwheel_trial{trialIndex});
%     Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
%     Screen('Flip', window);
% 
% end
%         

% Wait for a key press
KbStrokeWait;

% Clear the screen
sca;




















