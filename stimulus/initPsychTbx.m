% initPsychTbx.m
% 
% Uses default settings to set up Psychtoolbox. Detects screens linked to 
% the stimulus computer. Opens a screen window and sets it to the lowest
% luminance possible. Measures several parameters which are needed for
% stimulus presentation in playStimMain.
%
% INPUT: 
%  screenNumber - screen number corresponding to lightcrafter - check what
%   is displayed in PTB-INFO when this is run for the first time and what
%   Windows names the screen
%
% OUTPUTs:
%   window - pointer to the onscreen window
%   screenXpixels - total number of pixels along the X axis of the screen
%   screenYpixels - total number of pixels along the Y axis of the screen
%   ifi - inter-frame interval, aka screen refresh rate

function [window, windowRect, screenXpixels, screenYpixels, ...
    ifi] =  initPsychTbx(screenNumber)

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
% For help see: Screen Screens?
% screens = Screen('Screens');
% limits sync tests to critical ones
% Screen('Preference', 'SkipSyncTests', 0);
% suppress printout of all warnings
Screen('Preference', 'SuppressAllWarnings',1);
% Screen('Preference', 'VBLTimestampingMode', -1);
% makes initial screen black
Screen('Preference', 'VisualDebugLevel',1);

% If two screens are attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% that.
% screenNumber = max(screens);
% screenNumber = 1; % lightcrafter is screen 1 - 161223
% screenNumber = 2; % lightcrafter is now screen 2 - 170727 ?HHY

% Define black (white will be 1 and black 0). This is because
% luminace values are (in general) defined between 0 and 1.
% For help see: help BlackIndex
black = BlackIndex(screenNumber);
% white = WhiteIndex(screenNumber);

% Open an on screen window and color it black
% For help see: Screen OpenWindow?
% or go to http://docs.psychtoolbox.org/OpenWindow
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Measure the vertical refresh rate of the monitor
ifi = Screen('GetFlipInterval', window); % 80Hz

end 