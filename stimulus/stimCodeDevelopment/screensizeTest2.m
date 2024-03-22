% Clear the workspace and the screen
sca;
close all;
clear all;
clearvars;
clc;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
% For help see: Screen Screens?
screens = Screen('Screens');
Screen('Preference', 'VisualDebugLevel',1);

% Draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% this.
% For help see: help max
screenNumber = max(screens);

% Define black (white will be 1 and black 0). This is because
% luminace values are (in general) defined between 0 and 1.
% For help see: help BlackIndex
black = BlackIndex(screenNumber);

% Open an on screen window and color it black
% For help see: Screen OpenWindow?
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% screen
centeredRect = [160 100 680 1140];

% some parameters
epochDur = 2;
ifi = 1/80;
darkContrast = 0;
lightContrast = 1;

startTime = Screen('Flip', window);

frame = 0;
% keep playing the stimulus until the user presses a key on the keyboard
while ~KbCheck 
    currEpoch = floor((frame*ifi)/epochDur);
    % light epoch
    if (mod(currEpoch, 2))
        Screen('FillRect', window, lightContrast, centeredRect);
        Screen('Flip', window, startTime + frame*ifi); % flip window at the frame rate
    % dark epoch
    else
        Screen('FillRect', window, darkContrast, centeredRect);
        Screen('Flip', window, startTime + frame*ifi);  
    end 
    pause(0) % super important - without this NIDAQ won't scan
    frame = frame + 1;
end 

% Now we have drawn to the screen we wait for a keyboard button press (any
% key) to terminate the demo.
% For help see: help KbStrokeWait
KbStrokeWait;

% Clear the screen. "sca" is short hand for "Screen CloseAll". This clears
% all features related to PTB. Note: we leave the variables in the
% workspace so you can have a look at them if you want.
% For help see: help sca
sca;
