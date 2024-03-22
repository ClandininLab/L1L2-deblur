% screensizeTest_RGB.m
%
% INPUT: 
%  screenNumber - screen number corresponding to lightcrafter - check what
%   is displayed in PTB-INFO when this is run for the first time and what
%   Windows names the screen
%
function screensizeTest_RGB(screenNumber)
    % Clear the workspace and the screen
    sca;
    close all;
    % clear all;
    clear functions;
    clearvars -except screenNumber
    clc;

    % Here we call some default settings for setting up Psychtoolbox
    PsychDefaultSetup(2);

    % Get the screen numbers. This gives us a number for each of the screens
    % attached to our computer.
    % For help see: Screen Screens?
%     screens = Screen('Screens');
    Screen('Preference', 'VisualDebugLevel',1);


    % Draw we select the maximum of these numbers. So in a situation where we
    % have two screens attached to our monitor we will draw to the external
    % screen. When only one screen is attached to the monitor we will draw to
    % this.
    % For help see: help max
    % screenNumber = max(screens);
    % screenNumber = 1;
    % screenNumber = 2;

    % Open an on screen window and color it black
    % For help see: Screen OpenWindow?
    black = [0 0 0];
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
    % Screen('ColorRange',window,[],-1,[]);

    % screen
    centeredRect = [130 0 700 1140];
    pdRect = [30 0 60 60];
    val = 1; % present full intensity in screensizeTest
    colScreen = convertInt2RGB(repmat(val,1,3));

    Screen('FillRect', window, colScreen, centeredRect);
    Screen('FillRect', window, colScreen, pdRect);

    Screen('Flip', window);

    % Now we have drawn to the screen we wait for a keyboard button press (any
    % key) to terminate the demo.
    % For help see: help KbStrokeWait
    KbStrokeWait;

    % Clear the screen. "sca" is short hand for "Screen CloseAll". This clears
    % all features related to PTB. Note: we leave the variables in the
    % workspace so you can have a look at them if you want.
    % For help see: help sca
    sca;

end