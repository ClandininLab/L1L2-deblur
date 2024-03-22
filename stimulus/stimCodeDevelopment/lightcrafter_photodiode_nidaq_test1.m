% lightcrafter_photodiode_nidaq_test1.m
%
% Test script for running data acquisition in background while
% simultaneously presenting visual stimulus using psychtoolbox and
% lightcrafter
%
% Note to self (160724) - test 300 Hz frame presentation

% clean up
clear all
close all
clc
Screen('CloseAll');
clearvars;

% some parameters
numScans = 1000000;

% initialization for NIDAQ
s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1','ai0','Voltage');
s.Rate = 5000; % scans per second
s.IsContinuous = true;
s.NotifyWhenDataAvailableExceeds = 1000;

% listener to collect photodiode data
lh = addlistener(s,'DataAvailable', @collectData);

% initialization for data collection
global pdData pdTime;
pdData = zeros(numScans,1);
pdTime = zeros(numScans,1);

% some more initialization for stopping everything
keyCode = zeros(1,256);
escCode = zeros(1,256);
escCode(27) = 1;

% initialization for Psychtoolbox
PsychDefaultSetup(2);
screens = Screen('Screens');
% limits sync tests to critical ones
Screen('Preference', 'SkipSyncTests', 0);
% makes initial screen black
Screen('Preference', 'VisualDebugLevel',1);

% pre visual stimulus setup
screenNumber = max(screens);
% screenNumber = 1;
black = BlackIndex(screenNumber) * (64/256);
white = WhiteIndex(screenNumber) * (64/256);
white = WhiteIndex(screenNumber);
% Open an on screen window and color it black
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% Get the size of the on screen window in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
% get refresh rate of lightcrafter
ifi = Screen('GetFlipInterval', window);

% number of frames to wait (how often flip happens)
waitframes = 1;

% how often stimulus flashes back and forth
flashFrames = 2;

% set max priority
Priority(MaxPriority(window));


% start background acquisition
disp('Starting acquisition'); 
s.startBackground();

% frameCounter = 0;
% while s.IsRunning
%     if (mod(frameCounter,flashFrames*2) < flashFrames)
%         Screen('FillRect', window, white);
%         vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
%     else
%         Screen('FillRect', window, black);
%         vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
%     end
%     
%     % checks this every frame
%     pause(0)
%     [~, ~, keyCode, ~] = KbCheck;
%     if (isequal(escCode,keyCode))
%         disp('Stopping acquisition');
%         s.stop();
%         break; % may not be necessary, but force the break
%     end
%     frameCounter = frameCounter + 1;
% end

% prevFrame = 0;
% waitframes = 50;
% while s.IsRunning
%     if (prevFrame == 0)
%         prevFrame = 1;
%         Screen('FillRect', window, white);
%         vbl = Screen('Flip', window, vbl + (waitframes-1)*ifi + 0.001);
%     else
%         prevFrame = 0;
%         Screen('FillRect', window, black);
%         vbl = Screen('Flip', window, vbl + (waitframes-1)*ifi + 0.001);
%     end
%     
%     % checks this every frame
%     pause(0)
%     [~, ~, keyCode, ~] = KbCheck;
%     if (isequal(escCode,keyCode))
%         disp('Stopping acquisition');
%         s.stop();
%         break; % may not be necessary, but force the break
%     end
% end

numFlips = 400;
waitframes = 2;
% vbl = Screen('Flip', window); % gives a time stamp
startTime = GetSecs();
for i = 1:numFlips
    if (mod(i,2))
        Screen('FillRect', window, white);
%         vbl = Screen('Flip', window, vbl + (waitframes-1)*ifi + 0.001);
        Screen('Flip', window, startTime + i*ifi + (waitframes-1)*ifi + 0.001);
%          Screen('Flip',window);
%         Screen('Flip',window,0,2);
    else
        Screen('FillRect', window, black);
%         vbl = Screen('Flip', window, vbl + (waitframes-1)*ifi + 0.001);
        Screen('Flip', window, startTime + i*ifi + (waitframes-1)*ifi + 0.001);
%         Screen('Flip',window);
%         Screen('Flip',window,0,2);
    end
    pause(0) % this needs to be in here for correct data acquisition; unclear why  
%     % checks this every frame
%     pause(0)
%     [~, ~, keyCode, ~] = KbCheck;
%     if (isequal(escCode,keyCode))
%         disp('Stopping acquisition');
%         s.stop();
%         break; % may not be necessary, but force the break
%     end
end
s.stop();

Screen('CloseAll');
Priority(0);

pdData = pdData(1:s.ScansAcquired);
pdTime = pdTime(1:s.ScansAcquired);
plot(pdTime,pdData);

delete(lh);

% clear global variables at end
% clear global pdData
% clear global pdTime

