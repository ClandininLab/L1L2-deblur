% lightcrafter_photodiode_nidaq_test3.m
%
% Use to measure intensity at each 8-bit value
% Use OpenGL to present values

% clean up
clear all
close all
clc
Screen('CloseAll');
clearvars;

InitializeMatlabOpenGL();

% some parameters
numScans = 600*5000;

% initialization for NIDAQ
s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1','ai0','Voltage');
s.Rate = 5000; % scans per second
s.IsContinuous = true;
s.NotifyWhenDataAvailableExceeds = 1000;

% listener to collect photodiode data
lh = addlistener(s,'DataAvailable', @collectData2);

% initialization for data collection
global pdData pdTime;
pdData = zeros(numScans,1);
pdTime = zeros(numScans,1);

% some more initialization for stopping everything
keyCode = zeros(1,256);
escCode = zeros(1,256);
escCode(27) = 1;

% initialization for Psychtoolbox
PsychDefaultSetup(1);
screens = Screen('Screens');
% limits sync tests to critical ones
Screen('Preference', 'SkipSyncTests', 0);
% makes initial screen black
Screen('Preference', 'VisualDebugLevel',1);

% pre visual stimulus setup
screenNumber = max(screens);
% screenNumber = 1;
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
% Open an on screen window and color it black
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% Get the size of the on screen window in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
% get refresh rate of lightcrafter
ifi = Screen('GetFlipInterval', window);

% number of frames to wait (how often flip happens)
waitframes = 80;

% set max priority
Priority(MaxPriority(window));


% all values to test
numLevels = 256;
vals = 0:(numLevels-1);
vals = [vals 0];
% vals = [4 5 1 127 128 255 0];
% numLevels = length(vals);

% start background acquisition
disp('Starting acquisition'); 
s.startBackground();

epochDur = 1;
startTime = Screen('Flip', window);

frame = 0;
currEpoch = 0;
% keep playing the stimulus until the user presses a key on the keyboard
while (currEpoch < (length(vals)*2))
    currEpoch = floor((frame*ifi)/epochDur);
    % light epoch
    if (mod(currEpoch, 2))
        Screen('BeginOpenGL',window);
        col = vals((currEpoch-1)/2+1);
        col = col/255;
        glColor3f(col,col,col);
        glRectd(0,0,912,1140);
        Screen('EndOpenGL',window);
        Screen('Flip', window, startTime + frame*ifi); % flip window at the frame rate
    % dark epoch
    else
        Screen('BeginOpenGL',window);
        glColor3f(0,0,0);
        glRectd(0,0,912,1140);
        Screen('EndOpenGL',window);
        Screen('Flip', window, startTime + frame*ifi);  
    end 
    pause(0) % super important - without this NIDAQ won't scan
    frame = frame + 1;
end 

s.stop();

Screen('CloseAll');
Priority(0);

pdData = pdData(1:s.ScansAcquired);
pdTime = pdTime(1:s.ScansAcquired);
plot(pdTime,pdData);

delete(lh);

removeChannel(s,1);
release(s);


% compute mean intensities
startInd = 5500;
numToAvg = 4000;
sep = 10000;

meanInten = zeros(1,numLevels);

for i=1:numLevels
    stInd = startInd + (i-1)*sep;
    endInd = stInd + numToAvg-1;
    meanInten(i) = mean(pdData(stInd:endInd));
end

