% playStimMain.m 
%
% This is the main script for visual stimulus presentation. It initializes
% the NIDAQ device and Psychtoolbox. The user is prompted for a 
% stimulus .txt file and a Stimulus object is created based on
% the name of the stimulus in the file. The stimulus is displayed using the
% methods associated with the Stimulus object. Outputs stim.mat, which 
% is a struct containing the stimulus parameters, the raw stimulus values 
% presented, and the photodiode output.
%
% INPUT: 
%  screenNumber - screen number corresponding to lightcrafter - check what
%   is displayed in PTB-INFO when this is run for the first time and what
%   Windows names the screen
% 
% last updated 8/16/17 HHY
% Previous version of this code: playstim_v2.m

function playStimMain(screenNumber)

% Clean up
% clear all
clear functions;
close all
clc
Screen('CloseAll');
clearvars -except screenNumber

% set default figure position to make them display on this monitor
set(0,'defaultfigureposition',[-1500,300,600,500]);

% Make sure the path is correct for your computer to access all the
% functions in this code
homeDir = 'C:\Users\ClandininLab\Desktop'; % default directory for saving data
stimTxtDir = 'C:\Users\ClandininLab\Documents\2p-stim-code'; % location of stimulus txt files
stimCodeDir = 'C:\Users\ClandininLab\Documents\2p-stim-code'; % main scripts and everything else
addpath(stimCodeDir);

% Initialize NIDAQ 
[nidaq, lh] = initNidaq();

% Request stimulus .txt file name from user
% sFileName = input('What is the name of your stimulus .txt file (include the .txt part)? ', 's');
fprintf('\nSelect your stimulus .txt file. \n'); 
cd(stimTxtDir);
sFileName = uigetfile('*.txt', 'Select your stimulus .txt file.');

% Read the .txt file to find the appropriate stimulus class and to create 
% a stimulus object 
stim = getStim(sFileName);

fprintf(sFileName); % MMP added for sanity check
fprintf('\n');

% Ask the user for the intended total stimulus duration 
d = input('Enter duration of stimulus presentation (sec): ', 's');
stimDuration = str2double(d);

% save the NIDAQ scan rate
stim.setNIDAQScanRate(nidaq.Rate);
% NIDAQ scans/sec multiplied by seconds of total stimulus duration + some padding
numScans = nidaq.Rate*stimDuration + 5000; 

% initialization for data collection
global pdData pdTime imFrameTime; 
pdData = zeros(numScans,1);
pdTime = zeros(numScans,1);
imFrameTime = zeros(numScans,1);

% initialization for Psychtoolbox
[window, winRect, screenX, screenY, ifi] = initPsychTbx(screenNumber);
stim.IFI = ifi; % save the interframe interval

Priority(MaxPriority(window));

% initialize stimulus if stimulus has initialization code
if ismethod(stim,'initStim')
    disp('Initializing Stimulus');
    stim.initStim(window,ifi);
end

% External trigger to start stimulus 
% Note: moved to after PTB initialization, rewrote code to handle
disp('To start the stimulus, press ENTER');
enterKey = KbName('return');
while 1
    [~,~,keyCode] = KbCheck;
    if keyCode(enterKey)
        pause(0.1);
        break;
    end
end
display('To stop the stimulus manually, press any other key.');

% start background acquisition for NIDAQ
disp('Starting acquisition'); 
nidaq.startBackground();

% display stimulus 
stim.displayStim(window, ifi, stimDuration);

% stop NIDAQ background acquisition 
nidaq.stop();

% close the screen in Psychtoolbox
Screen('CloseAll');
Priority(0);

% keep the photodiode output up to the point stimulus presentation stops
pdData = pdData(1:nidaq.ScansAcquired);
pdTime = pdTime(1:nidaq.ScansAcquired);
imFrameTime = imFrameTime(1:nidaq.ScansAcquired);
plot(pdTime, pdData);

% Ask the user whether or not to save the data
saveCommand = input('Would you like to save this data? [n to cancel, any other key to save] \n', 's');
if ~strcmp(saveCommand, 'n')
    stim.saveData(pdData, pdTime, imFrameTime); % saves photodiode output and timing info from imaging computer 
    stim.exportData(homeDir, stimCodeDir); % exports stimulus parameters and output data to a folder
end 

delete(lh);

% clear global variables at end
clear global pdData
clear global pdTime
clear global imFrameTime

end 


