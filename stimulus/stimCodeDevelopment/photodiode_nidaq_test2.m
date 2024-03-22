% photodiode_nidaq_test2.m
%
% Test script for running data acquisition in background

% clean up
clear all
close all
clc

% some parameters
numScans = 10000;

% initialization for NIDAQ
s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1','ai0','Voltage');
% s.DurationInSeconds = 1;
s.IsContinuous = true;
% s.NotifyWhenDataAvailableExceeds = 1000;

% note to self: check how accurate timing is using defaults for analog
%  input channel

% initialization for data collection
global pdData pdTime;
pdData = zeros(numScans,1);
pdTime = zeros(numScans,1);

% How to save data
% s.ScansAcquired
% global variable - preallocate given user inputed estimated time or just
% use much larger than necessary and clip at end
% use s.ScansAcquired to index into preallocated array and save data on
% each 'DataAvailable' event

keyCode = zeros(1,256);
escCode = zeros(1,256);
spaceCode = zeros(1,256);
escCode(27) = 1;
spaceCode(32) = 1;


% lh = addlistener(s,'DataAvailable', @stopOnEscPress);
lh = addlistener(s,'DataAvailable', @collectData);
% lh = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps,event.Data));
% start background acquisition
disp('Starting acquisition'); 
s.startBackground();
  
% can wrap psychtoolbox into a loop like this
% while ~isequal(escCode,keyCode) 
%     [~, ~, keyCode, ~] = KbCheck;
%     if isequal(escCode,keyCode)
%         disp('Stopping acquisition');
%         s.stop();
%     end
% end

% while ~isequal(spaceCode,keyCode)   
%     [~, ~, keyCode, ~] = KbCheck;
%     if isequal(spaceCode,keyCode) 
%         disp('Stopping acquisition');
%         s.stop();
%     end
% end

while s.IsRunning
    pause(0)
    [~, ~, keyCode, ~] = KbCheck;
    if (isequal(spaceCode,keyCode))
        disp('Stopping acquisition');
        s.stop();
    end
end

pdData = pdData(1:s.ScansAcquired);
pdTime = pdTime(1:s.ScansAcquired);
plot(pdTime,pdData);

delete(lh);

% clear global variables at end
% clear global pdData
% clear global pdTime

