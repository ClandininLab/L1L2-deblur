% photodiode NiDAQ test 1
%
% Very simple, just get signal from photodiode
%
% 7/18/16

clear all
s = daq.createSession('ni');
% daq.getDevices

s.DurationInSeconds = 5;

addAnalogInputChannel(s,'Dev1','ai0','Voltage');
s.NotifyWhenDataAvailableExceeds = 1000;
lh = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps,event.Data));

% data = startForeground(s);
startBackground(s);
plot(data)