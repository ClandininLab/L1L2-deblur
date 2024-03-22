% initNidaq.m
% 
% Creates a new NIDAQ session. Adds input channels to count trigger pulses
% from the imaging computer and output channels to signal the start of
% imaging. Number of scans/second is current hard-coded as 5000.
% Ideally, we'd run this once before start of the day's imaging session,
% but that seems to create problems with saving the data. So we are running
% it for every experiment now (every time playStimMain is called)

function [d, lh] = initNidaq()

d = daq.createSession('ni');

% analog channel for the clock
aCh = addAnalogInputChannel(d,'Dev1','ai0','Voltage'); 
% this channel keeps a counter to count trigger pulses from the 
% imaging computer which indicates frame capture 
iCh = addCounterInputChannel(d, 'Dev1', 'ctr1', 'EdgeCount');
iCh.ActiveEdge = 'Rising';
% this channel sends triggers out to imaging computer to signal start
% of imaging 
% Sends voltage pulses at the specified frequency 
oCh = addCounterOutputChannel(d,'Dev1','ctr0','PulseGeneration'); 
% oCh.Frequency = 1; % why does this matter? 
% Time between start of background acquisition and sending of a pulse 
% to imaging computer 
oCh.InitialDelay = 0.050; % 50ms delay 

% change voltage range measured by analog input to -1 to +1 V
aCh.Range = [-1, 1];

d.Rate = 5000; % number of NIDAQ scans per second
d.IsContinuous = true;
d.NotifyWhenDataAvailableExceeds = 1000;

% listener to collect photodiode data
lh = addlistener(d,'DataAvailable', @collectData);

end