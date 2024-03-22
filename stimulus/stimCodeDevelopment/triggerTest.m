
% Code for setting trigger inputs and outputs 

clear all;
% release(d);

d = daq.createSession('ni');
ch = addCounterInputChannel(d, 'Dev1', 'ctr1', 'EdgeCount'); % read triggers from imaging computer 
aCh = addAnalogInputChannel(d,'Dev1','ai0','Voltage'); % clock
oCh = addCounterOutputChannel(d,'Dev1','ctr0','PulseGeneration'); % out trigger to start imaging
oCh.Frequency = 1;
oCh.InitialDelay = 3;
ch.ActiveEdge = 'Rising'; 
d.DurationInSeconds = 10;
d.Rate = 5000;
[data, timeStamps, triggerTime] = startForeground(d);
plot(timeStamps,data)

removeChannel(d,1:length(d.Channels));
release(d);


% add code to save all this information from the imaging computer - data,
% timeStamps, triggerTime - save one trigger Time data point to know what
% date/time the experiment was done. Save all in stim.Out

% Ask user where to save all the data into 
% Ask user for name of the stim.Out to save in a folder or something. So
% that the user doesn't have to physically write down the stimulus output
% file names 



% Leica can save images as .lif and not .tiff 
% .lif contains a bunch of LDM foldders 

% Have user type in the name of the LDM folder (e.g. 'LDM_005'). that name
% should be associated with the stimulus output file somehow - like in its
% name or separate folder that the anlaysis reads 

% associate specific stimulus ouput file with LDM Folder 


% To Do: save timestamps at which we received a pulse, not every 