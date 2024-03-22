% processPD.m
%
% Function that takes Stimulus object saved by the visual stimulus
%  presentation code during imaging and uses the photodiode timing and the
%  info saved about the stimulus to reconstruct what stimulus was presented
%  at what time.
%
% INPUTS:
%   stim - Stimulus object
%   imFrameStartTimes - times at which imaging frames were captured
%
% OUTPUTS:
%   stim - Stimulus object (can be modified wtih call to reconstructStim)
%   stimvalIF - stimulus value at each imaging frame, for prelim plotting
%   stimTimesChosen - stimulus times for each of those stimulus values,
%       again for prelim plotting
%   rcStim - reconstructed stimulus value on each epoch
%   rcStimInd - reconstructed order of epochs
%   stimEpochStartTimes - times of dark to light as well as light to dark
%       transitions (sorted concatenation of lightStartTimes and
%       darkStartTimes)
%   lightStartTimes - times of dark to light transitions
%   darkStartTimes - times of light to dark transitions
%   thresh - threshold voltage to separate light from dark in photodiode
%       output
%
% Last updated: 3/28/17 HHY (changed stimulus reconstruction)
%

function [stim, stimvalIF, stimTimesChosen, rcStim, rcStimInd, ...
        stimEpochStartTimes, lightStartTimes, darkStartTimes, thresh] ...
        = processPD(stim, imFrameStartTimes)
% check a few stim object properties
    if (isprop(stim.obj, 'NIDAQScanRate') && ~isempty(stim.obj.NIDAQScanRate))
        nidaqScanRate = stim.obj.NIDAQScanRate;
    else 
        nidaqScanRate = 5000;
    end 

    if (isprop(stim.obj, 'IFI') && ~isempty(stim.obj.IFI))
        stimIFI = stim.obj.IFI;
    else 
        stimIFI = 0.01;
    end 

    pdData = stim.obj.Out.pdData;
    pdTime = stim.obj.Out.pdTime;

%% Ask user to select threshold betwee light and dark by clicking on the
% photodiode plot. Set a default one. Give the option of accepting the
% default threshold upon pressing ENTER
    display('Zoom in on a dark-light transition in the photodiode signal.');
    display('Hit SPACE, then select a voltage threshold to separate dark from light epochs.');

    fig = figure; 
    plot(pdTime, pdData);
    ylabel('Voltage');
    xlabel('Time (seconds)');
    title('Photodiode signal');
    keypress = 0;

    while ~keypress
        keypress = waitforbuttonpress;
    end 
    [~, thresh] = ginput(1);
    close(fig);

    fprintf('Light-dark threshold is set at %f Volts. \n', thresh);

%% Get the timestamps at which stimulus changes 
% This step is common to all stimulus classes 
    [stimEpochStartTimes, lightStartTimes, darkStartTimes, order] = ...
        getStimEpochStartTimes(pdTime, pdData, imFrameStartTimes, ...
        thresh, nidaqScanRate, 1);

%% Stimulus reconstruction
% changed 3/7/17 to use method associated with each Stimulus subclass - HHY

    [rcStim, rcStimInd] = stim.obj.reconstructStim(stimEpochStartTimes, ...
        lightStartTimes, darkStartTimes, order);

    % At this point, epochStartTimes, rcStim, and imFrameStartTimes
    %  allow us to calculate fstimval.

%% Find the stimulus value at each imaging frame 
% (used only for prelim. dF/F visualization)

    % initialize
    stimvalIF = zeros(length(imFrameStartTimes), 1);
    % set the stimulus during the first imaging from to 0 - we're going
    % to throw out this frame anyway
    stimvalIF(1) = 0;
    
    % Save the stimulus timestamps corresponding to stimvalIF
    stimTimesChosen = zeros(length(imFrameStartTimes), 1);

    for i = 2:length(imFrameStartTimes)
%         i % debug
        prevStarts = find(stimEpochStartTimes < imFrameStartTimes(i));
        
        % if no stimulus onset before this imaging frame is detected, set
        % the stimulus value to 0
        if isempty(prevStarts)
            stimvalIF(i) = 0;
            continue
        end 
        
        candVals = rcStim(prevStarts);
        [m, subInd] = max(stimEpochStartTimes(prevStarts));
        stimvalIF(i) = candVals(subInd); 

        % debug
        stimTimesChosen(i) = m;
    end 
    
    % stimulus time stamps chose
    figure; plot(stimTimesChosen); 
    title('stimulus time stamps chosen');
    ylabel('stimulus time stamps (sec)'); xlabel('imaging frame');

    % plot stimvalIF
    figure; plot(stimvalIF); 
    title('Stimulus values at each imaging frame');
    ylabel('stimulus value'); xlabel('imaging frame');

    clear pdTime pdData
end 
