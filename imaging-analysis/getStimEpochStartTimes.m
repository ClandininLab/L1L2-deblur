% getStimEpochStartTimes.m
% 
% Timestamps (seconds) at which the stimulus started a new epoch 
% according to the photodiode on the NIDAQ timeline. An epoch transition
% is defined as an update in stimulus contrast value. 
%
% Clarification of confusing terminology: 
%      "start of light epoch" = dark-to-light transition (the timestamp
%      actually falls right between the dark and the light epochs, the point 
%      closest to threshold) 
%      "start of dark epoch" = light-to-dark transition 
% 
% Algorithm: 
% 1. Label each stimulus frame as either light (greater than threshold)
% or dark (less than threshold).
% 2. Take the difference between all adjacent light timestamps (should be a
% vector one element shorter than the number of light timestamps).
% 3. What distinguishes between the start of a new light epoch and the
% start of a new light frame within the same epoch is their distance from the
% most recent light timestamp. The distance from the start of a light epoch
% and the most recent light time is greater than the distance between two
% light frames. The inter-lightframe distance is currently hard-coded as 
% 0.0006s (see interFrameDist). 
% The code labels the indices of "true" light epoch start times. 
% Note that the first time isD2L is calculated, it labels the start of the
% dark epochs due to how the diff() function indexes. So shifting the
% indices to the right gives the light epoch start times. 
% 4. Get the timestamps of the true light start times.
% 5. To find the start of the dark epochs, we take the indices we used from
% the diff() function, and plug the indices into the dark timestamps.
% 
% The current algorithm fails record the very last 
% light-to-dark transition. However, this should not be a problem for 
% analysis since imaging always stops before the stimulus stops. If we are
% short on data, it's always on the imaging end. 
%
% ASSUMPTIONS
% 1. A frame is 3 patterns long.
% 2. The distance between frames is around 0.0006s. This value changes with
% the pattern duration set on the LightCrafter software!
% 
% PLOTS:
%   -times when photodiode signal is called as light or dark
%   -just the photodiode times labeled light
%   -plot light and dark start times
%   -plot times when stimulus changes 
%   -overlay of stimulus transitions and imaging frame times on top of
%   photodiode voltage signal over time, first 4 seconds
% 
%
% INPUTS:
%   pdTime - pdTime, from Stimulus object; times at which photodiode values
%       were captured
%   pdData - pdData, from Stimulus object; photodiode voltage values
%   imFrameStartTimes - times when imaging frames were captured
%   thresh - voltage threshold separating dark from light
%   NIDAQScanRate - as named, how often NIDAQ card captured data
%   makePlot - boolean for whether or not to display plots
%
% OUTPUTS:
%   epochStartTimes - times of dark to light as well as light to dark
%       transitions (sorted concatenation of lightStartTimes and
%       darkStartTimes)
%   lightStartTimes - times of dark to light transitions
%   darkStartTimes - times of light to dark transitions
%   order - indicies from before sorting concatenation of lightStartTimes 
%       and darkStartTimes, how are they ordered in epochStartTimes
%
%   220513 MMP edit b/c boolean(1) doesn't work in my MATLAB? lol

function [epochStartTimes, lightStartTimes, darkStartTimes, order] = ...
    getStimEpochStartTimes(pdTime, pdData, imFrameStartTimes, thresh, ...
    NIDAQScanRate, makePlot)

    nidaqScanIFI = 1/NIDAQScanRate;
    
    % indices in pData that are greater than threshold (aka light epochs). 
    % lights are labeled as 1s, darks labeled as 0s
    ldIndex = pdData > thresh; 
    
    buffer = 0.004; % higher the threshold means longer buffer (seconds)
    interFrameDist = 0.0006 + buffer; % distance between frames + a buffer (seconds)
%     patternLength = IFI/3; % 3 is the number of patterns in a frame
  
    % pdTimes labeled as light epoch
    timesLabeledLight = pdTime(ldIndex); 
%     timesLabeledDark = pdTime(~ldIndex);
    
    % Identify true light-dark-light transitions, defined as transitions 
    % longer than the interframe distance. The result of the following
    % operation is the indices for the start of dark epochs. 
    isD2LTrans = diff(timesLabeledLight) >= interFrameDist; 
    % a true dark period is longer than twice the length of 2 patterns 
%     isL2DTrans = diff(timesLabeledDark) >= patternLength*2; 

    % To find start of light epochs, shift one index one position
    % to right due to the diff() function and to account for the first
    % light epoch (thus the 'boolean(1)'). 
%    isLightStart = [boolean(1); isD2LTrans]; 
    isLightStart = [true; isD2LTrans]; 
    
    
    % pd times at start of light epoch
    lightStartTimes = timesLabeledLight(isLightStart); 
    darkStartTimes = timesLabeledLight(isD2LTrans) + nidaqScanIFI;
    
    % combine light and dark epoch start times to get times when the
    % stimulus changed. Keep the sorted order of the old indices. 
    [epochStartTimes, order] = sort([lightStartTimes; darkStartTimes]);
    
    % make plots if user said so
    if (makePlot)
     % label photodiode signals as light or dark
        figure; plot(ldIndex); 
        title('Photodiode signal labeled as light or dark');
        
     % plot just the photodiode times labeled light
        figure; plot(timesLabeledLight, 'b.'); 
        title('Timestamps of light epoch onsets');
        
     % plot light and dark start times
        figure; plot(lightStartTimes, 'r.'); hold on;
        plot(darkStartTimes, 'k.');
        legend('light phase start times', 'dark phase start times');
        title('Light and dark phase time stamps');
        
     % plot times when stimulus changes 
        figure; plot(epochStartTimes, 'bo');
        title('Onsets of stimulus change');

     % visualize epoch start times on pd plot
        xRange = 4/nidaqScanIFI; % 4 seconds
        plot(pdTime(1:xRange), pdData(1:xRange)); hold on;
        ylabel('voltage'); xlabel('time (sec)');
        ylim([-0.01 0.05]);

        for i = 1:length(find(epochStartTimes < pdTime(xRange)))
            line([epochStartTimes(i) epochStartTimes(i)], [0 thresh + .02], ...
                'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
            hold on;
        end
        
%         title('Stimulus transition times');

        % vizualize imaging frame times too
        for j = 1:length(find(imFrameStartTimes < pdTime(xRange)))
            line([imFrameStartTimes(j) imFrameStartTimes(j)], [0 thresh + .02], ...
                'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k'); 
            hold on;
        end 
        
        title('Transition times for stimulus epochs and imaging frames');
        
     % -------- % Plot the last transition times for stimulus and imaging % -------- %
     
%         % plot times when stimulus changes 
%         figure; plot(epochStartTimes, 'bo');
%         title('Onsets of stimulus change for last second of stimulus');
% 
%         % visualize epoch start times on pd plot
%         xRange = 100/nidaqScanIFI; % last 2 min
%         plot(pdTime(length(pdTime)-xRange:end), pdData(length(pdData)-xRange:end)); hold on;
%         ylabel('voltage'); xlabel('time (sec)');
%         
%         % last stimulus transition recorded 
%         line([epochStartTimes(end) epochStartTimes(end)],...
%             [0 thresh + .02], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
%         ylim([-0.01 0.05]);
%         hold on;
%         
%         % last imaging frame transition recorded 
%         line([imFrameStartTimes(end) imFrameStartTimes(end)], ...
%                 [0 thresh + .02], 'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'k'); 

    
    end 

    
end 

