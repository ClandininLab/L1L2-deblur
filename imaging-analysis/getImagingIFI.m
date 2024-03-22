% getImagingIFI.m
% 
% Calculates frame rate and times at which each frame occuredbased on 
%  trigger pulses received by the NIDAQ card from the imaging computer.
% Has the option to plot imaging frame time parameters:
%   -histogram of the distribution of interframe intervals
%   -all the IFIs
%   -counter for which frame was happening at which time
%   -time at which each frame occured
%
% INPUT: 
%   obj - Stimulus object
%   showPlots - boolean to control whether plots are displayed
%
% OUTPUT: 
%   ifi - imaging interframe interval
%   frameTimes - time at which each imaging frame started
%

function [ifi, frameTimes] = getImagingIFI(obj, showPlots)

    % frame count at the NIDAQ scan rate
    imft = obj.Out.imFrameTime; 
    % scans where a new frame occured
    indNewFrame = find(diff(imft) > 0) + 1; 
    % time at which new frame occured in seconds (onset)
    frameTimes = obj.Out.pdTime(indNewFrame); 
    allIFIs = diff(frameTimes); 
    ifi = mean(allIFIs);

    if showPlots
        % Plot interframe interval
        figure; 
        hist(allIFIs); 
        title('ifi distribution'); 
        xlabel('seconds'); 
        
        figure; 
        plot(allIFIs); 
        title('ifi'); 
        xlabel('frame number'); 
        ylabel('ifi (seconds)');
        
        figure; 
        plot(obj.Out.pdTime, imft, 'b.'); 
        title('frame count recorded by NIDAQ scan');
        ylabel('frame number'); 
        xlabel('time (sec)');

        figure;
        plot(frameTimes, 'bo'); 
        title('frame onset time (seconds)'); 
        ylabel('frame onset time (seconds)'); 
        xlabel('frame number');
    end 

end 