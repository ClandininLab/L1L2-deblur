% findStimFileBruker.m
%
% With all stim.mat files saved in a stim directory in the date folder and
%  not in the TSeries folders, take in a TSeries folder name and a date
%  folder path and return the full path to the stim.mat file
% Expects a stim folder in the date folder. Stim folder contains numbered
%  folders corresponding to the TSeries numbers
%
% INPUTS:
%   tSeriesName - name of TSeries folder
%   datePath - full path to date folder
%
% OUTPUTS:
%   stimFilepath - full path to stim.mat file for this TSeries
%
% CREATED: 11/14/21 - HHY
%
% UPDATED: 
%   11/16/21 - HHY
%
function stimFilepath = findStimFileBruker(tSeriesName, datePath)
    
    % get index number for this TSeries folder (last 3 characters of name)
    indStr = tSeriesName((end-2):end);
    
    % convert index string to number (removes leading zeros, if any)
    thisInd = str2num(indStr);
    
    % convert index back to string, without leading zeros, for stim folder
    %  names
    stimIndStr = num2str(thisInd);
    
    % full file path for stim.mat
    stimFilepath = [datePath filesep 'stim' filesep stimIndStr filesep ...
        'stim.mat'];
    
end