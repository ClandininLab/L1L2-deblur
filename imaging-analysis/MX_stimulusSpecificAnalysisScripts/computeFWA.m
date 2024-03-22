function [ normFWA, upStimVals, upStimTimes ] = computeFWA(tau, stimSampleRate, ...
    stimVals, stimTimes, respVals, respTimes)
% Computes the fluorescent-weighted average stimulus tau seconds into the
% past. Tip: flip the FWA vector before convolving with the stimulus. 
%
% INPUTS:
%   tau - length of stimulus window in seconds
%   stimSampleRate - stimulus sampling frequency (Hz). Should be at least
%       greater than stimulus contrast update rate
%   stimVals - vector of stimulus contrast values corresponding to 
%       timestamps in stimTimes 
%   stimTimes - vector of timestamps of contrast values in stimVals
%   respVals - vector of dFF values for a single ROI (column vector)
%   respTimes - vector of timestamps for dFF values in respVals
% 
% OUTPUTS:
%   normFWA - fluorescence-weighted average stimulus vector normalized to 
%       a unit vector. First index corresponds to time zero 
%   upStimVals - upsampled stim values based on stimSampleRate

%%
% Step 1: Upsample the stimulus to desired frequency
stimSampleIFI = 1/stimSampleRate;
winLength = tau/stimSampleIFI; % number of frames per stimulus window
stimUpdateTime = mean(diff(stimTimes)); % stim contrast duration
upSampleFactor = round(stimUpdateTime/stimSampleIFI); % number of repetitions per contrast value
upStimVals = repelem(stimVals, upSampleFactor); 
upStimTimes = interp1(stimTimes, 1:(1/upSampleFactor):length(stimTimes));

% debug
% upStimTimes = upStimTimes + 0.05;

% Step 2: Compute FWA stimulus vector of length tau seconds
sumFWS = zeros(winLength, 1); % array to store each FWS
        % vectorized implementation 
        % matrix = zeros(winLength, length(respVals));

% loop over each dF/F value (analogous to spikes in STA)
for w = 1:length(respVals)
    % find index of the front of the stimulus window, aka time 0, aka timestamp
    % of the most recent stimulus transition since this imaging frame
    % + 1 to account for time between stimulus and imaging timestamps
    iFront = find(upStimTimes < respTimes(w), 1, 'last');
    % index of back of the window
    iBack = iFront+1-winLength;
    
        % vectorized implementation 
        %     if (iBack > 0)
        %     matrix(:, w) = upStimVals(iBack:iFront);
        %     end 

    % calculate fluorescence-weighted stimulus (FWS):
    % weight stimulus vector by the response value 
    if (iBack > 0)
        fws = upStimVals(iBack:iFront).*respVals(w);
        sumFWS = sumFWS + fws; % cumulative sum
    end 
end 
        % sumFWS = matrix*respVals;

% Divide FWS by the number of dF/F values to get the FWA 
% (average fluorescence-weighted stimulus vector
FWA = sumFWS./length(respVals);

% Step 3: normalize FWA by its magnitude to shrink to unit vector 
norm = sqrt(sum(FWA.^2)); 
normFWA = FWA./norm;

% Step 4: flip the FWA such that time zero has the first index
% normFWA = flipud(normFWA); % kernel is just the FWA in reverse 

end

