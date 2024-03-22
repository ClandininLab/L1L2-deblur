% To Do

% 1. frame loop or stimulus update loop on the outside? 
%       WhiteNoise_1D_FullFieldFlash.m has a the frame loop on the outside
%       (method  A)
%       ShortFlashOntoGray_ files currently do the stimulus update loop on
%       the outside

% 2. How do I keep track of the the actual stimulus played into struct Out?
%           save the stimulus as you go. export to Out
%
% 3. Make saveData() take variable number of inputs 

% 4. For ShortFlashOntoGray, the epochs may be different lengths and when
% chosen randomly would make the total number of epochs estimate imprecise.
% should we draw a random contrast value as we play the sitmulus or set 
% the entire random sequence before?

% 5. How should we generate random values? In WhiteNoise_1d_FullFieldFlash,
% I implemented both and commented out one of them. Current method is draw 
% a random number each time the stimulus is updated 

% 6. For whitenoise, should the durations be the same for each epoch? can we
% just make it simple and have NumContrastValues, EpochDuration, and
% Totalduration, and assume that we are spacing the contrast values equally
% and equally in time? I just made it one duration value.