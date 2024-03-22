% White noise troubleshooting:
% WhiteNoise_1D_FullFieldFlash_playback.m takes a picture of the displayed
% stimulus at each frame, which is saved as stimMovie. Here, we compare 
% stimMovie to frameStimVal (the stimulus sequence generated prior to the
% display for loop) and rawStim (the stimulus values saved during the
% display period).

% play the stimulus back
stim = obj;
imageArray = obj.Out.stimMovie;
N = size(imageArray, 3);

%% Play PTB's screenshots of the stimulus as a movie
figure;
for f = 1:N;
    imshow(imageArray(1:50, 1:50, f)./63);
end 

%% Plot the stimulus
contrastArray = reshape(imageArray(1,1,:), N,1);
contrastArray = contrastArray./63;

figure;
hist(contrastArray, 20); % sanity check to see distribution of contrasts

%% use if generated stimulus prior to display
figure;
plot(contrastArray); hold on;
plot(stim.Out.rawStim, 'y', 'LineWidth', 2); 
plot(stim.Out.frameStimVals, 'r');
legend('PTB screenshots', 'saved stimulus', 'stimulus generated prior');
ylim([-0.5 1.5]);

%% use if stimulus is generate live
figure;
plot(contrastArray); hold on;
plot(stim.Out.rawStim);
legend('PTB screenshots', 'saved stimulus');
ylim([-0.5 1.5]);