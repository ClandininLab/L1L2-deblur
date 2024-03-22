% testStimulusDisplay.m
% 
% Tests visual stimulus display locally. Make sure Psychtoolbox is
% installed on your computer.
%

% close all;
clear all;

addpath('C:\Users\Marjo\Documents\GitHub\2p-stim-code');

% set the stimulus duration for testing 
stimDuration = 6; % seconds

% % Make an object out of one of the stimulus subclasses and a .txt file
% s = FullFieldFlashSubclass('2sFFF_testsub.txt');
% s = FullFieldFlashTestSplitWindow('2sFFF_testsplitwindow.txt');
% s = FullFieldFlash_BinaryContrast('2sFFF_template.txt');
% s = SearchStimulusFlash('2s_searchStimFlash.txt');
% s = FullFieldFlash_NEpochs('FFF_nepochs_template.txt');

% s = ShortFlashOntoGray_test1('FFF_shortflash_onto_gray.txt');
% s = ShortFlashOntoGray_test2('FFF_shortflash_onto_gray.txt');
% s = ShortFlashOntoGray_test3('FFF_shortflash_onto_gray.txt');
% s = FullFieldFlashOntoGray('FFF_onto_gray_template.txt');
% s = WhiteNoise_1D_FullFieldFlash('whitenoise_fullfieldflash_template.txt');
% s = WhiteNoise_1D_FullFieldFlash('whitenoise_1D_timingtest.txt');
s = WhiteNoise_1D_FullFieldFlash_testflip('whitenoise_1D_rand9_20ms_testflip.txt');

% load stimulus parameters from the .txt file into the stimulus object
s.setParams;

% initialize Psychtoolbox
[window, windowRect, screenX, screenY, ifi] =  initPsychTbx();

% display the stimulus 
s.displayStim(window, ifi, stimDuration);

% close all on the screen
sca; 

% plot stimulus timing info
figure; plot(diff(s.Out.flipStart));
ylabegitl('difference between start of screen flip times (sec)');











