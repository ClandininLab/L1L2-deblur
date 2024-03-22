% test script for stimulus to stimulate screen center

% Clear the workspace and the screen
sca;
close all;
clear all;
clearvars;
clc;

% screen dimensions
screenDim = [130 0 700 1140]; % lower left, upper right verticies
% screenDim = [0 130 1140 700]; % lower left, upper right verticies
pdDim = [30 0 60 60];

% physical location of screen relative to fly (in cm)
scrn_ll = 14.2; % lower left corner
scrn_lr = 12.8; % lower right corner
scrn_ul = 8.7; % upper left corner
scrn_ur = 6.2; % upper right corner
scrnX = 9; % x is vertical
scrnY = 9; % y is horizonal

% conversion between cm and pixels
x_cm2px = (screenDim(3)-screenDim(1))/scrnX;
y_cm2px = (screenDim(4)-screenDim(2))/scrnY;

% compute corners of restricted region

% number of degrees to clip from each edge
edgeDeg = 10;

% list of verticies
verticies = zeros(4,2); % allocate, trapezoid verticies

% amount to clip from each, in pixels
% lower left x 
llx = computeEdgeDistance(edgeDeg, scrnX, scrn_ll, scrn_ul) * x_cm2px;
% lower left y
lly = computeEdgeDistance(edgeDeg, scrnY, scrn_ll, scrn_lr) * y_cm2px;
% lower right x
lrx = computeEdgeDistance(edgeDeg, scrnX, scrn_lr, scrn_ur) * x_cm2px;
% lower right y
lry = computeEdgeDistance(edgeDeg, scrnY, scrn_lr, scrn_ll) * y_cm2px;
% upper left x
ulx = computeEdgeDistance(edgeDeg, scrnX, scrn_ul, scrn_ll) * x_cm2px;
% upper left y
uly = computeEdgeDistance(edgeDeg, scrnY, scrn_ul, scrn_ur) * y_cm2px;
% upper right x
urx = computeEdgeDistance(edgeDeg, scrnX, scrn_ur, scrn_lr) * x_cm2px;
% upper right y
ury = computeEdgeDistance(edgeDeg, scrnY, scrn_ur, scrn_ul) * y_cm2px;

% lower left, lower right, upper left, upper right)
verticies(1,:) = [screenDim(1) + llx, screenDim(2) + lly];
verticies(2,:) = [screenDim(1) + lrx, screenDim(4) - lry];
verticies(4,:) = [screenDim(3) - ulx, screenDim(2) + uly];
verticies(3,:) = [screenDim(3) - urx, screenDim(4) - ury];


% set up Psychtoolbox
PsychDefaultSetup(2);
screens = Screen('Screens');
Screen('Preference', 'VisualDebugLevel',1);
screenNumber = max(screens);

% define contrast values
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);

% open window, color black
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
Screen('ColorRange',window,[],-1,[]);

% to test, draw white box for screen
Screen('FillRect', window, [1 1 1], screenDim);
% draw white box for photodiode
Screen('FillRect', window, [1 1 1], pdDim);
% draw black box for screen center
Screen('FillPoly',window, [0 0 0], verticies,1);

Screen('Flip', window);

KbStrokeWait;

sca;


