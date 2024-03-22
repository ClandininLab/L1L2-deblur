% PerspectiveCorrectionTest_v1.m
%
% Function to test perspective correction algorithm: distorts in both
%  horizontal and vertical dimensions according to angle spanned, linear
%  interpolation between corners to distort central area
%
% Mapping arbitrary quadrilateral to unit square. Get coordinates. Then use
%  coordinates to extract intensity values to generate texture.
%
% Mapping done with bilinear mapping function. See reference:
%  https://www.particleincell.com/2012/quad-interpolation/
%
% Sampling done with (look up interpolation methods in Matlab)
%
% Last updated: 1/15/17
%

function PerspectiveCorrectionTest_v1
%% define parameters
    % screen dimensions
    screenDim = [130 0 700 1140]; % lower left, upper right verticies
    % screenDim = [0 130 1140 700]; % rotated 90 deg, for horiz screen
%     screenDim = [130 0 415 570]; % rescaled 1/2 but aspect ratio to match lightcrafter
%     screenDim = [130 200 700 485]; % rescaled 1/2, rotated 90 deg
    pdDim = [30 0 60 60];

    % physical location of screen relative to fly (in cm)
    scrn_ll = 14.2; % lower left corner
    scrn_lr = 12.8; % lower right corner
    scrn_ul = 8.7; % upper left corner
    scrn_ur = 6.2; % upper right corner
    scrnX = 9; % x is vertical
    scrnY = 9; % y is horizonal
    
    % fly horizontally out to vertical plane of screen, left and right
    %  sides
    horiz_lside = 8.2;
    horiz_rside = 5.7;
    % vertical distance between horizontal plane of fly and screen, left
    %  and right sides
    vert_lside = 2.5;
    vert_rside = 2.5;
    
    % texture parameters
    txtrSize = 90; % texture size in degrees; square
    txtrRes = 2; % texture resolution, in pixels per degree
    txtrX = txtrSize - 25; % upper right corner of texture, x
    txtrY = txtrSize - 10; % upper right corner of texture, y
 
%% compute static measurements for mapping texture to screen
    % compute vertices of persepctive quadrilateral
    [urx, ury, ulx, uly, llx, lly, lrx, lry] = ...
        computePerspecQuadVerticies(scrn_ul, scrn_ur, scrn_ll, scrn_lr, ...
        scrnX, scrnY, horiz_lside, horiz_rside, vert_lside, vert_rside, ...
        txtrX,txtrY);
    
    % temp test plot of quadrilateral; 11/27/16, looks good
%     figure; hold on;
%     plot([0 0 90 90 0],[0 90 90 0 0],'b'); % texture
%     plot([lly lry ury uly lly],[llx lrx urx ulx llx],'r'); % quadrilateral
    
    % convert verticies of quadrilateral to coefficients for mapping
    %  quadrilateral to unit square
    xVert = [llx; lrx; urx; ulx];
    yVert = [lly; lry; ury; uly];
    [xCoeffs, yCoeffs] = mapPerspecQuadVertToUnitSquare(xVert,yVert);
    
    % temp test of LMtoXY function, using coefficients calculated in
    %  previous section
    % generate matrix of l and m values 11x11, count from 0 to 1 in 0.1
    %  step increments 
%     l = repmat((0:0.1:1)',1,11);
%     m = repmat((0:0.1:1),11,1);
%     [m,l]=meshgrid((0:0.1:1),(0:0.1:1));
%     
%     [quadX, quadY] = LMtoXY(xCoeffs,yCoeffs,l,m);
%     
%     % test plot: all generated (x,y) values should fall within
%     %  quadrilateral
%     figure; hold on;
%     plot([0 0 90 90 0],[0 90 90 0 0],'b'); % texture
%     plot([lly lry ury uly lly],[llx lrx urx ulx llx],'r'); % quadrilateral
%     plot(quadY,quadX,'kx'); % plot all test points as black x's
    
    % Use the coefficients to get (x,y) coordinates for getting intensity 
    %  values of texture
    
    screenXDim = screenDim(3) - screenDim(1);
    screenYDim = screenDim(4) - screenDim(2);
    
    % get (l,m) at resolution of txtrRes
    numBinsX = screenXDim / txtrRes;
    numBinsY = screenYDim / txtrRes;
    
    % compute l and m: evenly spaced between 0 and 1
%     l = repmat((0:(1/(numBinsX-1)):1)',1,numBinsY);
%     m = repmat((0:(1/(numBinsY-1)):1),numBinsX,1);
    [m,l] = meshgrid((0:(1/(numBinsY-1)):1),(0:(1/(numBinsX-1)):1));
    
    % get (x,y) corresponding to (l,m)
    [quadX,quadY] = LMtoXY(xCoeffs,yCoeffs,l,m);
    
%     figure; hold on;
%     plot([0 0 90 90 0],[0 90 90 0 0],'b'); % texture
%     plot([lly lry ury uly lly],[llx lrx urx ulx llx],'r'); % quadrilateral
%     plot(quadY,quadX,'kx'); % plot all test points as black x's

%% use coordinates to get texture values
    % test texture - define matrix at resolution txtrRes for
    %  txtrSize x txtrSize texture
    % 10 deg grating
%     testTxtr = [repmat([ones(10*txtrRes,txtrSize*txtrRes); ...
%         zeros(10*txtrRes,txtrSize*txtrRes)],4,1); ...
%         ones(10*txtrRes,txtrSize*txtrRes)];
    
    % test texture - 10 deg white square somewhere 
%     testTxtr = zeros(txtrSize*txtrRes,txtrSize*txtrRes);
%     testTxtr(81:end,61:70) = 1; 

    % test texture - 10 deg grating, rotated x degrees
    startTxtrSize = ceil(2*txtrSize*txtrRes*sind(45));
    startTestTxtr = [repmat([ones(10*txtrRes,startTxtrSize*txtrRes); ...
        zeros(10*txtrRes,startTxtrSize*txtrRes)],6,1); ...
        ones(8*txtrRes,startTxtrSize*txtrRes)]; 
    % rotate 
    degRot = 45; % 45 deg CCW rotation
    rotTestTxtr = imrotate(startTestTxtr,degRot,'bilinear');
    % crop back to txtrSize
    rotTxtrSize = size(rotTestTxtr);
    startX = floor(rotTxtrSize(1)/2 - (txtrSize*txtrRes)/2);
    endX = startX + (txtrSize*txtrRes) - 1;
    startY = floor(rotTxtrSize(2)/2 - (txtrSize*txtrRes)/2);
    endY = startY + (txtrSize*txtrRes) - 1;
    testTxtr = rotTestTxtr(startX:endX,startY:endY);
%     testTxtr = imcrop(rotTestTxtr,...
%         [startX, startY, (txtrSize*txtrRes)-1, (txtrSize*txtrRes)-1]);
    
    % coordinates of test texture 
    [txtrX, txtrY] = meshgrid((1/txtrRes):(1/txtrRes):txtrSize,...
        (1/txtrRes):(1/txtrRes):txtrSize);
    
    % use interp2 function (probably) to convert coordiates and texture to
    %  perspective corrected texture
    
    pcTxtr = interp2(txtrX,txtrY,testTxtr,quadX,quadY);

%% psychtoolbox

    % set up Psychtoolbox
    PsychDefaultSetup(2);
    Screen('Preference', 'VisualDebugLevel',1);
    screenNumber = 1;

    % define contrast values
    black = BlackIndex(screenNumber);
    val = 63/255;
    colScreen = [0 val 0];

    % open window, color black
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
    Screen('ColorRange',window,[],-1,[]);
    
    txtrInd = Screen('MakeTexture',window,pcTxtr*val);
    Screen('DrawTexture',window,txtrInd,[],screenDim);
    
    % draw white box for photodiode
    Screen('FillRect', window, colScreen, pdDim);
    
    Screen('Flip', window);
    
    
    
    KbStrokeWait;

    sca;
end