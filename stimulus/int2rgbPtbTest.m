% int2rgbPtbTest.m
%
% Displays a static test pattern using psychtoolbox to check that 
%  intensity values are displaying correctly when presenting different 
%  values for different patterns to increase effective frame rate of 
%  lightcrafter.
% Lightcrafter
%
% INPUT: 
%  screenNumber - screen number corresponding to lightcrafter - check what
%   is displayed in PTB-INFO when this is run for the first time and what
%   Windows names the screen
%
function int2rgbPtbTest(screenNumber)
    % Clear the workspace and the screen
    sca;
    close all;
    clear functions;
    clearvars -except screenNumber
    clc;

    % Here we call some default settings for setting up Psychtoolbox
    PsychDefaultSetup(2);

    % Get the screen numbers. This gives us a number for each of the screens
    % attached to our computer.
    % For help see: Screen Screens?
%     screens = Screen('Screens');
    Screen('Preference', 'VisualDebugLevel',1);

    % Define black (white will be 1 and black 0). This is because
    % luminace values are (in general) defined between 0 and 1.
    % For help see: help BlackIndex
    black = BlackIndex(screenNumber);
    white = WhiteIndex(screenNumber);

    % Open an on screen window and color it black
    % For help see: Screen OpenWindow?
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
    % Screen('ColorRange',window,[],-1,[]);
    
    % screen
    centeredRect = [130 0 700 1140];
    pdRect = [30 0 60 60];
    pdCol = [1 1 1];
    
    % get coordinates of all rectangles
    numCols = 64; % number of colors, equal to num rect hz
    numRectVt = 7;
    numPxHz = centeredRect(3) - centeredRect(1);
    numPxVt = centeredRect(4) - centeredRect(2);
    rectSizeHz = floor(numPxHz/numCols);
    rectSizeVt = floor(numPxVt/numRectVt);
    hzStart = round((numPxHz - (rectSizeHz * numCols))/2) + ...
        centeredRect(1);
    vtStart = round((numPxVt - (rectSizeVt * numRectVt))/2) + ...
        centeredRect(2);
    hzRectLeft = hzStart:rectSizeHz:(hzStart + rectSizeHz*(numCols-1));
    vtRectBottom = vtStart:rectSizeVt:(vtStart + rectSizeVt*(numRectVt-1));
    hzRectRight = (hzStart+rectSizeHz):rectSizeHz:(hzStart + ...
        rectSizeHz*numCols);
    vtRectTop = (vtStart+rectSizeVt):rectSizeVt:(vtStart + ...
        rectSizeVt*numRectVt);    
    
    hzRectLeftAll = repmat(hzRectLeft, 1, numRectVt);
    vtRectBottomAll = repelem(vtRectBottom, 1, numCols);
    hzRectRightAll = repmat(hzRectRight, 1, numRectVt);
    vtRectTopAll = repelem(vtRectTop, 1, numCols);
    
    rectCoords = [hzRectLeftAll; vtRectBottomAll; hzRectRightAll; ...
        vtRectTopAll];
    
    % get colors of all rectangles
    % intensity sequence for hz row 1 (only pattern 1)
    intVal = 0:(1/(numCols-1)):1; % increasing values 0 to 1, 64 diff vals
    row1Seq = upsample(intVal,3);
    % intensity sequence for hz row 2 (only pattern 2)
    row2Seq = circshift(row1Seq,1);
    % intensity sequence for hz row 3 (only pattern 3)
    row3Seq = circshift(row1Seq,2);
    % intensity sequence for hz row 4 (patterns 1 and 2 only)
    row4Seq = row1Seq + row2Seq;
    % intensity sequence for hz row 5 (patterns 2 and 3 only)
    row5Seq = row2Seq + row3Seq;
    % intensity sequence for hz row 6 (patterns 1 and 3 only)
    row6Seq = row1Seq + row3Seq;
    % intensity sequence for hz row 7 (patterns 1, 2, and 3)
    row7Seq = row1Seq + row2Seq + row3Seq;
    
    % concatenate into 1 vector
    intSeqAll = [row1Seq row2Seq row3Seq row4Seq row5Seq row6Seq row7Seq];
    % convert intensity sequence to RGB sequence
    rgbAll = convertInt2RGB(intSeqAll);
    % transpose for FillRect function of PtB
    rgbAll = rgbAll';
    
    
    Screen('FillRect', window, rgbAll, rectCoords);
    Screen('FillRect', window, pdCol, pdRect);

    Screen('Flip', window);

    % Now we have drawn to the screen we wait for a keyboard button press 
    % (any key) to terminate the demo.
    % For help see: help KbStrokeWait
    KbStrokeWait;

    % Clear the screen. "sca" is short hand for "Screen CloseAll". This 
    % clears all features related to PTB. Note: we leave the variables in
    % the workspace so you can have a look at them if you want.
    % For help see: help sca
    sca;

end