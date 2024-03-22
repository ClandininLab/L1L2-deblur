% saveSeriesImData.m
%
% Saves times series data in the correct fly folder subdirectory 
%
% INPUT:
%   seriesName - name of time series
%   alignedSeries - image sequence, aligned time series
%   unalignedSeries - image sequence, original, unaligned time series
%   meanImageAligned - mean over time of the whole image sequence
%   xDim - size of imaging frame, x dimension
%   yDim - size of imaging frame, y dimension
%   datepath - path to folder that contains .lif file
%
% OUTPUT:
%   no output, saves time series data as imDat.mat file into correct folder
%

function saveSeriesImData(seriesName, alignedSeries, unalignedSeries, ...
        meanImageAligned, xDim, yDim, datepath)

    % open date folder
    cd(datepath);
    flyFolders = dir(datepath);
    % remove files or nondirectories
    isub = [flyFolders(:).isdir]; % indices of subdirectories 
    flyFolders = flyFolders(isub);
    % remove '.' and '..'
    nameFolds = {flyFolders.name}'; 
    flyFolders(ismember(nameFolds,{'.','..'})) = [];

    % Loop through fly folders. For each subdirectory, 
    % search the labels of one of the frames of 
    % each series for LDM_<number indicated in the stim folder>
    % save the series in the correct stimulus folder 
    for i = 1:length(flyFolders)
        flypath = [datepath filesep flyFolders(i).name];
        % go into the fly folder
        cd(flypath);
        % get list of folders and files from this fly folder 
        subDir = dir(flypath);
        % remove files or nondirectories
        isub = [subDir(:).isdir]; % indices of subdirectories 
        subDir = subDir(isub);
        % remove '.' and '..'
        nameFolds = {subDir.name}'; 
        subDir(ismember(nameFolds,{'.','..'})) = [];

        % For each LDM folder in the fly folder 
        % open stim.mat and look for the imaging folder LDM name in the struct
        % for each time series, check whether the label contains that LDM 
        % if so, save that time series in the LDM folder 
        for j = 1:length(subDir)
    %         j % debug
            seriespath = [flypath, filesep, subDir(j).name];
            cd(seriespath);
            stim = load('stim.mat');
            sName = stim.obj.Out.imDataName;

            % if the times series name in stim.mat matches this series name
            if (strcmp(sName, seriesName))
                 fprintf('Saving image data for %s ... \n', seriesName);

                % save image data for this time series in a struct
                imDat = struct; 
                imDat.fileloc = seriespath; % file location
                imDat.flyName = flyFolders(i).name;
                imDat.seriesName = seriesName; % series name
                imDat.alSeries = alignedSeries; % aligned time series
                imDat.unSeries = unalignedSeries; % unaligned series
                imDat.avgIm = meanImageAligned; % average image
                imDat.nFrames = size(alignedSeries, 3);
                imDat.xPixels = xDim; % pixels in x dimension of FOV
                imDat.yPixels = yDim; % pixels in y dimension of FOV

                save('imDat.mat', '-struct', 'imDat');
                cd(datepath);
                return; 
            end
        end 
    end 
    % if we've searched the last LDM folder and still haven't found
    % the right stim.mat file, print error message
    fprintf('No stim.mat file found for time series %s \n', seriesName); 
end 