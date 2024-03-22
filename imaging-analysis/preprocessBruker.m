% preprocessBruker.m
%
% First function to run on imaging data collected on the Bruker scope.
%
% Accepts date folder. Reads each imaging time series (tiff file) frame by 
%  frame. 
% Aligns each time series separately. Saves the aligned series, 
%  unaligned series, mean image, and field of view dimensions in a struct 
%  in the fly folder subdirectory that contains the stimulus metadata. 
%
% NOTE: assumes folder organization of date folder a series of TSeries
% folders, each corresponding to one trial and containing a .tiff file for
% the 
%
% Last update: 11/16/21 - HHY
%

function preprocessBruker()

    numRefFrames = 30; % default to 30 reference frames for alignment

    disp('Select a date folder to preprocess');
    datepath = uigetdir;
    curDir = pwd;
    cd(datepath);
    
    % get TSeries folder names and paths
    tSeriesFolders = dir('TSeries*');
    
    % number of TSeries folder
    nSeries = length(tSeriesFolders);
    
    % loop through all TSeries folders, read in metadata (frame times,
    %  xDim, yDim), image files
    for s = 1:nSeries
        
        % name for this TSeries
        tSeriesName = tSeriesFolders(s).name;
        
        % full path for this TSeries folder
        tSeriesFolderPath = [tSeriesFolders(s).folder filesep ...
            tSeriesFolders(s).name];
        
        % full path of xml file with frame times, image dimensions
        % assume it has the same name as the TSeries folder it is in
        xmlFilePath = [tSeriesFolderPath filesep tSeriesName '.xml'];
        
        fprintf('Reading XML file for %s\n', tSeriesName);
        
        % get image dimensions and frame times from xml file
        [frameTimes, xDim, yDim] = parseBrukerXML(xmlFilePath);
        
        
        % read in image file for this TSeries
        
        % total number of frames, from xml file
        totFrames = length(frameTimes);
        
        % image file full path, assume file name is 'stack.tiff'
        imgFilepath = [tSeriesFolderPath filesep 'stack.tiff'];
        
%         % get info for tiff file
%         imgInfo = imfinfo(imgFilepath);
        
%         % tiff file as Tiff Obj - note, this is much faster than imread
        tiffObj = Tiff(imgFilepath,'r');

        % use TIFFStack class to read in time series - 2/16/22

        % create TIFFStack object
%         tsStack = TIFFStack(imgFilepath);
        
        % series must have enough frames to align
        if (totFrames > numRefFrames) 
            fprintf('\n Reading series %s \n', tSeriesName);
            
            % allocate memory, uint16 format (for 13 bit images)
            unalignedSeries = uint16(zeros(yDim, xDim, totFrames));
            
            % read in time series, frame at a time
            for i = 1:totFrames
                % imread is way too slow
%                 unalignedSeries(:,:,i) = imread(imgFilepath, ...
%                     'Index', i, 'Info', imgInfo);

                % reading using Tiff Object
                unalignedSeries(:,:,i) = read(tiffObj);

                % go to next image, but only if there is another image
                if (i ~= totFrames)
                    tiffObj.nextDirectory;
                end

                % reading using TIFFStack object
%                 unalignedSeries(:,:,i) = tsStack(:,:,i);

                
                if ~(mod(i,1000))
                    fprintf('Reading img: %d \n',i);
                end
            end

%             % read in tiff using imread_big - doesn't have 2^16 overflow
%             % issue of Tiff Obj - 1/5/22
%             % doesn't seem to read time series appropriately 2/4/22 - data
%             %  looks like noise
%             [unalignedSeries, nFrames] = imread_big(imgFilepath);

            
            % 1/5/22 note: nFrames != totFrames, unclear why - figure out
            % Bruker output
            % 2/4/22 - figured out - issue with python code generating time
            % series
            
            % Create reference frame for image alignment 
            % reference stack using first numRefFrames frames 
            refStack = unalignedSeries(:,:,1:numRefFrames); 
            % max intensity projection of ref stack.
            refFrame = max(refStack, [], 3); 

            fprintf('Aligning series %s \n', tSeriesName);
            
            % Align the time series
            alignedSeries = fccAlignment(unalignedSeries, refFrame, 'xml');
            alignedSeries = uint16(alignedSeries); % save into cell array
            fprintf('%s alignment completed! \n', tSeriesName);

            % Average over time of the whole aligned time series
            meanImageAligned = mean(alignedSeries, 3);
            
            % save info into imDat.mat
            saveBrukerImgSeriesData(tSeriesName, alignedSeries, ...
                unalignedSeries, meanImageAligned, xDim, yDim, ...
                frameTimes,  tSeriesFolderPath);
            
            disp('Saved!');
            
            clear unalignedSeries alignedSeries refStack refFrame meanImageAligned

        % skip time series with too few images
        else
            fprintf('\n Skipping series %s. Less than %d frames.  \n', ...
                tSeriesName, numRefFrames);
        end
                
    end
end