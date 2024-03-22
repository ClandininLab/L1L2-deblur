% preprocess.m
%
% First function to run on imaging data. 
% Accepts a .lif file. Reads each imaging time series frame by frame. 
% Aligns each time series separately. Saves the aligned series, 
% unaligned series, mean image, and field of view dimensions in a struct 
% in the fly folder subdirectory that contains the stimulus metadata. 
%
% NOTE: assumes folder organization of date folder with one .lif file and
%   date folder has fly folders that have time series folders
%
% Last update: 3/8/17 HHY
%
function preprocess()

    numRefFrames = 30; % default to 30 reference frames for alignment

    display('Select the date folder that contains your .lif file.');
    datepath = uigetdir;
    curDir = pwd;
    cd(datepath);
    lifFile = dir('*.lif');

    % Initialize logging
    bfInitLogging();

    % initialize bioformat reader
    r = bfGetReader(lifFile.name); % assumes there is only one .lif file in the date directory
    nSeries = r.getSeriesCount(); % index starts from 0

    % loop through each time series in the .lif file
    for s = 1:nSeries
        r.setSeries(s - 1); % tells reader to extract the (s-1)th series
        totFrames = r.getImageCount();

        if totFrames > numRefFrames % series must have enough frames to align   
            % Retrieve metadata about this time series)
            meta = r.getMetadataStore();
            label = char(meta.getImageName(s - 1));
            % if LDM part of the series name is followed by a '/', keep just the
            % LDM part
            if ~isempty(regexp(label, '/', 'ONCE'))
                seriesName = label(1:regexp(label, '/') - 1);
            else 
                seriesName = label;
            end 
            
            xDim = meta.getPixelsSizeX(s - 1).getValue();
            yDim = meta.getPixelsSizeY(s - 1).getValue();
        %     lineScanRate = meta.getLaserRepetitionRate(int instrumentIndex, int lightSourceIndex);

            fprintf('\n Reading series %d \n', s - 1);

            % allocate memory
            unalignedSeries = zeros(yDim, xDim, totFrames);
            % Read the time series frame by frame
            for iFrame = 1:totFrames
                % read each frame of time series one by one
                unalignedSeries(:,:,iFrame) = uint16(bfGetPlane(r, iFrame)); 
                % class uint16 by default
                % commented out to stop printing of all the .'s; slows 
                %  preprocess down
        %         if mod(iFrame, 72) == 1
        %             fprintf('\n    ');
        %         end
        %         fprintf('.');
            end 

            % Create reference frame for image alignment 
            % reference stack using first numRefFrames frames 
            refStack = unalignedSeries(:,:,1:numRefFrames); 
            % max intensity projection of ref stack.
            refFrame = max(refStack, [], 3); 

            % Align the time series
            alignedSeries = fccAlignment(unalignedSeries, refFrame, 'xml');
            alignedSeries = uint16(alignedSeries); % save into cell array
            fprintf('%s alignment completed! \n', seriesName);

            % Average over time of the whole aligned time series
            meanImageAligned = mean(alignedSeries, 3);

        %     % visualize mean image
        %     figure; imshow(meanImageAligned, []); title('mean aligned');
        %     figure; imagesc((meanImageAligned)); colormap('gray'); title('mean aligned');
        %     figure; imagesc((mean(unalignedSeries,3))); colormap('gray'); title('mean unaligned');
        %     figure; imagesc((refFrame)); colormap('gray'); title('refFrame');

            saveSeriesImData(seriesName, alignedSeries, unalignedSeries, ...
                meanImageAligned, xDim, yDim, datepath);

            fprintf('Saved! \n');

            clear unalignedSeries alignedSeries refStack refFrame meanImageAligned
        else
            fprintf('\n Skipping series %d. Less than %d frames.  \n', ...
                s - 1, numRefFrames);
        end

    end 

    r.close(); % close the reader

end