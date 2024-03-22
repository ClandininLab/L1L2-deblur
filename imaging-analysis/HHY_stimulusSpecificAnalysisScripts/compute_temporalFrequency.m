% compute_temporalFrequency.m
%
% Script to load impulse responses from full field flash onto gray saved 
%  data and compute temporal frequency representation of each cell's
%  impulse response
%
% Updated: HHY 9/27/17
%

savePath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170927/';

% get directory
dirName = uigetdir('/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/');

cd(dirName);

files = dir;

% run processing on every .mat file in directory
for i = 3:length(files)
    % load .mat file
    load(files(i).name);
    
    % compute temporal frequency representation on each ROI
    for j = 1:size(roiDataMat,1)
        % compute for each impulse response in rats
        for k = 1:length(roiDataMat(j,2).rats)
            timeTrace = roiDataMat(j,2).rats{k};
%             nSamp = 2^nextpow2(length(timeTrace)); 
            nSamp = 1024;
            % normalize impulse response
%             timeTrace = timeTrace / max(abs(timeTrace));
            % smooth impulse response (Savitzky-Golay filter)
            smoothTrace = smooth(timeTrace',7,'sgolay',5);
            
            % perform fft
            y = fft(smoothTrace,nSamp);
            % compute power
%             power = abs(y).^2/nSamp;
%             power = abs(y/nSamp);
            % compute magitude
            magnitude = abs(y);
            % compute phase
            phase = unwrap(angle(y));
            % temporal frequencies of fft
            tempFreq = (0:length(y)-1)*interpFrameRate/length(y);
            
            % indicies of positive freq half
            p1 = 1:(floor(length(y)/2)+1);
            
            % only the positive frequency half of the frequency spectrum
            posMagnitude = magnitude(p1);
            posPhase = phase(p1);
            posTempFreq = tempFreq(p1);
            
%             figure; loglog(tempFreq,posPowerSpec)

            % save temp freq power spectrum into roiDataMat
            roiDataMat(j,2).tfMag{k} = posMagnitude;
            roiDataMat(j,2).tfPhase{k} = posPhase;
            roiDataMat(j,2).tempFreq{k} = posTempFreq;
            
        end
    end
    
    % resave .mat file
    save([savePath filesep files(i).name], 'roiDataMat', 'roiMetaMat', ...
        'iResp', 'iInv','binWidthMult','interpFrameRate','inv',...
        'yScale', '-v7.3');

end
