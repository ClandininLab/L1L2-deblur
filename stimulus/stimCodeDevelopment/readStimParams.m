function stim = readStimParams(sFileName)
% Reads stimulus parameters from a .txt file and saves into a structure
% named 'stim'
% Called by playstim_v2.m

% stim.type = cell(1,7);
% stim.contrast = cell(1,7);
stim.funcName = '';

% Open the file
fileID = fopen(sFileName);
line1 = textscan(fileID, '%s %f', 1); % read line 1 to get number of params
line2 = textscan(fileID, '%s %f', 1); % read line 2 to get number of epochs
line3 = textscan(fileID, '%s %s', 1); % read line 3 to get stimulus function name

numParams = line1{2}; 
numEpochs = line2{2};
stim.funcName = line3{2};

% Read the rest of the file to store values across epochs
% for the various parameters
% First loop through each parameter, reading the file line by line
for p = 1:numParams-1
    paramName = textscan(fileID,'%s', 1); % name of paramter
    values = textscan(fileID,'%f %f %f %f %f %f %f', 1); % parameter values
    
%     % Loop through all the epochs and save the values into the correct
%     % parameter field
%     for epoch = 1:numEpochs
%         if strcmp(paramName{1}, 'Stimulus.stimtype');
%             stim.type{epoch} = values{epoch};
%         elseif strcmp(paramName{1}, 'Stimulus.contrast');
%             stim.contrast{epoch} = values{epoch};
%         elseif strcmp(paramName{1}, 'Stimulus.duration');
%             stim.duration{epoch} = values{epoch};
%         elseif strcmp(paramName{1}, 'Stimulus.contrast');
%             stim.contrast{epoch} = values{epoch};
%         end 
%     end 

% % faster way to save all the parameters:
%     if strcmp(paramName{1}, 'Stimulus.stimtype');
%            stim.type = cell(1,numEpochs);
%            stim.type(1:numEpochs) = values(1:numEpochs);
%     elseif strcmp(paramName{1}, 'Stimulus.duration');
%            stim.duration = cell(1,numEpochs);
%            stim.duration(1:numEpochs) = values(1:numEpochs);
%     elseif strcmp(paramName{1}, 'Stimulus.contrast');
%            stim.contrast = cell(1,numEpochs);
%            stim.contrast(1:numEpochs) = values(1:numEpochs);
%     end 
    
% even faster way  
   pname = char(paramName{1});
   field = matlab.lang.makeValidName(pname);
   stim.(field) = cell(1,numEpochs);
   stim.(field)(1:numEpochs) = values(1:numEpochs);


    
end 

fclose(fileID);

% % good to know: how to treat each column in the text file as a vector
% fileID = fopen('stimparams.txt');
% names = textscan(fileID,'%s %s %s %s %s %s %s %s');
% fclose(fileID);
end 

