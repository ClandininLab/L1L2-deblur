% getStim.m
% 
% Reads the first few lines of the stimulus text file provided by the user 
% to determine the appropriate stimulus class to use. Creates an object of
% that stimulus class and calls its setParams() function to set the
% paramter values. 
%
% INPUT:
%   txtFile - string containing the name of the Stimulus .txt file
%
% OUTPUT:
%   s - a Stimulus object (an instantiation of the Stimulus class) 
% 
% Called by playStimMain.m

function s = getStim(txtFile)

% Open the file
fileID = fopen(txtFile);
fgetl(fileID); % line 1 is just the number of params
fgetl(fileID); % line 2 is just the number of epochs
line3 = textscan(fileID, '%s %s', 1); % read line 3 to get stimulus class name
fclose(fileID);

% Convert the stimulus class name to a function. str2func() returns a 
% function handle that can be used to call the constructor method in the 
% stimulus class. The constructor method is what initializes an object of
% that stimulus class.
stimClass = str2func(char(line3{2}));

% Create a stimulus object from the function handle
s = stimClass(txtFile); 
s.setParams;

end 