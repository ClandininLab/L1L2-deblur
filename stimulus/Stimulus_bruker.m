% Stimulus is an abstract class for stimulus presentation that provides 
%  methods used by all of its subclasses (e.g. FullFieldFlash) for 
%  stimulus presentation, saving timing info, and subsequent analysis.
% It takes in a text file of stimulus parameters as input. 
%
% See Matlab documentation on object-oriented programming and abstract
% classes:
%   https://www.mathworks.com/help/matlab/object-oriented-programming.html
%   https://www.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html
%
% Constructor: obj = Stimulus(txtFile)
%
% Abstract Properties
%   ParamList - list of stimulus parameters
%   Out - stimulus output
%
% Properties
%   TxtFileName - name of text file that contains all the parameters that
%       define the stimulus
%   nParams - number of parameters defined in text file
%   nEpochs - number of different epochs defined in text file, where each
%       epoch has different parameters
%   IFI - interframe interval of visual stimulus
%   NIDAQScanRate - NIDAQ scans per second (i.e. data points per second)
%   screenDim - dimensions and location of screen, defined as top L and
%       bottom R coordinates 
%   pdDim - dimensions and location of photodiode, defined as top L and
%       bottom R coordinates
%   physDist - struct of physical distances bewteen fly and screen
%
% Abstract Methods
%   displayStim - displays stimulus using psychtoolbox
%   reconstructStim - reconstructs stimulus for analysis
%   selectDFFbaseline - selects imaging frames to be used for F0 in
%       computing dF/F
%
% Methods
%   exportData - saves data as .mat file
%   saveData - record into this object the photodiode output and imaging
%       frame time data (as well as timing for both)
%   setParams - reads .txt file of parameters and imports them into the
%       object
%   setNIDAQScanRate - saves NIDAQScanRate
%
% Updates:
%   3/7/17 - added reconstructStim abstract method
%   3/8/17 - added selectDFFbaseline abstract method
%   8/5/21 - changed physical distances for use with Bruker; change file
%   name to Stimulus.m when using (and change orig stim file to other name)
%
% Last update: 8/5/21 MMP
%

classdef (Abstract) Stimulus < handle
    % properties implemented in the subclass 
    properties (Abstract)
        ParamList
        Out % stimulus output
    end 
    % properties implemented here
    properties
        TxtFileName
        nParams
        nEpochs   
        IFI % interframe interval
        NIDAQScanRate % NIDAQ scans per second 
        screenDim = [230 200 700 1140]; % top L & bottom R coordinates
        pdDim = [50 980 140 1160]; % photodiode coordinates
        RGBscale = 63/255; % correction for 6-bit color coding
        
        % struct that defines physical distances (in cm)
        physDist = struct('fly2scrnLR', 11.7,... % fly to screen lower left corner
            'fly2scrnUR', 11.4,... % fly to screen lower right corner
            'fly2scrnLL', 6.6,... % fly to screen upper left corner
            'fly2scrnUL', 6.1,... % fly to screen upper right corner
            'fly2horizL', 10.2,... % fly horizontally to vertical plane of screen, left side
            'fly2horizR', 10.6,... % fly horizontally to vertical plane of screen, right side
            'fly2topL', -4,... % horizontal plane of fly to screen top, vertically, on left
            'fly2topR', -4,... % horizontal plane of fly to screen top, vertically, on right
            'scrnX', 9,... % screen side length, x is vertical
            'scrnY', 9); % screen side length, y is horizonal
%         physDist = struct('fly2scrnLL', 11.7,... % fly to screen lower left corner
%             'fly2scrnLR', 11.4,... % fly to screen lower right corner
%             'fly2scrnUL', 6.6,... % fly to screen upper left corner
%             'fly2scrnUR', 6.1,... % fly to screen upper right corner
%             'fly2horizL', 6.5,... % fly horizontally to vertical plane of screen, left side
%             'fly2horizR', 6.0,... % fly horizontally to vertical plane of screen, right side
%             'fly2topL', 0.7,... % horizontal plane of fly to screen top, vertically, on left
%             'fly2topR', 0.7,... % horizontal plane of fly to screen top, vertically, on right
%             'scrnX', 9,... % screen side length, x is vertical
%             'scrnY', 9); % screen side length, y is horizonal
    end 
    
    % methods implemented in subclass
    methods (Abstract)
        % Displays the stimulus using Psychtoolbox. Runs stimulus for
        %  duration specified by user OR until user presses the keyboard.
        % Inputs:
        %   obj - this object, always first input
        %   window - psychtoolbox pointer to window visual stimuli are
        %       presented to
        %   ifi - visual stimulus presentation interframe interval
        %   stimDuration - how long to present stimulus for
        displayStim(obj, window, ifi, stimDuration) 
        
        % Reconstructs the stimulus presented, for analysis
        % Inputs:
        %   obj - this object, always first input
        %   stimEpochStartTimes - lightStartTimes and darkStartTimes,
        %       sorted
        %   lightStartTimes - times when photodiode changed from dark to
        %       light
        %   darkStartTimes - times when photodiode changed from light to
        %       dark
        %   order - order of entries in stimEpochStartTimes
        % Outputs:
        %   rcStim - stimulus value at each stimEpochStartTime
        %   rcStimInd - epoch index at each stimEpochStartTime
        [rcStim, rcStimInd] = reconstructStim(obj, stimEpochStartTimes, ...
            lightStartTimes, darkStartTimes, order)
        
        % selects imaging frames to be used as F0 in computing dF/F during
        %  analysis
        % Inputs:
        %   obj- this object, always first input
        %   imgFrameTimes - start times of each imaging frame
        %   bksSignal - background-subtracted fluorescence trace
        %   lightStartTimes - times when photodiode changed from dark to
        %       light
        %   darkStartTimes - times when photodiode changed from light to
        %       dark
        % Outputs:
        %   baselineFrameTimes - times to use to calculate dF/F baseline
        %   baselineSignals - background-subtracted fluorescence values
        %       corresponding to baselineFrameTimes
        [baselineFrameTimes, baselineSignals] = selectDFFbaseline(obj,...
            imgFrameTimes, bksSignal, lightStartTimes, darkStartTimes)
    end 
    
    % methods implemented here but are called from the subclass object 
    methods
        % Constructor
        function obj = Stimulus(txtFile)
            obj.TxtFileName = txtFile;
        end 
        
        % Exports data as a .mat file in folder specified by user. Assumes
        % current directory is the correct fly folder (e.g. 160825_fly1)
        function exportData(obj, homeDir, stimCodeDir)
            fprintf('Navigate to and select the fly folder for saving your data in. \n');
            saveDir = uigetdir(homeDir, '');
            cd(saveDir);
%             saveFolder = input('Enter a name for the time series folder to save your stimulus data: ', 's');
            saveFolder = obj.Out.imDataName; % time series folder name
            
            % use default folder name if user doesn't specify one or writes
            % a folder name that already exists (avoids overwriting)
            if  isempty(saveFolder) || isdir(saveFolder)
                saveFolder = datestr(now, 'yymmdd_HH_MM_SS');
                fprintf('Invalid folder name. Default folder created: %s \n', saveFolder);
            end 
         
            mkdir(saveFolder);
            cd(saveFolder);
            save stim obj % save the stimulus object 
            cd(stimCodeDir);
            display('Saved!');
        end 
        
        % Saves the photodiode output and other timing info (specify!)
        function saveData(obj, photodiodeData, photodiodeTime, ...
                imagingFrameTime)
            obj.Out.date = datestr(now);
            obj.Out.pdData = photodiodeData; % photodiode voltage reading
            obj.Out.pdTime = photodiodeTime; 
            obj.Out.imFrameTime = imagingFrameTime; % imaging frame times
%             obj.Out.zdepth = input('Enter the z-depth you used: ');
            obj.Out.imDataName = input(...
                'Enter the name of this time series in the .lif file: ', 's');
        end 
        
        % Reads the .txt file and sets the parameters of this stimulus
        % class to the corresponding values in the .txt file.
        % How this works: Because the subclass object is also an object 
        % of the superclass, this superclass method recognizes that it is 
        % operating on the subclass object
        function setParams(obj)
            fileID = fopen(obj.TxtFileName); % open the txt file
            line1 = textscan(fileID, '%s %f', 1); % number of params
            line2 = textscan(fileID, '%s %f', 1); % number of epochs
            line3 = textscan(fileID, '%s %s', 1); % stimulus function name
           
            % save the values
            obj.nParams = line1{2};  
            obj.nEpochs = line2{2};   
%             className = line3{2};
            
%             % check to ensure stimulus class in text file matches this one
%             if ~strcmp(className, 'FullFieldFlash')
%                 fprintf('Error! Class name mismatch');
%                 return;
%             end 
            
            % loop through txt file and list of parameters to assign values 
            % to the expected parameters 
            while ~feof(fileID)
                pstr = textscan(fileID,'%s', 1); % name of parameter
                paramName = char(pstr{1}); % access the cell to save the name 
                formatStr = repmat('%f ', 1, obj.nEpochs); % assumes user set "EPOCHS" correctly
                values = textscan(fileID, formatStr, 1); % parameter values
                
                % loop through parameter list
                for i = 1:length(obj.ParamList) 
                    % if parameter name in txt file matches name in class's
                    % parameter list, assign the values of the paramters to
                    % the corresponding property in this class
                    if strcmp(paramName, obj.ParamList{i})
                        % get rid of empty cell elements and NaNs
                        pcell = values(1:end);
                        pcell = pcell(~cellfun('isempty', pcell));
                        pmat = cell2mat(pcell);
                        pcell = pcell(~isnan(pmat)); 
                        % save the non-Nan param values
                        pname = matlab.lang.makeValidName(paramName); 
                        obj.(pname) = pcell;
                    end 
                end 
            end 
            fclose(fileID);
        end 
        
        % Saves the NIDAQ scan rate 
        function setNIDAQScanRate(obj, rate)
            obj.NIDAQScanRate = rate;
        end 
    end 
end 
