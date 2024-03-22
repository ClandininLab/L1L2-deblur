% FullFieldFlash is a subclass of the Stimulus superclass.
% last update: 08.19.16

classdef FullFieldFlash < handle
    properties 
        TxtFileName
        nParams
        nEpochs
        ParamList = {'Duration', 'Contrast'};
        Duration
        Contrast
        Out
    end 
    methods
        % Constructor
        function obj = FullFieldFlash(txtFile)
            obj.TxtFileName = txtFile;
        end 
        
        % Displays the stimulus using the parameters specified in the
        % .txt file
        function displayStim(obj, window, ifi)
            % Contrast values for each epoch. Currently either 1 or 0,
            % we'll make this more rigorous in coming versions
            darkContrast = obj.Contrast{1};
            lightContrast = obj.Contrast{2};
            
            % Seconds per stimulus epoch, assuming that light and dark
            % durations are of equal length
            epochDur = obj.Duration{1};
            
            % --- Psychtoolbox code --- *
            startTime = Screen('Flip', window);
            
            frame = 0;
            % keep playing the stimulus until the user presses a key on the keyboard
            while ~KbCheck 
                currEpoch = floor((frame*ifi)/epochDur);
                % light epoch
                if (mod(currEpoch, 2))
                    Screen('FillRect', window, lightContrast); % draw the stimulus
                    Screen('Flip', window, startTime + frame*ifi); % flip window at the frame rate
                % dark epoch
                else
                    Screen('FillRect', window, darkContrast);
                    Screen('Flip', window, startTime + frame*ifi);  
                end 
                pause(0) % super important - without this NIDAQ won't scan
                frame = frame + 1;
            end 
            
            sca;
        end 
        
        % Saves the data and timing outputs of the photodiode into this
        % class's Out struct
        function saveData(obj, photodiodeData, photodiodeTime)
            obj.Out.pdData = photodiodeData;
            obj.Out.pdTime = photodiodeTime;
        end 
        
        % Reads the .txt file and sets the parameters of this stimulus
        % class to the corresponding values in the .txt file
        function obj = setParams(obj)
            fileID = fopen(obj.TxtFileName); % open the txt file
            line1 = textscan(fileID, '%s %f', 1); % read line 1 to get number of params
            line2 = textscan(fileID, '%s %f', 1); % read line 2 to get number of epochs
            line3 = textscan(fileID, '%s %s', 1); % read line 3 to get stimulus function name
           
            % save the values
            obj.nParams = line1{2}; 
            obj.nEpochs = line2{2};   
            className = line3{2};
            
            % check to ensure stimulus class in text file matches this one
            if ~strcmp(className, 'FullFieldFlash')
                fprintf('Error! Class name mismatch');
                return;
            end 
            
            % loop through txt file and list of parameters to assine values 
            % to the expected parameters 
            while ~feof(fileID)
                pstr = textscan(fileID,'%s', 1); % name of paramter
                paramName = char(pstr{1}); % access the cell to save the name 
                values = textscan(fileID,'%f %f %f %f %f %f %f', 1); % parameter values
                
                % loop through parameter list
                for i = 1:length(obj.ParamList) 
                    % if parameter name in txt file matches name in class's
                    % parameter list, assign the values of the paramters to
                    % the corresponding property in this class
                    if strcmp(paramName, obj.ParamList{i})
                        pname = matlab.lang.makeValidName(paramName); 
                        obj.(pname) = cell(1:obj.nEpochs);  % assumes user set "EPOCHS" correctly 
                        obj.(pname) = values(1:obj.nEpochs); % assign values to parameter property 
                    end 
                end 
            end 
        end 
    end 
end 