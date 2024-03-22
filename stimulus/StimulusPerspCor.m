% StimulusPerspCor is an abstract class that is a subclass of the abstract
%  class Stimulus.
% Adds to Stimulus methods and properties for perspective correction.
%
% Abstract Properties
%   Init - output variables of stimulus initialization
%
% Properties
%   txtrParams - struct of texture parameters for perspective correction
%
% Abstract Methods
%   initStim - function to initalize stimulus
%
% Last updated: 1/15/17 HHY
%

classdef (Abstract) StimulusPerspCor < Stimulus
    % properties implemented in subclass
    properties (Abstract)
        Init
    end
    % properties implemented here
    properties 
        % strct of texture parameters for perspective correction
        txtrParams = struct(...
            'txtrSize',90,... % size of texture, in degrees of fly's visual field
            'txtrRes', 2,... % texture resolution, in pixels per degree
            'xCorner', 65,... % upper right corner of screen in texture, x
            'yCorner', 85); % upper right corner of screen in texture, y
    end
    
    % methods implemented in subclass
    methods (Abstract)
        % function to initalize stimulus
        initStim(obj,window,ifi)
    end
    
    % methods implemented here
    methods
        % constructor; inherited from superclass but must be redefined here
        function obj = StimulusPerspCor(txtFile)
            obj@Stimulus(txtFile);
        end
    end
end