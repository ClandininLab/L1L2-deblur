% returnMetricSAF.m
%
% Function that takes in metricStrctArray (from
%  compute_bootstrappedMetrics), specified field, and contrast (light/dark)
%  and returns specified metric value for every element in array
%
% INPUT
%   metricStrctArry - saved metric struct array from
%       compute_bootstrappedMetrics.m
%   field - field of array to return
%   whichContrast - dark = 1; light = 2
%
% OUPUT
%   returnMetric - vector (or matrix) of returned metric nElements x size
%       of field; i.e. column vector for single values; matrix with each
%       row corresponding to different element for vector values
%
% Updated: HHY 10/12/17
%

function returnMetric = returnMetricSAF(metricStrctArry,...
    field, whichContrast)
    if (isfield(metricStrctArry,field)) % only if field is one specified
        % preallocate
        returnMetric = zeros(length(metricStrctArry),...
            size(metricStrctArry(1).(field),2));
        % get value for each element
        for i=1:length(metricStrctArry)
            if (size(metricStrctArry(1).(field),2) == 1)
                returnMetric(i) = ...
                    metricStrctArry(i).(field)(whichContrast);
            else
                returnMetric(i,:) = ...
                    metricStrctArry(i).(field)(whichContrast,:);
            end
        end
    else
        % return empty vector if field is not a valid
        returnMetric = []; 
    end
end