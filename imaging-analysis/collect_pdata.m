% collect_pdata.m
%
% Function that asks for a fly folder (that contains folders for each
%  individual time series and copies pData files into pData_dir. Appends
%  name of fly folder to the pData file name.
%
% No input or output, but asks user for fly folder through GUI and saves
%  pData file into pData_dir as side effect.
%
function collect_pdata

[fname]=uigetdir;
d=dir(fname);
pData_dir = '/Users/hyang/Documents/New Imaging Analysis/pData/';


num=0;
for i=3:length(d)
    if d(i).isdir
        cd([fname filesep d(i).name]);
        ind = strfind(fname,filesep);
        flyID = fname(ind(end)+1:end);
        e = dir('*_pData.mat');
        if(~isempty(e))
            copyfile(e.name,[pData_dir flyID '_' e.name]);
        end
        cd ..
    end
end


