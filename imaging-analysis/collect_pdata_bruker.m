% collect_pdata_bruker.m
%
% Function that asks for a fly folder (that contains folders for each
%  individual time series and copies pData files into pData_dir. Appends
%  name of fly folder to the pData file name.
%
% No input or output, but asks user for fly folder through GUI and saves
%  pData file into pData_dir as side effect.
%
function collect_pdata_bruker

[fname]=uigetdir;
d=dir(fname);
pData_dir = '/Volumes/LEGTRACKING/L1L2_Michelle/220220/pdata_220328';


for i=3:length(d)
    if d(i).isdir
        cd([fname filesep d(i).name]);
        ind = strfind(fname,filesep);
        e = dir('TSeries*_pData.mat');
        if(~isempty(e))
            copyfile(e.name,[pData_dir filesep e.name]);
        end
        cd ..
    end
end

end


