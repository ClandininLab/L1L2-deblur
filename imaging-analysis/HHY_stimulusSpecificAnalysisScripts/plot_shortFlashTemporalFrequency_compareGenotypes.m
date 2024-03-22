% plot_shortFlashTemporalFrequency_compareGenotypes.m
%
% quick and dirty script to plot temporal frequency power spectra from 
% multiple genotypes on top of each other
%
% for L2 project
%
% 9/25/17

%% L2 ort rescue

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 60];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, ort1, Df(3R)BSC809, UAS-ort
saveName1 = 'L2_ASAP2f_UASort_ort1_Df3RBSC809';
load(saveName1, 'roiDataMat', 'iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, ort[1]/Df(3R)BSC809, UAS-ort experimental';

% L2>ASAP2f, Df(3R)BSC809, UAS-ort
saveName2 = 'L2_ASAP2f_UASort_Df3RBSC809';
load(saveName2, 'roiDataMat', 'iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, Df(3R)BSC809/+, UAS-ort control';

% L2>ASAP2f, ort1, UAS-ort
saveName3 = 'L2_ASAP2f_UASort_ort1';
load(saveName3, 'roiDataMat', 'iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f, ort[1]/+, UAS-ort control';

% L2>ASAP2f, ort1, Df(3R)BSC809
saveName4 = 'L2_ASAP2f_ort1_Df3RBSC809';
load(saveName4, 'roiDataMat');
respROIMat4 = roiDataMat(:,2);
titleStr4 = 'L2>ASAP2f, ort[1]/Df(3R)BSC809 control';

%% L2>>ASAP2f and TNT (GMR16H03-lexA) - silence L2, img L2
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 60];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f and TNT experimental
saveName1 = 'L2lexA_ASAP2f_lexAopTNT';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, TNT experimental';

% L2>ASAP2f control
saveName2 = 'L2lexA_ASAP2f';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f control';

%% L2>>ASAP2f and TNT (21Dhh>>ASAP2f; GMR16H03-lexA>>lexAop-TNT) 
% silence L2 image L2
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 60];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f and TNT experimental
saveName1 = '21Dhh_ASAP2f_GMR16H03_TNT';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, TNT experimental';

% L2>ASAP2f GMR16H03-lexA control
saveName2 = '21Dhh_ASAP2f_GMR16H03';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f GMR16H03-lexA control';

% L2>ASAP2f lexAop-TNT control
saveName3 = '21Dhh_ASAP2f_TNT_200pwm';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f lexAop-TNT control';

%% L2>>ASAP2f, 200 vs 20 pwm, no CDM
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 60];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, 200 pwm
saveName1 = '21Dhh_ASAP2f_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, 200pwm';

% L2>ASAP2f, 20 pwm
saveName2 = '21Dhh_ASAP2f_20pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, 20pwm';

%% L2>>ASAP2f, +/- CDM, 200 pwm
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 60];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, 200 pwm
saveName1 = '21Dhh_ASAP2f_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, 200pwm';

% L2>ASAP2f, 200 pwm + CDM
saveName2 = '21Dhh_ASAP2f_CDM_200pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, 200pwm, CDM';

%% L2>>ASAP2f, 200 vs 20 pwm, with CDM
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, 200 pwm, + CDM
saveName1 = '21Dhh_ASAP2f_CDM_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, 200pwm, CDM';

% L2>ASAP2f, 20 pwm, + CDM
saveName2 = '21Dhh_ASAP2f_CDM_20pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, 20pwm, CDM';

%% L2>>ASAP2f, +/- CDM, 20 pwm
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, 20 pwm
saveName1 = '21Dhh_ASAP2f_20pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, 20pwm';

% L2>ASAP2f, 20 pwm + CDM
saveName2 = '21Dhh_ASAP2f_CDM_20pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, 20pwm, CDM';

%% L2>>ASAP2f, 200 vs 20 pwm, +/- CDM

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, 200 pwm
saveName1 = '21Dhh_ASAP2f_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, 200pwm';

% L2>ASAP2f, 200 pwm, + CDM
saveName2 = '21Dhh_ASAP2f_CDM_200pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, 200pwm, CDM';

% L2>ASAP2f, 20 pwm
saveName3 = '21Dhh_ASAP2f_20pwm';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f, 20pwm';

% L2>ASAP2f, 20 pwm, + CDM
saveName4 = '21Dhh_ASAP2f_CDM_20pwm';
load(saveName4, 'roiDataMat','iResp');
respROIMat4 = roiDataMat(iResp,2);
titleStr4 = 'L2>ASAP2f, 20pwm, CDM';

%% L1>>ASAP2f, 200 vs 20 pwm, no CDM
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L1>ASAP2f, 200 pwm
saveName1 = 'GMR37E04_M1_ASAP2f_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L1>ASAP2f, 200pwm';

% L1>ASAP2f, 20 pwm
saveName2 = 'GMR37E04_M1_ASAP2f_20pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L1>ASAP2f, 20pwm';

%% L1>>ASAP2f, +/- CDM, 200 pwm

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L1>ASAP2f, 200 pwm
saveName1 = 'GMR37E04_M1_ASAP2f_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L1>ASAP2f, 200pwm';

% L1>ASAP2f, 200 pwm, + CDM
saveName2 = 'GMR37E04_M1_ASAP2f_CDM_200pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L1>ASAP2f, 200pwm, CDM';

%% L1>>ASAP2f, 200 vs 20 pwm, with CDM
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L1>ASAP2f, 200 pwm, + CDM
saveName1 = 'GMR37E04_M1_ASAP2f_CDM_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L1>ASAP2f, 200pwm, CDM';

% L1>ASAP2f, 20 pwm, + CDM
saveName2 = 'GMR37E04_M1_ASAP2f_CDM_20pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L1>ASAP2f, 20pwm, CDM';

%% L1>>ASAP2f, +/- CDM, 20 pwm

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L1>ASAP2f, 20 pwm
saveName1 = 'GMR37E04_M1_ASAP2f_20pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L1>ASAP2f, 20pwm';

% L1>ASAP2f, 20 pwm, + CDM
saveName2 = 'GMR37E04_M1_ASAP2f_CDM_20pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L1>ASAP2f, 20pwm, CDM';

%% L1>>ASAP2f, 200 vs 20 pwm, +/- CDM

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L1>ASAP2f, 200 pwm
saveName1 = 'GMR37E04_M1_ASAP2f_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L1>ASAP2f, 200pwm';

% L1>ASAP2f, 200 pwm, + CDM
saveName2 = 'GMR37E04_M1_ASAP2f_CDM_200pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L1>ASAP2f, 200pwm, CDM';

% L1>ASAP2f, 20 pwm
saveName3 = 'GMR37E04_M1_ASAP2f_20pwm';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L1>ASAP2f, 20pwm';

% L1>ASAP2f, 20 pwm, + CDM
saveName4 = 'GMR37E04_M1_ASAP2f_CDM_20pwm';
load(saveName4, 'roiDataMat','iResp');
respROIMat4 = roiDataMat(iResp,2);
titleStr4 = 'L1>ASAP2f, 20pwm, CDM';

%% L2>>ASAP2f, ort silencing, 200 pwm

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, ort[C1-3]-lexA>lexAop-TNT, 200 pwm
saveName1 = '21Dhh_ASAP2f_ortC1-3_TNT_200pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'ort[C1-3]-lexA>lexAop-TNT experimental, 200 pwm';

% L2>ASAP2f, ort[C1-3]-lexA control, 200 pwm
saveName2 = '21Dhh_ASAP2f_ortC1-3_200pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'ort[C1-3]-lexA control, 200 pwm';

% L2>ASAP2f, lexAop-TNT control, 200 pwm
saveName3 = '21Dhh_ASAP2f_TNT_200pwm';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'lexAop-TNT control, 200 pwm';


%% L2>>ASAP2f, ort silencing, 20 pwm

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170925';
cd(dataPath)

xScale = [0, 30];
yScale = [10E-3, 10];
pairedEpochs = [1 2];

% L2>ASAP2f, ort[C1-3]-lexA>lexAop-TNT, 20 pwm
saveName1 = '21Dhh_ASAP2f_ortC1-3_TNT_20pwm';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'ort[C1-3]-lexA>lexAop-TNT experimental, 20 pwm';

% L2>ASAP2f, ort[C1-3]-lexA control, 20 pwm
saveName2 = '21Dhh_ASAP2f_ortC1-3_20pwm';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'ort[C1-3]-lexA control, 20 pwm';

% L2>ASAP2f, lexAop-TNT control, 20 pwm
saveName3 = '21Dhh_ASAP2f_TNT_20pwm';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'lexAop-TNT control, 20 pwm';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compare 4

cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        tfPs1 = zeros(length(respROIMat1), length(respROIMat1(1).tfPs{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            tfPs1(n, :) = respROIMat1(n).tfPs{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(tfPs1,1);
        e1 = std(tfPs1,[],1)./sqrt(size(tfPs1,1));
        h1 = semilogy_err_patch_v2(...
            respROIMat1(1).tempFreq{pairedEpochs(i,j)},...
            m1 ,...
            e1, ...
            cm(1,:),(cm(1,:)+1)/2);
        
        tfPs2 = zeros(length(respROIMat2), length(respROIMat2(1).tfPs{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            tfPs2(n, :) = respROIMat2(n).tfPs{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(tfPs2,1);
        e2 = std(tfPs2,[],1)./sqrt(size(tfPs2,1));
        h2 = semilogy_err_patch_v2(...
            respROIMat2(1).tempFreq{pairedEpochs(i,j)},...
            m2,...
            e2, ...
            cm(2,:),(cm(2,:)+1)/2);
        
        tfPs3 = zeros(length(respROIMat3), length(respROIMat3(1).tfPs{pairedEpochs(i,j)}));
        flyID3 = zeros(1, length(respROIMat3));
        for n = 1:length(respROIMat3)
            tfPs3(n, :) = respROIMat3(n).tfPs{pairedEpochs(i,j)};
            flyID3(n) = respROIMat3(n).flyID;
        end         
        m3 = mean(tfPs3,1);
        e3 = std(tfPs3,[],1)./sqrt(size(tfPs3,1));        
        h3 = semilogy_err_patch_v2(...
            respROIMat3(1).tempFreq{pairedEpochs(i,j)},...
            m3,...
            e3, ...
            cm(3,:),(cm(3,:)+1)/2);
        
        tfPs4 = zeros(length(respROIMat4), length(respROIMat4(1).tfPs{pairedEpochs(i,j)}));
        flyID4 = zeros(1, length(respROIMat4));
        for n = 1:length(respROIMat4)
            tfPs4(n, :) = respROIMat4(n).tfPs{pairedEpochs(i,j)};
            flyID4(n) = respROIMat4(n).flyID;
        end         
        m4 = mean(tfPs4,1);
        e4 = std(tfPs4,[],1)./sqrt(size(tfPs4,1));        
        h4 = semilogy_err_patch_v2(...
            respROIMat4(1).tempFreq{pairedEpochs(i,j)},...
            m4,...
            e4, ...
            cm(4,:),(cm(4,:)+1)/2);
       
        % plotting parameters
        xlabel('Temporal Frequency (Hz)');
        ylabel('Amplitude');
        xlim(xScale);
        ylim(yScale);

        
        legend([h1 h2 h3 h4], ...
            [titleStr1 ', N cells = ' num2str(size(tfPs1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(tfPs2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        [titleStr3 ', N cells = ' num2str(size(tfPs3,1)) ...
        ', N flies = ' num2str(length(unique(flyID3)))],...
        [titleStr4 ', N cells = ' num2str(size(tfPs4,1)) ...
        ', N flies = ' num2str(length(unique(flyID4)))],...        
        'Location','northeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% compare 3
cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        tfPs1 = zeros(length(respROIMat1), length(respROIMat1(1).tfPs{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            tfPs1(n, :) = respROIMat1(n).tfPs{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(tfPs1,1);
        e1 = std(tfPs1,[],1)./sqrt(size(tfPs1,1));
        h1 = semilogy_err_patch_v2(...
            respROIMat1(1).tempFreq{pairedEpochs(i,j)},...
            m1 ,...
            e1, ...
            cm(1,:),(cm(1,:)+1)/2);
        
        tfPs2 = zeros(length(respROIMat2), length(respROIMat2(1).tfPs{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            tfPs2(n, :) = respROIMat2(n).tfPs{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(tfPs2,1);
        e2 = std(tfPs2,[],1)./sqrt(size(tfPs2,1));
        h2 = semilogy_err_patch_v2(...
            respROIMat2(1).tempFreq{pairedEpochs(i,j)},...
            m2,...
            e2, ...
            cm(2,:),(cm(2,:)+1)/2);
        
        tfPs3 = zeros(length(respROIMat3), length(respROIMat3(1).tfPs{pairedEpochs(i,j)}));
        flyID3 = zeros(1, length(respROIMat3));
        for n = 1:length(respROIMat3)
            tfPs3(n, :) = respROIMat3(n).tfPs{pairedEpochs(i,j)};
            flyID3(n) = respROIMat3(n).flyID;
        end         
        m3 = mean(tfPs3,1);
        e3 = std(tfPs3,[],1)./sqrt(size(tfPs3,1));        
        h3 = semilogy_err_patch_v2(...
            respROIMat3(1).tempFreq{pairedEpochs(i,j)},...
            m3,...
            e3, ...
            cm(3,:),(cm(3,:)+1)/2);
       
        % plotting parameters
        xlabel('Temporal Frequency (Hz)');
        ylabel('Amplitude');
        xlim(xScale);
        ylim(yScale);

        
        legend([h1 h2 h3], ...
            [titleStr1 ', N cells = ' num2str(size(tfPs1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(tfPs2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        [titleStr3 ', N cells = ' num2str(size(tfPs3,1)) ...
        ', N flies = ' num2str(length(unique(flyID3)))],...       
        'Location','northeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% compare 2
cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        tfPs1 = zeros(length(respROIMat1), length(respROIMat1(1).tfPs{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            tfPs1(n, :) = respROIMat1(n).tfPs{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(tfPs1,1);
        e1 = std(tfPs1,[],1)./sqrt(size(tfPs1,1));
        h1 = semilogy_err_patch_v2(...
            respROIMat1(1).tempFreq{pairedEpochs(i,j)},...
            m1 ,...
            e1, ...
            cm(1,:),(cm(1,:)+1)/2);
        
        tfPs2 = zeros(length(respROIMat2), length(respROIMat2(1).tfPs{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            tfPs2(n, :) = respROIMat2(n).tfPs{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(tfPs2,1);
        e2 = std(tfPs2,[],1)./sqrt(size(tfPs2,1));
        h2 = semilogy_err_patch_v2(...
            respROIMat2(1).tempFreq{pairedEpochs(i,j)},...
            m2,...
            e2, ...
            cm(2,:),(cm(2,:)+1)/2);
       
        % plotting parameters
        xlabel('Temporal Frequency (Hz)');
        ylabel('Amplitude');
        xlim(xScale);
        ylim(yScale);

        
        legend([h1 h2], ...
            [titleStr1 ', N cells = ' num2str(size(tfPs1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(tfPs2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...       
        'Location','northeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% compare 4

cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        tfPs1 = zeros(length(respROIMat1), length(respROIMat1(1).tfPs{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            tfPs1(n, :) = respROIMat1(n).tfPs{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(tfPs1,1);
        e1 = std(tfPs1,[],1)./sqrt(size(tfPs1,1));
        h1 = loglog_err_patch_v2(...
            respROIMat1(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m1(2:end),...
            e1(2:end), ...
            cm(1,:),(cm(1,:)+1)/2);
        
        tfPs2 = zeros(length(respROIMat2), length(respROIMat2(1).tfPs{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            tfPs2(n, :) = respROIMat2(n).tfPs{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(tfPs2,1);
        e2 = std(tfPs2,[],1)./sqrt(size(tfPs2,1));
        h2 = loglog_err_patch_v2(...
            respROIMat2(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m2(2:end),...
            e2(2:end), ...
            cm(2,:),(cm(2,:)+1)/2);
        
        tfPs3 = zeros(length(respROIMat3), length(respROIMat3(1).tfPs{pairedEpochs(i,j)}));
        flyID3 = zeros(1, length(respROIMat3));
        for n = 1:length(respROIMat3)
            tfPs3(n, :) = respROIMat3(n).tfPs{pairedEpochs(i,j)};
            flyID3(n) = respROIMat3(n).flyID;
        end         
        m3 = mean(tfPs3,1);
        e3 = std(tfPs3,[],1)./sqrt(size(tfPs3,1));        
        h3 = loglog_err_patch_v2(...
            respROIMat3(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m3(2:end),...
            e3(2:end), ...
            cm(3,:),(cm(3,:)+1)/2);
        
        tfPs4 = zeros(length(respROIMat4), length(respROIMat4(1).tfPs{pairedEpochs(i,j)}));
        flyID4 = zeros(1, length(respROIMat4));
        for n = 1:length(respROIMat4)
            tfPs4(n, :) = respROIMat4(n).tfPs{pairedEpochs(i,j)};
            flyID4(n) = respROIMat4(n).flyID;
        end         
        m4 = mean(tfPs4,1);
        e4 = std(tfPs4,[],1)./sqrt(size(tfPs4,1));        
        h4 = loglog_err_patch_v2(...
            respROIMat4(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m4(2:end),...
            e4(2:end), ...
            cm(4,:),(cm(4,:)+1)/2);
       
        % plotting parameters
        xlabel('Temporal Frequency (Hz)');
        ylabel('Amplitude');
        xlim(xScale);
        ylim(yScale);

        
        legend([h1 h2 h3 h4], ...
            [titleStr1 ', N cells = ' num2str(size(tfPs1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(tfPs2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        [titleStr3 ', N cells = ' num2str(size(tfPs3,1)) ...
        ', N flies = ' num2str(length(unique(flyID3)))],...
        [titleStr4 ', N cells = ' num2str(size(tfPs4,1)) ...
        ', N flies = ' num2str(length(unique(flyID4)))],...        
        'Location','northeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% compare 3 - loglog
cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        tfPs1 = zeros(length(respROIMat1), length(respROIMat1(1).tfPs{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            tfPs1(n, :) = respROIMat1(n).tfPs{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(tfPs1,1);
        e1 = std(tfPs1,[],1)./sqrt(size(tfPs1,1));
        h1 = loglog_err_patch_v2(...
            respROIMat1(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m1(2:end),...
            e1(2:end), ...
            cm(1,:),(cm(1,:)+1)/2);
        
        tfPs2 = zeros(length(respROIMat2), length(respROIMat2(1).tfPs{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            tfPs2(n, :) = respROIMat2(n).tfPs{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(tfPs2,1);
        e2 = std(tfPs2,[],1)./sqrt(size(tfPs2,1));
        h2 = loglog_err_patch_v2(...
            respROIMat2(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m2(2:end),...
            e2(2:end), ...
            cm(2,:),(cm(2,:)+1)/2);
        
        tfPs3 = zeros(length(respROIMat3), length(respROIMat3(1).tfPs{pairedEpochs(i,j)}));
        flyID3 = zeros(1, length(respROIMat3));
        for n = 1:length(respROIMat3)
            tfPs3(n, :) = respROIMat3(n).tfPs{pairedEpochs(i,j)};
            flyID3(n) = respROIMat3(n).flyID;
        end         
        m3 = mean(tfPs3,1);
        e3 = std(tfPs3,[],1)./sqrt(size(tfPs3,1));        
        h3 = loglog_err_patch_v2(...
            respROIMat3(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m3(2:end),...
            e3(2:end), ...
            cm(3,:),(cm(3,:)+1)/2);
       
        % plotting parameters
        xlabel('Temporal Frequency (Hz)');
        ylabel('Amplitude');
        xlim(xScale);
        ylim(yScale);

        
        legend([h1 h2 h3], ...
            [titleStr1 ', N cells = ' num2str(size(tfPs1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(tfPs2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        [titleStr3 ', N cells = ' num2str(size(tfPs3,1)) ...
        ', N flies = ' num2str(length(unique(flyID3)))],...       
        'Location','northeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end


%% compare 2 - loglog
cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        tfPs1 = zeros(length(respROIMat1), length(respROIMat1(1).tfPs{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            tfPs1(n, :) = respROIMat1(n).tfPs{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(tfPs1,1);
        e1 = std(tfPs1,[],1)./sqrt(size(tfPs1,1));
        h1 = loglog_err_patch_v2(...
            respROIMat1(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m1(2:end) ,...
            e1(2:end), ...
            cm(1,:),(cm(1,:)+1)/2);
        
        tfPs2 = zeros(length(respROIMat2), length(respROIMat2(1).tfPs{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            tfPs2(n, :) = respROIMat2(n).tfPs{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(tfPs2,1);
        e2 = std(tfPs2,[],1)./sqrt(size(tfPs2,1));
        h2 = loglog_err_patch_v2(...
            respROIMat2(1).tempFreq{pairedEpochs(i,j)}(2:end),...
            m2(2:end),...
            e2(2:end), ...
            cm(2,:),(cm(2,:)+1)/2);
       
        % plotting parameters
        xlabel('Temporal Frequency (Hz)');
        ylabel('Amplitude');
        xlim(xScale);
        ylim(yScale);

        
        legend([h1 h2], ...
            [titleStr1 ', N cells = ' num2str(size(tfPs1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(tfPs2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...       
        'Location','northeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% light vs dark 1 genotype
figure;
for j = 1:size(pairedEpochs, 2)
    tfPs1 = zeros(length(respROIMat1), length(respROIMat1(1).tfPs{pairedEpochs(i,j)}));
    flyID1 = zeros(1, length(respROIMat1));
    for n = 1:length(respROIMat1)
        tfPs1(n, :) = respROIMat1(n).tfPs{pairedEpochs(i,j)};
        flyID1(n) = respROIMat1(n).flyID;
    end 
    m1 = mean(tfPs1,1);
    e1 = std(tfPs1,[],1)./sqrt(size(tfPs1,1));
    h(j) = loglog_err_patch_v2(...
        respROIMat1(1).tempFreq{pairedEpochs(i,j)}(2:end),...
        m1(2:end) ,...
        e1(2:end), ...
        cm(j,:),(cm(j,:)+1)/2);
end
% plotting parameters
xlabel('Temporal Frequency (Hz)');
ylabel('Amplitude');
xlim(xScale);
ylim(yScale);
        
legend(h,'Dark','Light','Location','northeast');