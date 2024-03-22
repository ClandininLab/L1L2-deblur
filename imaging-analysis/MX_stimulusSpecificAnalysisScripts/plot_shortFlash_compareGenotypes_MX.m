% plot_shortFlash_compareGenotypes.m
%
% quick and dirty script to plot short flash responses from multiple
% genotypes on top of each other
%
% for L2 cell type specific silencing
%
% 5/30/17

clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170518';
cd(dataPath)

yScale = [-0.04 0.08];
pairedEpochs = [1 2];
inv = 1;

%% L2 ort rescue

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

%% L2 ort rescue, L2 silencing

% L2>ASAP2f, ort1, Df(3R)BSC809, UAS-ort, UAS-shi[ts] experimental
saveName1 = 'L2_ASAP2f_UASort_UASshits_ort1_Df3RBSC809_1h37';
load(saveName1, 'roiDataMat', 'iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, ort[1]/Df(3R)BSC809, UAS-ort, UAS-shit[ts] experimental';

% L2>ASAP2f, ort1, Df(3R)BSC809, UAS-ort
saveName2 = 'L2_ASAP2f_UASort_ort1_Df3RBSC809_1h37';
load(saveName2, 'roiDataMat', 'iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, ort[1]/Df(3R)BSC809, UAS-ort control';

% L2>ASAP2f, ort1, UAS-ort, UAS-shi[ts]
saveName3 = 'L2_ASAP2f_UASort_UASshits_ort1_1h37';
load(saveName3, 'roiDataMat', 'iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f, ort[1]/+, UAS-ort, UAS-shi[ts] control';

% L2>ASAP2f, Df(3R)BSC809, UAS-ort, UAS-shi[ts]
saveName4 = 'L2_ASAP2f_UASort_UASshits_Df3RBSC809_1h37';
load(saveName4, 'roiDataMat', 'iResp');
respROIMat4 = roiDataMat(iResp,2);
titleStr4 = 'L2>ASAP2f, Df(3R)BSC809/+, UAS-ort, UAS-shi[ts] control';


%% L2 ort TNT

% L2>ASAP2f, ort>TNT
saveName1 = 'L2_ASAP2f_ort-TNT';
load(saveName1, 'respROIMat');
respROIMat1 = respROIMat(:,2);
titleStr1 = 'L2>ASAP2f, ort>TNT';

% L2>ASAP2f, ort ctrl
saveName2 = 'L2_ASAP2f_ort';
load(saveName2, 'respROIMat');
respROIMat2 = respROIMat(:,2);
titleStr2 = 'L2>ASAP2f, ort control';

% L2>ASAP2f, TNT control
saveName3 = 'L2_ASAP2f_TNT';
load(saveName3, 'respROIMat');
respROIMat3 = respROIMat(:,2);
titleStr3 = 'L2>ASAP2f, TNT control';

%% L2 ort TNT +/- TTX
% L2>ASAP2f, ort>TNT no TTX
saveName1 = 'L2_ASAP2f_ort-TNT';
load(saveName1, 'respROIMat');
respROIMat1 = respROIMat(:,2);
titleStr1 = 'no TTX';

% L2>ASAP2f, ort>TNT + TTX
dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170227';
cd(dataPath)
saveName2 = 'L2_ASAP2f_ort-TNT_TTX';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = '+TTX';

titleStr = 'L2>ASAP2f, ort>TNT, +/- TTX';

%% for ort>shi[ts] w/L2-lexA
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170117';
cd(dataPath)

yScale = [-0.04 0.04];
pairedEpochs = [1 2];
inv = 1;

% L2>ASAP2f, ort>shi[ts]
saveName1 = 'L2lexA-ASAP2f_ort-shits';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, ort>shi[ts]';

% L2>ASAP2f, ort ctrl
saveName2 = 'L2lexA-ASAP2f_ort';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, ort control';

% L2>ASAP2f, ort>TNT
saveName3 = 'L2lexA-ASAP2f_shits';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f, shi[ts] control';

%% for Lawf2>TNT
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170518';
cd(dataPath)

yScale = [-0.04 0.06];
pairedEpochs = [1 2];
inv = 1;

% L2>ASAP2f, Lawf2>TNT
saveName1 = 'L2lexA_ASAP2f_Lawf2-TNT';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, Lawf2>TNT';

% L2>ASAP2f, Lawf2 Gal4 control
saveName2 = 'L2lexA_ASAP2f_Lawf2';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, Lawf2-Gal4 control';

% L2>ASAP2f, TNT control
saveName3 = 'L2lexA_ASAP2f_TNT';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f, TNT control';

%% for Lawf1>TNT
% clear all
% close all
% clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/170518';
cd(dataPath)

yScale = [-0.04 0.06];
pairedEpochs = [1 2];
inv = 1;

% L2>ASAP2f, Lawf1>TNT
saveName1 = 'L2lexA_ASAP2f_Lawf1-TNT';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, Lawf1>TNT';

% L2>ASAP2f, Lawf1 Gal4 control
saveName2 = 'L2lexA_ASAP2f_Lawf1';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f, Lawf1-Gal4 control';

% L2>ASAP2f, TNT control
saveName3 = 'L2lexA_ASAP2f_TNT';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f, TNT control';

%% for ort>TNT
clear all
close all

yScale = [-0.06 0.06];
pairedEpochs = [1 2];
inv = 1;

% L2>ASAP2f, ort>TNT
dataPath1 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\L2-ort-silencing\GMR16H03-ort-TNT';
cd(dataPath1);
saveName1 = 'L2_ASAP2f_ort_TNT_170527';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'L2>ASAP2f, ort>TNT';


% L2>ASAP2f ort Gal4, No-TNT control
dataPath2 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\L2-ort-silencing\GMR16H03_ort-Gal4_No-TNT-control';
cd(dataPath2);
saveName2 = 'L2-ASAP2f_ort-Gal4_170527';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'L2>ASAP2f ort-Gal4 (no-TNT control)';

% L2>ASAP2f, TNT control
dataPath3 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\L2-ort-silencing\GMR16H03_TNT_No-Gal4-Control';
cd(dataPath3);
saveName3 = 'L2_GMR16H03_ASAP2f_TNT_170419';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'L2>ASAP2f UAS-TNT (no-Gal4 control)';

%% 21Dhh>ASAP2f Different mean luminances 
clear all
close all

yScale = [-0.06 0.06];
pairedEpochs = [1 2];
inv = 1;

% 200 pwm
dataPath1 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\21Dhh-ASAP2f';
cd(dataPath1);
saveName1 = '21Dhh_ASAP2f';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = '21Dhh>ASAP2f 200pwm';

% 75 pwm
dataPath2 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\L2-meanLuminances\3flies';
cd(dataPath2);
saveName2 = 'L2_ASAP2f_meanLum';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,3);
titleStr2 = '21Dhh>ASAP2f 75pwm';
% 20 pwm
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = '21Dhh>ASAP2f 20pwm';

%% octopamine
clear all
close all

yScale = [-0.06 0.06];
pairedEpochs = [1 2];
inv = 1;

dataPath1 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\21Dhh-ASAP2f';
cd(dataPath1);

% L2>ASAP2f
saveName1 = '21Dhh_ASAP2f';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = '21Dhh>ASAP2f';

dataPath2 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\L2-octopamine\21Dhh-ASAP2f+octopamine';
cd(dataPath2);

% L2>ASAP2f + octopamine
saveName2 = '21Dhh_ASAP2f_octopamine';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = '21Dhh>ASAP2f octopamine';

%% 21Dhh>>ASAP2f + CDM 20uM
clear all
close all

yScale = [-0.06 0.06];
pairedEpochs = [1 2];
inv = 1;

dataPath1 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\21Dhh-ASAP2f';
cd(dataPath1);

% L2>ASAP2f
saveName1 = '21Dhh_ASAP2f';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = '21Dhh>ASAP2f';

dataPath2 = 'D:\Marjorie\ClandininLabStanford\Imaging\AnalyzedData\L2-octopamine\21Dhh-ASAP2f+CDM20uM';
cd(dataPath2);

% L2>ASAP2f + CDM 20uM
saveName2 = '21Dhh_ASAP2f_CDM20uM';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = '21Dhh>ASAP2f + CDM 20uM';

% % L2>ASAP2f + CDM 20uM wait 15+min
% saveName2 = '21Dhh_ASAP2f_CDM20uM_wait15+min';
% load(saveName2, 'roiDataMat','iResp');
% respROIMat2 = roiDataMat(iResp,2);
% titleStr2 = '21Dhh>ASAP2f + CDM 20uM wait 15+min';

%% compare 4

flashDurs = cell2mat(respROIMat1(1).stimDat.FlashDuration);
grayDur = respROIMat1(1).stimDat.GrayDuration{1}; 
seqDur = flashDurs + grayDur;

cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    iDur = pairedEpochs(i);
    seqDur = flashDurs(iDur) + grayDur;
    flashDur = flashDurs(iDur);

    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        hold on;
        rats1 = zeros(length(respROIMat1), length(respROIMat1(1).rats{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            rats1(n, :) = respROIMat1(n).rats{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(rats1,1);
        e1 = std(rats1,[],1)./sqrt(size(rats1,1));
        h1 = plot_err_patch_v2(...
            respROIMat1(1).t{pairedEpochs(i,j)},...
            m1 ,...
            e1, ...
            cm(1,:),(cm(1,:)+1)/2);
        
        rats2 = zeros(length(respROIMat2), length(respROIMat2(1).rats{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            rats2(n, :) = respROIMat2(n).rats{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(rats2,1);
        e2 = std(rats2,[],1)./sqrt(size(rats2,1));
        h2 = plot_err_patch_v2(...
            respROIMat2(1).t{pairedEpochs(i,j)},...
            m2,...
            e2, ...
            cm(2,:),(cm(2,:)+1)/2);
        
        rats3 = zeros(length(respROIMat3), length(respROIMat3(1).rats{pairedEpochs(i,j)}));
        flyID3 = zeros(1, length(respROIMat3));
        for n = 1:length(respROIMat3)
            rats3(n, :) = respROIMat3(n).rats{pairedEpochs(i,j)};
            flyID3(n) = respROIMat3(n).flyID;
        end         
        m3 = mean(rats3,1);
        e3 = std(rats3,[],1)./sqrt(size(rats3,1));        
        h3 = plot_err_patch_v2(...
            respROIMat3(1).t{pairedEpochs(i,j)},...
            m3,...
            e3, ...
            cm(3,:),(cm(3,:)+1)/2);
        
        rats4 = zeros(length(respROIMat4), length(respROIMat4(1).rats{pairedEpochs(i,j)}));
        flyID4 = zeros(1, length(respROIMat4));
        for n = 1:length(respROIMat4)
            rats4(n, :) = respROIMat4(n).rats{pairedEpochs(i,j)};
            flyID4(n) = respROIMat4(n).flyID;
        end         
        m4 = mean(rats4,1);
        e4 = std(rats4,[],1)./sqrt(size(rats4,1));        
        h4 = plot_err_patch_v2(...
            respROIMat4(1).t{pairedEpochs(i,j)},...
            m4,...
            e4, ...
            cm(4,:),(cm(4,:)+1)/2);
       
        % plotting parameters
        xScale = [0, seqDur];
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        xlim(xScale);
        ylim(yScale);
        line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line
        
        if (inv)
            patchY = yScale(1);
        else
            patchY = yScale(2);
        end 
        
        if (inv)
            set(gca,'YDir','reverse');
        end
        
        if (j==2)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [1 1 1]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        elseif (j==1)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [0 0 0]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        end
        line([flashDurs(iDur) flashDurs(iDur)], yScale, 'color',[0 0 0],...
            'linestyle','--');
        legend([h1 h2 h3 h4], ...
            [titleStr1 ', N cells = ' num2str(size(rats1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(rats2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        [titleStr3 ', N cells = ' num2str(size(rats3,1)) ...
        ', N flies = ' num2str(length(unique(flyID3)))],...
        [titleStr4 ', N cells = ' num2str(size(rats4,1)) ...
        ', N flies = ' num2str(length(unique(flyID4)))],...        
        'Location','southeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% compare 3
flashDurs = cell2mat(respROIMat1(1).stimDat.FlashDuration);
grayDur = respROIMat1(1).stimDat.GrayDuration{1}; 
seqDur = flashDurs + grayDur;

cm = colormap('lines');
%close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    iDur = pairedEpochs(i);
    seqDur = flashDurs(iDur) + grayDur;
    flashDur = flashDurs(iDur);

    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        hold on;
        rats1 = zeros(length(respROIMat1), length(respROIMat1(1).rats{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            rats1(n, :) = respROIMat1(n).rats{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(rats1,1);
        e1 = std(rats1,[],1)./sqrt(size(rats1,1));
        h1 = plot_err_patch_v2(...
            (0:length(m1))* respROIMat1(1).BIN_SHIFT,...
            [m1 m1(1)],...
            [e1 e1(1)], ...
            cm(1,:),(cm(1,:)+1)/2);
        
        rats2 = zeros(length(respROIMat2), length(respROIMat2(1).rats{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            rats2(n, :) = respROIMat2(n).rats{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(rats2,1);
        e2 = std(rats2,[],1)./sqrt(size(rats2,1));
        h2 = plot_err_patch_v2(...
            (0:length(m2))* respROIMat2(1).BIN_SHIFT,...
            [m2 m2(1)],...
            [e2 e2(1)], ...
            cm(2,:),(cm(2,:)+1)/2);
        
        rats3 = zeros(length(respROIMat3), length(respROIMat3(1).rats{pairedEpochs(i,j)}));
        flyID3 = zeros(1, length(respROIMat3));
        for n = 1:length(respROIMat3)
            rats3(n, :) = respROIMat3(n).rats{pairedEpochs(i,j)};
            flyID3(n) = respROIMat3(n).flyID;
        end         
        m3 = mean(rats3,1);
        e3 = std(rats3,[],1)./sqrt(size(rats3,1));        
        h3 = plot_err_patch_v2(...
            (0:length(m3))* respROIMat3(1).BIN_SHIFT,...
            [m3 m3(1)],...
            [e3 e3(1)], ...
            cm(3,:),(cm(3,:)+1)/2);
       
        % plotting parameters
        xScale = [0, seqDur];
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        xlim(xScale);
        ylim(yScale);
        line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line
        
        if (inv)
            patchY = yScale(1);
        else
            patchY = yScale(2);
        end 
        
        if (inv)
            set(gca,'YDir','reverse');
        end
        
        if (j==2)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [1 1 1]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        elseif (j==1)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [0 0 0]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        end
        line([flashDurs(iDur) flashDurs(iDur)], yScale, 'color',[0 0 0],...
            'linestyle','--');
        legend([h1 h2 h3], ...
            [titleStr1 ', N cells = ' num2str(size(rats1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(rats2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        [titleStr3 ', N cells = ' num2str(size(rats3,1)) ...
        ', N flies = ' num2str(length(unique(flyID3)))],...
        'Location','southeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% for comparing 2
flashDurs = cell2mat(respROIMat1(1).stimDat.FlashDuration);
grayDur = respROIMat1(1).stimDat.GrayDuration{1}; 
seqDur = flashDurs + grayDur;

cm = colormap('lines');
close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    iDur = pairedEpochs(i);
    seqDur = flashDurs(iDur) + grayDur;
    flashDur = flashDurs(iDur);

    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        hold on;
        rats1 = zeros(length(respROIMat1), length(respROIMat1(1).rats{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            rats1(n, :) = respROIMat1(n).rats{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(rats1,1);
        e1 = std(rats1,[],1)./sqrt(size(rats1,1));
        h1 = plot_err_patch_v2(...
            (0:length(m1))* respROIMat1(n).BIN_SHIFT,...
            [m1 m1(1)],...
            [e1 e1(1)], ...
            cm(1,:),(cm(1,:)+1)/2);
        
        rats2 = zeros(length(respROIMat2), length(respROIMat2(1).rats{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            rats2(n, :) = respROIMat2(n).rats{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(rats2,1);
        e2 = std(rats2,[],1)./sqrt(size(rats2,1));
        h2 = plot_err_patch_v2(...
            (0:length(m2))* respROIMat2(n).BIN_SHIFT,...
            [m2 m2(1)],...
            [e2 e2(1)], ...
            cm(2,:),(cm(2,:)+1)/2);
       
        % plotting parameters
        xScale = [0, seqDur];
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        xlim(xScale);
        ylim(yScale);
        line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line
        
        if (inv)
            patchY = yScale(1);
        else
            patchY = yScale(2);
        end 
        
        if (inv)
            set(gca,'YDir','reverse');
        end
        
        if (j==2)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [1 1 1]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        elseif (j==1)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [0 0 0]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        end
        line([flashDurs(iDur) flashDurs(iDur)], yScale, 'color',[0 0 0],...
            'linestyle','--');
        legend([h1 h2], ...
            [titleStr1 ', N cells = ' num2str(size(rats1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(rats2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        'Location','southeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end
%%
flashDurs = cell2mat(respROIMat1(1).stimDat.FlashDuration);
grayDur = respROIMat1(1).stimDat.GrayDuration{1}; 
seqDur = flashDurs + grayDur;

cm = colormap('lines');
close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    iDur = pairedEpochs(i);
    seqDur = flashDurs(iDur) + grayDur;
    flashDur = flashDurs(iDur);

    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        hold on;
        rats1 = zeros(length(respROIMat1), length(respROIMat1(1).rats{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            rats1(n, :) = respROIMat1(n).rats{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(rats1,1);
        e1 = std(rats1,[],1)./sqrt(size(rats1,1));
        h1 = plot_err_patch_v2(...
            (0:length(m1))* respROIMat1(n).BIN_SHIFT,...
            [0 m1],...
            [e1(1) e1], ...
            cm(1,:),(cm(1,:)+1)/2);
        
        rats2 = zeros(length(respROIMat2), length(respROIMat2(1).rats{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            rats2(n, :) = respROIMat2(n).rats{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(rats2,1);
        e2 = std(rats2,[],1)./sqrt(size(rats2,1));
        h2 = plot_err_patch_v2(...
            (0:length(m2))* respROIMat2(n).BIN_SHIFT,...
            [0 m2],...
            [e2(1) e2], ...
            cm(2,:),(cm(2,:)+1)/2);
        
        rats3 = zeros(length(respROIMat3), length(respROIMat3(1).rats{pairedEpochs(i,j)}));
        flyID3 = zeros(1, length(respROIMat3));
        for n = 1:length(respROIMat3)
            rats3(n, :) = respROIMat3(n).rats{pairedEpochs(i,j)};
            flyID3(n) = respROIMat3(n).flyID;
        end         
        m3 = mean(rats3,1);
        e3 = std(rats3,[],1)./sqrt(size(rats3,1));        
        h3 = plot_err_patch_v2(...
            (0:length(m3))* respROIMat3(n).BIN_SHIFT,...
            [0 m3],...
            [e3(1) e3], ...
            cm(3,:),(cm(3,:)+1)/2);
       
        % plotting parameters
        xScale = [0, seqDur];
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        xlim(xScale);
        ylim(yScale);
        line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line
        
        if (inv)
            patchY = yScale(1);
        else
            patchY = yScale(2);
        end 
        
        if (inv)
            set(gca,'YDir','reverse');
        end
        
        if (j==2)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [1 1 1]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        elseif (j==1)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [0 0 0]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        end
        line([flashDurs(iDur) flashDurs(iDur)], yScale, 'color',[0 0 0],...
            'linestyle','--');
        legend([h1 h2 h3], ...
            [titleStr1 ', N cells = ' num2str(size(rats1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(rats2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        [titleStr3 ', N cells = ' num2str(size(rats3,1)) ...
        ', N flies = ' num2str(length(unique(flyID3)))],...
        'Location','southeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end

%% compare 2
flashDurs = cell2mat(respROIMat1(1).stimDat.FlashDuration);
grayDur = respROIMat1(1).stimDat.GrayDuration{1}; 
seqDur = flashDurs + grayDur;

cm = colormap('lines');
close all

for i = 1:size(pairedEpochs, 1) % number of flash lengths
    iDur = pairedEpochs(i);
    seqDur = flashDurs(iDur) + grayDur;
    flashDur = flashDurs(iDur);

    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        figure; 
        hold on;
        rats1 = zeros(length(respROIMat1), length(respROIMat1(1).rats{pairedEpochs(i,j)}));
        flyID1 = zeros(1, length(respROIMat1));
        for n = 1:length(respROIMat1)
            rats1(n, :) = respROIMat1(n).rats{pairedEpochs(i,j)};
            flyID1(n) = respROIMat1(n).flyID;
        end 
        m1 = mean(rats1,1);
        e1 = std(rats1,[],1)./sqrt(size(rats1,1));
        h1 = plot_err_patch_v2(...
            (0:length(m1))* respROIMat1(n).BIN_SHIFT,...
            [0 m1],...
            [e1(1) e1], ...
            cm(1,:),(cm(1,:)+1)/2);
        
        rats2 = zeros(length(respROIMat2), length(respROIMat2(1).rats{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            rats2(n, :) = respROIMat2(n).rats{pairedEpochs(i,j)};
            flyID2(n) = respROIMat2(n).flyID;
        end 
        m2 = mean(rats2,1);
        e2 = std(rats2,[],1)./sqrt(size(rats2,1));
        h2 = plot_err_patch_v2(...
            (0:length(m2))* respROIMat2(n).BIN_SHIFT,...
            [0 m2],...
            [e2(1) e2], ...
            cm(2,:),(cm(2,:)+1)/2);
        
       
        % plotting parameters
        xScale = [0, seqDur];
        xlabel('time (sec)');
        ylabel('response (dF/F)');
        xlim(xScale);
        ylim(yScale);
        line([0 seqDur],[0 0],'color',[0 0 0]); % x-axis line
        
        if (inv)
            patchY = yScale(1);
        else
            patchY = yScale(2);
        end 
        
        if (inv)
            set(gca,'YDir','reverse');
        end
        
        if (j==2)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [1 1 1]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        elseif (j==1)
            patch([0 0 flashDur flashDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)],...
                [0 0 0]);
            patch([flashDur flashDur flashDur+grayDur flashDur+grayDur],...
                [(patchY*0.85) (patchY*0.9) (patchY*0.9) (patchY*0.85)]...
                ,[0.5 0.5 0.5]);
        end
        line([flashDurs(iDur) flashDurs(iDur)], yScale, 'color',[0 0 0],...
            'linestyle','--');
        legend([h1 h2], ...
            [titleStr1 ', N cells = ' num2str(size(rats1,1)) ...
        ', N flies = ' num2str(length(unique(flyID1)))],...
        [titleStr2 ', N cells = ' num2str(size(rats2,1)) ...
        ', N flies = ' num2str(length(unique(flyID2)))],...
        'Location','southeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end