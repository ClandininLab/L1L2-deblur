% plot_shortFlash_compareGenotypes.m
%
% quick and dirty script to plot short flash responses from multiple
% genotypes on top of each other
%
% for L2 project - cell type specific silencing
%
% 10/22/17


%% L2 cell type-specific silencing - R11D03AD, R19C10DBD, UAS-TNT (Lawf2)
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/171021_celltypespecific';
cd(dataPath)

yScale = [-0.05 0.07];
pairedEpochs = [1 2];
inv = 1;

% Lawf2>>TNT
saveName1 = 'L2lexA_ASAP2f_Lawf2-TNT';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'Lawf2>>TNT';

% Lawf2-Gal4
saveName2 = 'L2lexA_ASAP2f_Lawf2';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'Lawf2-Gal4 control';

% TNT
saveName3 = 'L2lexA_ASAP2f_TNT';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'TNT control';

%% L2 cell type-specific silencing - R52H01AD, R17C11DBD, UAS-TNT (Lawf1)
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/171021_celltypespecific';
cd(dataPath)

yScale = [-0.05 0.07];
pairedEpochs = [1 2];
inv = 1;

% Lawf1>>TNT
saveName1 = 'L2lexA_ASAP2f_Lawf1-TNT';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'Lawf1>>TNT';

% Lawf1-Gal4
saveName2 = 'L2lexA_ASAP2f_Lawf1';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'Lawf1-Gal4 control';

% TNT
saveName3 = 'L2lexA_ASAP2f_TNT';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'TNT control';

%% L2 cell type-specific silencing - R20C11AD, R48D11DBD UAS-TNT (C2/C3)
clear all
close all
clc

dataPath = '/Users/hyang/Documents/New Imaging Analysis/AnalyzedData/171021_celltypespecific';
cd(dataPath)

yScale = [-0.05 0.07];
pairedEpochs = [1 2];
inv = 1;

% C2/C3>>TNT
saveName1 = 'L2lexA_ASAP2f_C2C3-TNT';
load(saveName1, 'roiDataMat','iResp');
respROIMat1 = roiDataMat(iResp,2);
titleStr1 = 'C2C3>>TNT';

% C2/C3-Gal4
saveName2 = 'L2lexA_ASAP2f_C2C3';
load(saveName2, 'roiDataMat','iResp');
respROIMat2 = roiDataMat(iResp,2);
titleStr2 = 'C2C3-Gal4 control';

% TNT
saveName3 = 'L2lexA_ASAP2f_TNT';
load(saveName3, 'roiDataMat','iResp');
respROIMat3 = roiDataMat(iResp,2);
titleStr3 = 'TNT control';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% compare 2, with individual cells
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
%         for n = 1:length(respROIMat1)
%             rats1(n, :) = respROIMat1(n).rats{pairedEpochs(i,j)};
%             plot(respROIMat1(n).t{pairedEpochs(i,j)},rats1(n,:),...
%                 'Color', (cm(1,:)+1)/2);
%             flyID1(n) = respROIMat1(n).flyID;
%         end 
%         m1 = mean(rats1,1);
%         e1 = std(rats1,[],1)./sqrt(size(rats1,1));
%         h1 = plot_err_patch_v2(...
%             respROIMat1(1).t{pairedEpochs(i,j)},...
%             m1 ,...
%             e1, ...
%             cm(1,:),(cm(1,:)+1)/2);
        
        rats2 = zeros(length(respROIMat2), length(respROIMat2(1).rats{pairedEpochs(i,j)}));
        flyID2 = zeros(1, length(respROIMat2));
        for n = 1:length(respROIMat2)
            rats2(n, :) = respROIMat2(n).rats{pairedEpochs(i,j)};
            plot(respROIMat2(n).t{pairedEpochs(i,j)},rats2(n,:),...
                'Color', (cm(2,:)+1)/2);
            flyID2(n) = respROIMat2(n).flyID;
        end 
%         m2 = mean(rats2,1);
%         e2 = std(rats2,[],1)./sqrt(size(rats2,1));
%         h2 = plot_err_patch_v2(...
%             respROIMat2(1).t{pairedEpochs(i,j)},...
%             m2,...
%             e2, ...
%             cm(2,:),(cm(2,:)+1)/2);
       
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
%         legend([h1 h2], ...
%             [titleStr1 ', N cells = ' num2str(size(rats1,1)) ...
%         ', N flies = ' num2str(length(unique(flyID1)))],...
%         [titleStr2 ', N cells = ' num2str(size(rats2,1)) ...
%         ', N flies = ' num2str(length(unique(flyID2)))],...
%         'Location','southeast');
    
        if (j==2)
            title('Light flash');
        elseif (j==1)
            title('Dark flash');
        end
    end 

end