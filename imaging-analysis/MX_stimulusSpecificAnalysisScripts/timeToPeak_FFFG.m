% timeToPeak_FFFG.m
% Time to peak calculations for L2 FFFG plots 
%
% Time to peak is defined as time to either the minimum or maximum value in
% the response. Note that this does not always yield the underlying time to
% peak when SNR is poor.
% 
% Marjorie Xie
% 6/9/2017

%% Single Genotype 
for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        TTP1 = zeros(length(respROIMat(:,2)),2); % times to peak matrix
        
        for r = 1:length(respROIMat(:,2)) % number of ROIs in the data matrix
            rats1 = [0, respROIMat(r,2).rats{pairedEpochs(i,j)}]; 
            timeArray = (0:length(rats1))*respROIMat(r,2).BIN_SHIFT;
            
            % Hyperpol peak 
            [minVal, iMin] = min(rats1);
            tMin = timeArray(iMin); 
            TTP1(r,1) = tMin;
            
            % Depol peak
            [maxVal, iMax] = max(rats1);
            tMax = timeArray(iMax);
            TTP1(r,2) = tMax;            
        end 
        
        if (j==1)
            lightTTP1 = TTP1;
        end 
        if (j==2)
            darkTTP1 = TTP1;
        end 
    end
    
    % distribution of time to peaks for light response
    figure;
    histogram(lightTTP1(:,1), 'BinWidth', respROIMat(r,2).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]); hold on; 
    histogram(lightTTP1(:,2), 'BinWidth', respROIMat(r,2).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]);  
    xlabel('time(sec)'); ylabel('n cells');
    legend('hyperpol peak', 'depol peak');
    title('light response time to peaks');
    
    % distribution of time to peaks for dark response
    figure;
    histogram(darkTTP1(:,1), 'BinWidth', respROIMat(r,2).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]); hold on; 
    histogram(darkTTP1(:,2), 'BinWidth', respROIMat(r,2).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]);  
    xlabel('time(sec)'); ylabel('n cells');
    legend('hyperpol peak', 'depol peak');
    title('dark response time to peaks');
end 

%% Genotype comparisons 
% Run necessary parts of plot_shortFlash_compareGenotypes.m first

%% Compare between 2 genotypes 
for i = 1:size(pairedEpochs, 1) % number of flash lengths
    for j = 1:size(pairedEpochs, 2) % number of contrast values per flash length
        TTP1 = zeros(length(respROIMat1), 2); % times to peak matrix, light on L, dark on R
        TTP2 = zeros(length(respROIMat2), 2); 
        for r = 1:length(respROIMat1) % number of ROIs in the data matrix
            rats1 = [0, respROIMat1(r).rats{pairedEpochs(i,j)}]; 
            timeArray = (0:length(rats1))*respROIMat1(r).BIN_SHIFT;
            % Hyperpol peak 
            [minVal, iMin] = min(rats1);
            tMin = timeArray(iMin); 
            TTP1(r,1) = tMin;
            % Depol peak
            [maxVal, iMax] = max(rats1);
            tMax = timeArray(iMax);
            TTP1(r,2) = tMax;            
        end 
        for r = 1:length(respROIMat2) % number of ROIs in the data matrix
            rats2 = [0, respROIMat2(r).rats{pairedEpochs(i,j)}]; 
            timeArray = (0:length(rats2))*respROIMat2(r).BIN_SHIFT;
            % Hyperpol peak 
            [minVal, iMin] = min(rats2);
            tMin = timeArray(iMin); 
            TTP2(r,1) = tMin;
            % Depol peak
            [maxVal, iMax] = max(rats2);
            tMax = timeArray(iMax);
            TTP2(r,2) = tMax;            
        end 
        
        if (j==1)
            lightTTP1 = TTP1;
            lightTTP2 = TTP2;
        end 
        if (j==2)
            darkTTP1 = TTP1;
            darkTTP2 = TTP2;
        end 
    end
    
    % hyperpol light response time to peak distribution 
    figure;
    histogram(lightTTP1(:,1), 'BinWidth', respROIMat1(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]); hold on; 
    histogram(lightTTP2(:,1), 'BinWidth', respROIMat2(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]);  
    xlabel('time(sec)'); ylabel('n cells');
    legend(titleStr1, titleStr2);
    title('light response time to peak hyperpol');
    
    % depol light response time to peak distribution 
    figure;
    histogram(lightTTP1(:,2), 'BinWidth', respROIMat1(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]); hold on; 
    histogram(lightTTP2(:,2), 'BinWidth', respROIMat2(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]);  
    xlabel('time(sec)'); ylabel('n cells');
    legend(titleStr1, titleStr2);
    title('light response time to peak depol');
    
    % hyperpol dark response time to peak distribution 
    figure;
    histogram(darkTTP1(:,1), 'BinWidth', respROIMat1(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]); hold on; 
    histogram(darkTTP2(:,1), 'BinWidth', respROIMat2(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]);  
    xlabel('time(sec)'); ylabel('n cells');
    legend(titleStr1, titleStr2);
    title('dark response time to peak hyperpol');
    
    % depol dark response time to peak distribution 
    figure;
    histogram(darkTTP1(:,2), 'BinWidth', respROIMat1(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]); hold on; 
    histogram(darkTTP2(:,2), 'BinWidth', respROIMat2(1).BIN_SHIFT, ...
        'BinLimits',[0, timeArray(end)]);  
    xlabel('time(sec)'); ylabel('n cells');
    legend(titleStr1, titleStr2);
    title('dark response time to peaks depol');
end 

