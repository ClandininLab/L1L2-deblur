startInd = 5500;
numToAvg = 4000;
sep = 10000;

meanInten = zeros(1,numLevels);

for i=1:numLevels
    stInd = startInd + (i-1)*sep;
    endInd = stInd + numToAvg-1;
    meanInten(i) = mean(pdData(stInd:endInd));
end