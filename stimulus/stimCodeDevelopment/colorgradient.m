numCol = 64;
col = fliplr(linspace(0,63/255,numCol))';
% gcol = [zeros(size(col)); ones(size(col))*64; ones(size(col))*128; ...
%     ones(size(col))*192]/255;
% col = [col; col; col; col];
myColMap = horzcat(zeros(size(col)), col, zeros(size(col)));
% myColMap = horzcat(col, zeros(size(col)), zeros(size(col)));
% myColMap = horzcat(col,gcol,zeros(size(col)));
% myColMap = horzcat(zeros(size(col)),zeros(size(col)),col);
imagesc(0:(numCol-1));
colormap(myColMap)


x = uint8(63:-1:0);
y = bitshift(x,2);



colb = double(y)' / 255;
numCol = 64;
myColMap = horzcat(zeros(size(col)),col,colb);
figure; imagesc(0:(numCol-1));
colormap(myColMap)