% fft_test.m
% script to test fft 
t = respROIMat1(1).t{1};
x = respROIMat1(1).rats{1};
x = m1;

sampRate = 120;     % sampling rate in Hz
n = 1024;
y = fft(x,n);                               % Compute DFT of x
m = abs(y);                               % Magnitude
p = unwrap(angle(y));                     % Phase

f = (0:length(y)-1)*sampRate/length(y);        % Frequency vector

p1 = 1:(floor(length(y)/2)+1);             % indicies of first half

% sgf = sgolayfilt(m,3,21);

% f = smooth(x',0.05,'lowess');
% f = smooth(x',3);
smoX = smooth(x',7,'sgolay',5);
figure; 
plot(t,x)
hold on;
plot(t,smoX)

windowSize = 3;
boxFn = (1/windowSize)*ones(1,windowSize);
filtX = filter(boxFn,1,x);
filtY = fft(filtX,n);
filtM = abs(filtY);                               % Magnitude
filtP = unwrap(angle(filtY));                     % Phase
f2 = (0:length(filtY)-1)*120/length(filtY);  


figure;
subplot(2,1,1)
% loglog(f(p1),m(p1));
plot(f(p1),m(p1));
hold on;
% loglog(f2(p1),filtM(p1),'r');
title('Magnitude')



subplot(2,1,2)
plot(f(p1),p(p1)*180/pi)
hold on;
% loglog(f2(p1),filtP(p1)*180/pi,'r');
title('Phase')

