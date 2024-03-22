% test key press
clear all

counter = 0;
keyCode = zeros(1,256);
escCode = zeros(1,256);
escCode(27) = 1;
keydown = 0;

while ~isequal(escCode,keyCode)
    [~, ~, keyCode, ~] = KbCheck;
    pause(0.01);
    counter = counter + 1;
end

while ~keydown
    [keydown, ~, keyCode, ~] = KbCheck;
    pause(0.01);
    counter = counter + 1;
end

