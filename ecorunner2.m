function ret = ecorunner2()
clear;
clc;
num = 2;
N = [10 25 50 75 100 150 200 250 300 350 400 450 500];
for i = 1:length(N)
    eco2(N(i), num);
end

ret = true;
end