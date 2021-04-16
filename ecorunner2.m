%call this file in any batch files on the supercomputer
function ret = ecorunner2()
clear;
clc;
num = 1; %number of groups in simulations
N = [10 15 20 25 50 75 100 150 200 300 400]; %number of individuals in each group
for i = 1:length(N)
    eco2(N(i), num);
end

ret = true;
end