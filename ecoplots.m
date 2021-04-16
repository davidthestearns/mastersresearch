%this file uses the same variables as eco2.m and produces plots of each of
%the functions whose gradients can be found in eco2.m
clear;
clc;
x = [-3:.01:3];
numgroups =2;
time = 1000;
step =.05 ;
N = 100;
mobs = [1 1];
maxR = 8;
attMean = 1.0;
attNoise = .01*randi([90 110],numgroups);
attWeight = 4;
attSpread = .3;
segWeight =5;
segSpread = 1;
envWeight = 0;
gamma = [.2 .4 .6 .8 1];
minSep = 1;
avDist = 1;
j=1;
figure(1)
hold on
for i=1:5
   plot(x,segPot(x,segWeight,segSpread,N,gamma(i))); 
end
legend('.2','.4','.6','.8','1')
title('Segregation function G for different values of gamma')
hold off
figure(2)
hold on
for i=1:5
   plot(x,-attPot(x,segWeight,segSpread,N,gamma(i),attWeight,attSpread,attNoise,j)); 
end
legend('.2','.4','.6','.8','1')
title('Attractive function V - G for different values of gamma')
hold off
gamma = .5;
attspread = [.1 .2 .3 .4 .5];
figure(3)
hold on
for i=1:5
    plot(x,-attPot(x,segWeight,segSpread,N,gamma,attWeight,attspread(i),attNoise,j));
end
legend('.1','.2','.3','.4','.5')
title('Attractive function V - G for different values of b')
hold off
N = [20 40 60 80 100];
figure(4)
hold on
for i=1:5
  plot(x,-attPot(x,segWeight,segSpread,N(i),gamma,attWeight,attSpread,attNoise,j));  
end

legend('20','40','60','80','100');
title('Attractive function V - G for different values of N')

hold off
figure(5)
hold on
for i=1:5
  plot(x,segPot(x,segWeight,segSpread,N(i),gamma));  
end
legend('20','40','60','80','100');
title('Segregative function G for different values of N')
hold off


function seg = segPot(x,segWeight,segSpread,N,gamma)
    seg = N^(.5*gamma)*segWeight*exp(-(N^(gamma/2))*segSpread*(x.^2));
        
end
function att = attPot(x,segWeight,segSpread,N,gamma,attWeight,attSpread,attNoise,j)
    att = -segPot(x,segWeight,segSpread,N,gamma) + N^(.5*gamma)*attWeight*attNoise(j)*exp(-(N^(gamma/2))*attSpread*(x.^2));
        
end