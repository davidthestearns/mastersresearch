function ret = eco2(n,num)
%initializing variables with global scope
global numgroups N maxR attMean attNoise attWeight segWeight envWeight gamma minSep avDist time step mobs thresh sted numsted i attSpread segSpread
numgroups =num;
time = 3000;
step =.04 ;
N = n;
mobs = [1 1];
maxR = 7;
attMean = 1.0;
attNoise = .01*randi([90 110],numgroups);
attWeight = 4;
attSpread = .3;
segWeight =5;
segSpread = 1;
envWeight = 3;
gamma = 0.5;
minSep = 1;
avDist = 1;
%thresh = 0.119;
thresh = .001;
sted = 0;
numsted = 3;
i = 1;
%reading in data and converting it to useable formats. Also defining output
%arrays.
xData = table2array(readtable('bigX.csv'));
yData = table2array(readtable('bigY.csv'));
allTimeX = zeros(time,N,numgroups);
allTimeY = zeros(time,N,numgroups);
allTimeX(1,:,:)=xData(1:N,:);
allTimeY(1,:,:)=yData(1:N,:);
%iterating code for calculating motion of particles
while sted<numsted
    oldX = allTimeX(i,:,:);
    oldY = allTimeY(i,:,:);
    potsX = zeros(size(oldX));
    potsY = zeros(size(oldY));
    if i>5*time/6 && i<5*time/6 +1
        step = step -.01;
    end
    for j = 1:numgroups
        for k = 1:N
            attGrad = zeros(2);
            segGrad = zeros(2);
            %environmental potential
            envGrad = envPot(oldX(1,k,j),oldY(1,k,j),envWeight,maxR);
            %segregation potential
            for aj = 1:numgroups
               if aj == j
                   continue;
               end
               for ak = 1:N
                   %includes periodic boundary conditions
                   num = segPot(min([(maxR-max(oldX(1,k,j),oldX(1,ak,aj))+min(oldX(1,k,j),oldX(1,ak,aj))) (oldX(1,k,j)-oldX(1,ak,aj))]), min([(maxR-max(oldY(1,k,j),oldY(1,ak,aj))+min(oldY(1,k,j),oldY(1,ak,aj))) (oldY(1,k,j)-oldY(1,ak,aj))]),segWeight,segSpread,N,gamma);
                   segGrad= segGrad + num;                   
               end
            end
            segGrad = -1*segGrad/N;
            %attractive potential
            for ak =1:N
                if ak ==k
                    continue
                end
                num = attPot(min([(maxR-max(oldX(1,k,j),oldX(1,ak,j))+min(oldX(1,k,j),oldX(1,ak,j))) (oldX(1,k,j)-oldX(1,ak,j))]), min([(maxR-max(oldY(1,k,j),oldY(1,ak,j))+min(oldY(1,k,j),oldY(1,ak,j))) (oldY(1,k,j)-oldY(1,ak,j))]),segWeight,segSpread,N,gamma,attWeight,attSpread,attNoise,j);
                attGrad = attGrad + num;
            end
            attGrad = -1*attGrad/N;
            totPot = attGrad + segGrad - mobs(j)*envGrad;
            if isnan(totPot) ==true
                totPot = rand();
            end
            potsX(1,k,j) = totPot(1);
            potsY(1,k,j)= totPot(2);        

        end
    end
    newX = oldX + step .*potsX;
    newY = oldY + step .*potsY;
    for j=1:numgroups
        for k=1:N
            xval = newX(1,k,j);
            yval = newY(1,k,j);
            if xval<0
                newX(1,k,j)=maxR-.5*rand();
            end
            if xval>maxR
                newX(1,k,j)=0+.5*rand();
            end
            if yval<0
                newY(1,k,j)=maxR-.5*rand();
            end
            if yval>maxR
                newY(1,k,j)=0+.5*rand();
            end
        end
    end
    i=i+1;
    allTimeX(i,:,:)=newX;
    allTimeY(i,:,:)=newY;
    %calculation of steady state
    dist= sqrt(sum(sum((allTimeX(i,:,:)-allTimeX(i-1,:,:)))).^2+sum(sum((allTimeY(i,:,:)-allTimeY(i-1,:,:)).^2)));

    if dist<thresh
        sted=sted+1;
    end
    if i>time
        break
    end
end
filename = 'twogroupenv_' + string(N) + '_individuals.mp4';
anim = simData(time,allTimeX,allTimeY,filename);
hold on
scatter(allTimeX(1,:,1),allTimeY(1,:,1),'r')
scatter(allTimeX(1,:,2),allTimeY(1,:,2),'g')
%title('Initial Positions of One Group')
hold off
%functions for visualizing data and calculating potentials
function sim = simData(time,allTimeX,allTimeY,filename)
    ani = [];
    colors = ['r','g','b','y','m'];
    v = VideoWriter(filename,'Motion JPEG 2000');
    open(v)
    for i= 1:time

        if allTimeX(i,3,1) ~= 0
            hold on;
            xlim([0 7]);
            ylim([0 7]);
            scatter(allTimeX(i,:,1),allTimeY(i,:,1),36,colors(1));
            scatter(allTimeX(i,:,2),allTimeY(i,:,2),36,colors(2));


            hold off;

            fra = getframe();
            writeVideo(v,fra);
            ani = horzcat(ani,fra);
            
            cla();
        end
    end
    close(v);
    sim = ani;
end
function seg = segPot(x,y,segWeight,segSpread,N,gamma)
    seg = -[N^(1.5*gamma)*segWeight*(2*x)*segSpread*exp(-(N^(gamma/2))*segSpread*(x^2+y^2));
        N^(1.5*gamma)*segWeight*(2*y)*segSpread*exp(-(N^(gamma/2))*segSpread*(x^2+y^2))];
end
function att = attPot(x,y,segWeight,segSpread,N,gamma,attWeight,attSpread,attNoise,j)
    att = segPot(x,y,segWeight,segSpread,N,gamma) + [N^(1.5*gamma)*attWeight*attNoise(j)*(2*x)*attSpread*exp(-(N^(gamma/2))*attSpread*(x^2+y^2));
        N^(1.5*gamma)*attWeight*attNoise(j)*(2*y)*attSpread*exp(-(N^(gamma/2))*attSpread*(x^2+y^2))];
end
function env = envPot(x,y,envWeight,maxR)
    env = envWeight*[2*(x-maxR/6)*exp(-.8*((x-maxR/6)^2+(y-maxR/6)^2)); 2*(y-maxR/6)*exp(-.8*((x-maxR/6)^2+(y-maxR/6)^2))]; 
end
end
