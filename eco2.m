function ret =eco2(n,num)
%initializing variables with global scope to distinguish them
global numgroups N maxR attMean attNoise attWeight segWeight envWeight gamma1 gamma2 minSep avDist time step mobs thresh sted numsted i attSpread segSpread
numgroups =num;
time = 6000; %time*step =total time
step =.1 ;
N = n;
mobs = [1 1]; %environmental mobility
maxR = 7; %size of the territory, don't alter
attMean = 1.0;
attNoise = .01*randi([90 110],numgroups); %sigma from paper
attWeight = 7; %R from paper
attSpread = .4; %b from paper
segWeight =4; %S from paper
segSpread = 1; %a from paper
envWeight = 3; %weight this agains R and S so forces are comproble
gamma1 = 0;
gamma2 = .5; %nonlinear aggregation coefficient
minSep = 1;
avDist = 1;
%thresh = 0.119;
thresh = .0001; %threshold for steady state
sted = 0;
numsted = 3;
i = 1;
eBool = false; %false for contin env, true for data based env
%reading in data and converting it to useable formats. Also defining output
%arrays.
%The included files bigX and bigY are 1000x2 arrays. Can use the tool
%mockaroo on the internet to create differently sized datasets. Make the
%number of rows be number of individuals, and number of columns numgroups
xData = table2array(readtable('bigX.csv'));
yData = table2array(readtable('bigY.csv'));
sand = readtable('SAND.csv');
edge = readtable('EDGE.csv');
grass = readtable('GRASS.csv');
dem = readtable('DEM.csv');
sand = table2array(sand(:,2:end));
edge = table2array(edge(:,2:end));
grass = table2array(grass(:,2:end));
dem = table2array(dem(:,2:end));
[sX,sY] = gradient(sand);
[eX,eY] = gradient(edge);
[gX,gY] = gradient(grass);
[dX,dY] = gradient(dem);
envX = sX+eX+gX+dX;
envY = sY+eY+gY+dY;
%output arrays below
allTimeX = zeros(time,N,numgroups);
allTimeY = zeros(time,N,numgroups);
if numgroups>1
allTimeX(1,:,:)=xData(1:N,:);
allTimeY(1,:,:)=yData(1:N,:);
else
allTimeX(1,:,:)=xData(1:N,1);
allTimeY(1,:,:)=yData(1:N,1);
end
%iterating code for calculating motion of particles
while sted<numsted
    oldX = allTimeX(i,:,:);
    oldY = allTimeY(i,:,:);
    potsX = zeros(size(oldX));
    potsY = zeros(size(oldY));
    %variable timestep, optional
    if i>5*time/6 && i<5*time/6 +1
        step = step -.01;
    end
    for j = 1:numgroups
        for k = 1:N
            attGrad = zeros(2);
            segGrad = zeros(2);
            %environmental potential
            if eBool == false
                envGrad = envPot(oldX(1,k,j),oldY(1,k,j),envWeight,maxR);
            end
            if eBool ==true
                xpart = round(11*oldX(1,k,j));
                ypart = round(11*oldY(1,k,j));
                if xpart<=0
                    xpart = 1;
                end
                if ypart<=0
                    ypart = 1;
                end
                if xpart>77
                    xpart = 77;
                end
                if ypart>77
                    ypart = 77;
                end
                envGrad = envWeight*[envX(xpart,ypart);envY(xpart,ypart)];
            end
            %segregation potential
            for aj = 1:numgroups
               if aj == j
                   continue;
               end
               for ak = 1:N
                   %includes periodic boundary conditions
                   %mins and max compares dist across boundary and within
                   %boundary
                   num = segPot(min([(maxR-max(oldX(1,k,j),oldX(1,ak,aj))+min(oldX(1,k,j),oldX(1,ak,aj))) (oldX(1,k,j)-oldX(1,ak,aj))]), min([(maxR-max(oldY(1,k,j),oldY(1,ak,aj))+min(oldY(1,k,j),oldY(1,ak,aj))) (oldY(1,k,j)-oldY(1,ak,aj))]),segWeight,segSpread,N,gamma2);
                   segGrad= segGrad + num;                   
               end
            end
            segGrad = -1*segGrad/N;
            %attractive potential
            for ak =1:N
                if ak ==k
                    continue
                end
                num1 = attPot(min([(maxR-max(oldX(1,k,j),oldX(1,ak,j))+min(oldX(1,k,j),oldX(1,ak,j))) (oldX(1,k,j)-oldX(1,ak,j))]), min([(maxR-max(oldY(1,k,j),oldY(1,ak,j))+min(oldY(1,k,j),oldY(1,ak,j))) (oldY(1,k,j)-oldY(1,ak,j))]),N,gamma1,attWeight,attSpread,attNoise,j);
                num2 = segPot(min([(maxR-max(oldX(1,k,j),oldX(1,ak,j))+min(oldX(1,k,j),oldX(1,ak,j))) (oldX(1,k,j)-oldX(1,ak,j))]), min([(maxR-max(oldY(1,k,j),oldY(1,ak,j))+min(oldY(1,k,j),oldY(1,ak,j))) (oldY(1,k,j)-oldY(1,ak,j))]),segWeight,segSpread,N,gamma2);
                num = -num1+num2;
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
    %this keeps all the data inside the data set
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
%animation code
filename = 'onegroupnoenvdiffgamcontenv' + string(N) + '_individuals';
anim = simData(time,allTimeX,allTimeY,filename);
hold on
%have a number of scatter calls here equal to numgroups
scatter(allTimeX(1,:,1),allTimeY(1,:,1),'r')
scatter(allTimeX(1,:,2),allTimeY(1,:,2),'g')
%title('Initial Positions of One Group')
hold off
%functions for visualizing data and calculating potentials
function sim = simData(time,allTimeX,allTimeY,filename)
    ani = [];
    colors = ['r','g','b','y','m'];
    v = VideoWriter(filename,'Motion JPEG AVI');
    open(v)
    for i= 1:time

        if allTimeX(i,3,1) ~= 0
            hold on;
            xlim([0 7]);
            ylim([0 7]);
            %have a number of scatter calls here equal to numgroups
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
ret = anim;
%potential functions, already expressed as gradients so one doesn't have to
%be computed
function seg = segPot(x,y,segWeight,segSpread,N,gamma)
    seg = -[N^(1.5*gamma)*segWeight*(2*x)*segSpread*exp(-(N^(gamma/2))*segSpread*(x^2+y^2));
        N^(1.5*gamma)*segWeight*(2*y)*segSpread*exp(-(N^(gamma/2))*segSpread*(x^2+y^2))];
end
function att = attPot(x,y,N,gamma,attWeight,attSpread,attNoise,j)
    att = -[N^(1.5*gamma)*attWeight*attNoise(j)*(2*x)*attSpread*exp(-(N^(gamma/2))*attSpread*(x^2+y^2));
        N^(1.5*gamma)*attWeight*attNoise(j)*(2*y)*attSpread*exp(-(N^(gamma/2))*attSpread*(x^2+y^2))];
end
function env = envPot(x,y,envWeight,maxR)
    env = envWeight*[2*(x-maxR/3)*exp(-.8*((x-maxR/3)^2+(y-maxR/3)^2)); 2*(y-maxR/3)*exp(-.8*((x-maxR/3)^2+(y-maxR/3)^2))]; 
end
end
