clear()
%variables. Most are passed into a struct to make passing them into the
%potential functions easier (see bottom of code)
p=struct('numgroups',{},'time',{},'step',{},'N',{},'attMean',{},'attNoise',{},'attWeight',{},'segWeight',{},'envWeight',{},'gamma',{},'minSep',{},'avDist',{},'maxR',{});
p(1).numgroups =2;
p(2).time =1000;
p(3).step=.1 ;
p(4).N=100;
mobs=[0 0];
maxR=8;
p(5).attMean=1.0;
p(6).attNoise=.1;
p(7).attWeight=0;
p(8).segWeight=7;
p(9).envWeight=1;
p(10).gamma=.1;
p(11).minSep=1;
p(12).avDist = 1;
p(13).maxR = 8;
thresh=0.05;
sted = 0;
numsted = 3;
i=1;
numgroups=2;
%need to implement periodic stuff, where distance is calculated wacky style
%close to the boundary

%reading in data and converting it to useable formats. Also defining output
%arrays.
xData = table2array(readtable('bigX.csv'));
yData = table2array(readtable('bigY.csv'));
allTimeX = zeros(p(2).time,length(xData(:,1)),length(xData(1,:)));
allTimeY = zeros(p(2).time,length(xData(:,1)),length(xData(1,:)));
allTimeX(1,:,:)=xData;
allTimeY(1,:,:)=yData;
%min(abs(maxR-intX(1,k,j)-intX(1,xk,xj)),abs(intX(1,k,j)-intX(1,xk,xj)))
%can parfor in the following loops: gamma values, group ID's, individuals
%this is the iterator through the potential functions. Once the individuals
%stabilize, or the time threshold is broken, the final values will output
%and animate. Can also be changed to output into a csv file later.
while sted<numsted
    intX=allTimeX(i,:,:);
    intY=allTimeY(i,:,:);
    newX = zeros(size(allTimeX(i,:,:)));
    newY = zeros(size(allTimeY(i,:,:)));
    for j=1:numgroups

        gsize = p(4).N;
        g=gsize;
        for k = 1:p(4).N
            segGrad = zeros(2);
            attGrad = zeros(2);
            envGrad = zeros(2);
            %totPot = zeros(2);
            for xj = 1:numgroups
                if xj==j
                    continue
                else
                    for xk = 1:p(4).N
                        %gseg = segPot((intX(1,k,j)-intX(1,xk,xj)),(intY(1,k,j)-intY(1,xk,xj)),p);
                        gseg = segPot(min([(maxR-intX(1,k,j)-intX(1,xk,xj)) (intX(1,k,j)-intX(1,xk,xj))]),min([(maxR-intY(1,k,j)-intY(1,xk,xj)) (intY(1,k,j)-intY(1,xk,xj))]),p);
                        if isnan(gseg(1)) == true
                            gseg(1)=100;
                        end
                        if isnan(gseg(2))==true
                            gseg(2)=100;
                        end
                        segGrad = segGrad +gseg;
                    end
                end
            end
            segGrad = -1*segGrad*(1/gsize);
            for ak=1:p(4).N
                if ak==k
                    continue
                else
                    %attG = attPot((intX(1,ak,j)-intX(1,k,j)),(intY(1,ak,j)-intY(1,k,j)),p);
                    attG = attPot(min([(maxR-intX(1,k,j)-intX(1,ak,j)) (intX(1,k,j)-intX(1,ak,j))]),min([(maxR-intY(1,k,j)-intY(1,ak,j)) (intY(1,k,j)-intY(1,ak,j))]),p);
                    if isnan(attG(1)) ==true
                        attG(1)=1;
                    end
                    if isnan(attG(2))==true
                        attG(2)=1;
                    attGrad = attGrad +attG;
                    end
                end
            end

            attGrad=-1*attGrad*(1/gsize);
            epo = envPot(intX(1,k,j),intY(1,k,j),p);
            if isnan(epo(1)) ==true
                epo(1)=0;
            end
            if isnan(epo(2))==true
                epo(2)= 0;
            end
            envGrad = envGrad+epo;
            totPot = segGrad+attGrad+mobs(j).*envGrad;
            checker = isnan(totPot);
            newX(1,k,j)=intX(1,k,j)+p(3).step*totPot(1);
            newY(1,k,j)=intY(1,k,j)+p(3).step*totPot(2);
            teX = newX(1,k,j);
            teY = newY(1,k,j);
            if teX>=maxR
                newX(1,k,j)=0+.5*rand();
            end
            if teX<=0
                newX(1,k,j)=maxR-.5*rand();
            end
            if teY>=maxR
                newY(1,k,j)=0+.5*rand();
            end
            if teY<=0
                newY(1,k,j)=maxR-.5*rand();
            end
        end
    end
    i=i+1;
    allTimeX(i,:,:)=newX;  
    allTimeY(i,:,:)=newY;
    dist= sqrt(sum(sum((allTimeX(i,:,:)-allTimeX(i-1,:,:)))).^2+sum(sum((allTimeY(i,:,:)-allTimeY(i-1,:,:)).^2)))
    if dist<thresh
        sted=sted+1;
    end
    if i>p(2).time
        break
    end
end
anim = simData(p(2).time,allTimeX,allTimeY);
%the code below is for animating the data. To view it in its entirety at 12
%FPS, write movie(ani) into the command line once the calculations are
%done. The animation will also play very slowly immediately.

 
%potential functions for the code above. These are prewritten as
%gradients, so just putting coordinates in them will output the proper
%values with no need for differentiation.
%check these things, probably messed up
function seg = segPot(x,y,p) 
    seg =-[(p(4).N)^(p(10).gamma)*p(12).avDist*p(8).segWeight*(p(4).N)^(.5*p(10).gamma)*(x)/abs(x)*exp((p(4).N)^(p(10).gamma/2)*(-abs(x+y)));
        (p(4).N)^(p(10).gamma)*p(12).avDist*p(8).segWeight*(p(4).N)^(.5*p(10).gamma)*(y)/abs(y)*exp((p(4).N)^(p(10).gamma/2)*(-abs(x+y)))]; 
end
function att = attPot(x,y,p)
    seg = [(p(4).N)^(p(10).gamma)*p(12).avDist*p(8).segWeight*(p(4).N)^(.5*p(10).gamma)*(x)/abs(x)*exp((p(4).N)^(p(10).gamma/2)*(-abs(x+y)));
        (p(4).N)^(p(10).gamma)*p(12).avDist*p(8).segWeight*(p(4).N)^(.5*p(10).gamma)*(y)/abs(y)*exp((p(4).N)^(p(10).gamma/2)*(-abs(x+y)))]; 
    att = -seg+ [(p(4).N)^(p(10).gamma)*p(7).attWeight*(p(4).N)^(.5*p(10).gamma)*0.3*(x)/abs(x)*exp(((p(4).N)^(.5*p(10).gamma))*(-0.3*abs(x+y))); 
    (p(4).N)^(p(10).gamma)*p(5).attMean*p(7).attWeight*(p(4).N)^(.5*p(10).gamma)*0.3*(y)/abs(y)*exp(((p(4).N)^(.5*p(10).gamma))*(-0.3*abs(x+y)))];
end
function env =envPot(x,y,p) 
    env = -[2*p(9).envWeight*.6*(x-.5*p(13).maxR)*exp(-.6*((x-.5*p(13).maxR)^2+(y-.5*p(13).maxR)^2)); 2*p(9).envWeight*.6*(y-.5*p(13).maxR)*exp(-.6*((x-.5*p(13).maxR)^2+(y-.5*p(13).maxR)^2))];
end
function sim = simData(time,allTimeX,allTimeY)
    ani = [];
    colors = ['r','g','b','y','m'];

    for i= 1:time


        if allTimeX(i,3,1) ~= 0
            hold on;
            %fig1 = axes('Xlim',[0 8],'Ylim',[0 8]);
            scatter(allTimeX(i,:,1),allTimeY(i,:,1),36,colors(1));
            scatter(allTimeX(i,:,2),allTimeY(i,:,2),36,colors(2));

            hold off;

            fra = getframe();
            ani = horzcat(ani,fra);
            cla();
        end
    end
    sim = ani;
end
    


 