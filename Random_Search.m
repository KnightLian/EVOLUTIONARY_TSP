clc;
clear all;
tic
filedir='.\TSP2.txt'; 

datain = textread(filedir, '', 'delimiter', ',');
[rows,columns]=size(datain);

initialdist = twopoint(datain,rows,columns);
longdist=initialdist;
longset=datain;
shortdist=initialdist;
shortset=datain;

run = 0;
totalrun=1000000;
summarytemp=zeros(totalrun+2,2);
summarytemp(1,1)=summarytemp(1,1)+run;
summarytemp(1,2)=summarytemp(1,2)+initialdist;

summaryshort=zeros(totalrun+2,2);
summaryshort(1,1)=summaryshort(1,1)+run;
summaryshort(1,2)=summaryshort(1,2)+initialdist;

summarylong=zeros(totalrun+2,2);
summarylong(1,1)=summarylong(1,1)+run;
summarylong(1,2)=summarylong(1,2)+initialdist;

while run<=totalrun       
    data_random=datain(randperm(rows),:);
    [rows,columns]=size(data_random);
    randomdist = twopoint(data_random,rows,columns);    
    if randomdist>longdist %largest distance
        longdist=randomdist;
        longset=data_random;
    end    
    if randomdist<shortdist %short distance
        shortdist=randomdist;
        shortset=data_random;
    end                   
    run=run+1;     
    summarytemp(run+1,1)=summarytemp(run+1,1)+run;
    summarytemp(run+1,2)=summarytemp(run+1,2)+randomdist;         
    summaryshort(run+1,1)=summaryshort(run+1,1)+run;
    summaryshort(run+1,2)=summaryshort(run+1,2)+shortdist;      
    summarylong(run+1,1)=summarylong(run+1,1)+run;
    summarylong(run+1,2)=summarylong(run+1,2)+longdist;    
end

goth=4; %run 4 times
savefilename=[filedir(1:end-4) 'randomgo' num2str(goth) '.mat'];
save(savefilename)

recorddistx=summarytemp(:,1);
recorddisty=summarytemp(:,2);

shortx=summaryshort(:,1);
shorty=summaryshort(:,2);
figure(1)%plot shortest distance
plot(shortx,shorty,'b')

longx=summarylong(:,1);
longy=summarylong(:,2);
figure(2)%plot longest distance
plot(longx,longy,'b')

toc

function sumdist=twopoint(data,datar,datac)
    temp=zeros(datar,datac);
    tempdist=zeros(datar,1);
    temp(1:datar-1,:)=temp(1:datar-1,:)+data([2:end],:);
    temp(datar,:)=temp(datar,:)+data(1,:);
    deltadata=temp-data;
    tempdist=tempdist+((deltadata(:,1).^2)+(deltadata(:,2).^2)).^0.5;
    sumdist=sum(tempdist);
end 
