clc;
clear all;
tic
filedir='.\TSP2.txt'; 

datain = textread(filedir, '', 'delimiter', ',');
[rows,columns]=size(datain);

totalrun=1000000;
%%%%%%%%%%%FindShortDistance%%%%%%%%%%%%%%%%%%%%%%%%%
run = 0;
initialdist = twopoint(datain,rows,columns);
shortdist=initialdist;
shortset=datain;
endtempshort=zeros(totalrun+2,2);
endtempshort(1,1)=endtempshort(1,1)+run;
endtempshort(1,2)=endtempshort(1,2)+initialdist;

summaryshort=zeros(totalrun+2,2);
summaryshort(1,1)=summaryshort(1,1)+run;
summaryshort(1,2)=summaryshort(1,2)+initialdist;

matahead=datain;

while run<=totalrun    
    HCdistahead = twopoint(matahead,rows,columns);       
    matnow=swape2row(matahead,rows);
    HCdistnow = twopoint(matnow,rows,columns);       
    if HCdistnow<HCdistahead %shortest distance
        shortdist=HCdistnow;
        shortset=matnow;
    else
        shortdist=HCdistahead;
        shortset=matahead;          
    end     
    run=run+1;    
    endtempshort(run+1,1)=endtempshort(run+1,1)+run;
    endtempshort(run+1,2)=endtempshort(run+1,2)+HCdistnow;          
    summaryshort(run+1,1)=summaryshort(run+1,1)+run;
    summaryshort(run+1,2)=summaryshort(run+1,2)+shortdist;      
    matahead=shortset;   %matafter    
end

%%%%%%%%%%%FindLongDistance%%%%%%%%%%%%%%%%%%%%%%%%%
run = 0;
initialdist = twopoint(datain,rows,columns);
longdist=initialdist;
longset=datain;
endtemplong=zeros(totalrun+2,2);
endtemplong(1,1)=endtemplong(1,1)+run;
endtemplong(1,2)=endtemplong(1,2)+initialdist;

summarylong=zeros(totalrun+2,2);
summarylong(1,1)=summarylong(1,1)+run;
summarylong(1,2)=summarylong(1,2)+initialdist;

matahead=datain;

while run<=totalrun    
    HCdistahead = twopoint(matahead,rows,columns);       
    matnow=swape2row(matahead,rows);
    HCdistnow = twopoint(matnow,rows,columns);       
    if HCdistnow>HCdistahead %longest distance
        longdist=HCdistnow;
        longset=matnow;
    else
        longdist=HCdistahead;
        longset=matahead;
    end     
    run=run+1;    
    endtemplong(run+1,1)=endtemplong(run+1,1)+run;
    endtemplong(run+1,2)=endtemplong(run+1,2)+HCdistnow;          
    summarylong(run+1,1)=summarylong(run+1,1)+run;
    summarylong(run+1,2)=summarylong(run+1,2)+longdist;  
    matahead=longset;   %matafter
end

goth=4; %run 4 times
savefilename=[filedir(1:end-4) 'HCgo' num2str(goth) '.mat'];
save(savefilename)

shortx=summaryshort(:,1);
shorty=summaryshort(:,2);
figure(1)%plot shortest distance
plot(shortx,shorty,'b')

longx=summarylong(:,1);
longy=summarylong(:,2);
figure(2)%plot longest distance
plot(longx,longy,'b')

toc

function doutput=swape2row(dinput,datar)
    rand2num=randperm(datar,2);
    temp=dinput(rand2num(1),:);
    dinput(rand2num(1),:)=dinput(rand2num(2),:);
    dinput(rand2num(2),:)=temp;
    doutput=dinput;
end 

function sumdist=twopoint(data,datar,datac)
    temp=zeros(datar,datac);
    tempdist=zeros(datar,1);
    temp(1:datar-1,:)=temp(1:datar-1,:)+data([2:end],:);
    temp(datar,:)=temp(datar,:)+data(1,:);
    deltadata=temp-data;
    tempdist=tempdist+((deltadata(:,1).^2)+(deltadata(:,2).^2)).^0.5;
    sumdist=sum(tempdist);
end 
