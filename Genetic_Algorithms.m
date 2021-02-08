clc;
clear all;
tic
filedir='.\TSP2.txt'; 
% filedir='.\Test.txt'; 

datain = textread(filedir, '', 'delimiter', ',');
[rows,columns]=size(datain);

totalrun=100000;

run = 0;
initialdist = twopoint(datain,rows,columns);
shortdist=initialdist;
shortset=datain;

summaryshort=zeros(totalrun+2,2);
summaryshort(1,1)=summaryshort(1,1)+run;
summaryshort(1,2)=summaryshort(1,2)+initialdist;

%%%%%%%%%%%CreateEvalutionMatrix%%%%%%%%%%%%%%%%%%%%%%%%%
evasize=50;%Must Even
ii=1;
evamatrix=zeros(evasize,rows+1);%+1 is distance
while ii<=evasize
    gene=randperm(rows);    
    evamatrix(ii,1:rows)=evamatrix(ii,1:rows) + gene;
    datanew=datain(evamatrix(ii,1:rows),:);
    datanewdist = twopoint(datanew,rows,columns);
    evamatrix(ii,rows+1)=evamatrix(ii,rows+1) + datanewdist;    
    ii=ii+1;   
end
evamatrix=sortrows(evamatrix,rows+1,'ascend');

endtempshort=zeros(evasize,totalrun+1);
endtempshort(:,1)=evamatrix(:,rows+1);
% endtempshort(1,1)+run;
% endtempshort(1,2)=endtempshort(1,2)+initialdist;

%%%%%%%%%%%FindShortDistance%%%%%%%%%%%%%%%%%%%%%%%%%
pickrows=round2even(randperm(30,1));%randperm(5,1)+5

GAdistahead = initialdist;
shortseq=1:rows;

while run<=totalrun   
    
    matbegin=evamatrix(1:pickrows,:);    
    GAmatrix=crossover(datain,matbegin,pickrows,rows,columns);
        
    evamatrix((evasize-pickrows+1):end,:)=GAmatrix;
    evamatrix=sortrows(evamatrix,rows+1,'ascend');    
    GAdistnow=evamatrix(1,rows+1);
    
    if GAdistnow<GAdistahead %shortest distance
        shortdist=GAdistnow;
        shortseq=evamatrix(1,1:rows);
    else
        shortdist=GAdistahead;      
    end     
    run=run+1;       
    endtempshort(:,run+1)=evamatrix(:,rows+1); 
    summaryshort(run+1,1)=summaryshort(run+1,1)+run;
    summaryshort(run+1,2)=summaryshort(run+1,2)+shortdist;         
    GAdistahead = shortdist;
end
shortset=datain(shortseq,:);  
% plot(shortset(:,1),shortset(:,2)); %plot the graph

goth=9; %run 4 times
savefilename=[filedir(1:end-4) 'GAShort' num2str(goth) '.mat'];
save(savefilename)

shortx=summaryshort(:,1);
shorty=summaryshort(:,2);
figure(1)%plot shortest distance
plot(shortx,shorty,'b')

toc

function GAmatrix=crossover(datain,evamatrix,defsize,rows,columns)    
    crosspoint=randperm(998,1);
    mutrange=randperm(998,2);
    mutstart=mutrange(1);
    mutend=mutrange(2);
    evacrossmat=zeros(defsize,rows+1);%+1 is distance
    ii=1;
    while ii<=defsize
        gene1=evamatrix(ii,1:rows);
        gene2=evamatrix(ii+1,1:rows); 
        temp1=gene1(1:crosspoint);    
        temp2=(ismember(gene2,temp1)-1)*(-1).*gene2;
        temp2=temp2(find(temp2~=0));
        newgene1=[temp1,temp2];        
        newgene1(:,mutstart:mutend)=fliplr(newgene1(:,mutstart:mutend));                        
        datanewgene1=datain(newgene1,:);           
        newgene1dist = twopoint(datanewgene1,rows,columns);           
        evacrossmat(ii,1:rows)=evacrossmat(ii,1:rows) + newgene1;
        evacrossmat(ii,rows+1)=evacrossmat(ii,rows+1) + newgene1dist;         
        temp3=gene1(crosspoint+1:end);
        temp4=(ismember(gene2,temp3)-1)*(-1).*gene2;
        temp4=temp4(find(temp4~=0));
        newgene2=[temp4,temp3];         
        newgene2=swape2row(newgene2,rows);
        datanewgene2=datain(newgene2,:);  
        newgene2dist = twopoint(datanewgene2,rows,columns);    
        evacrossmat(ii+1,1:rows)=evacrossmat(ii+1,1:rows) + newgene2;
        evacrossmat(ii+1,rows+1)=evacrossmat(ii+1,rows+1) + newgene2dist;    
        ii=ii+2;   
    end    
    GAmatrix=evacrossmat;
end

function doutput=swape2row(dinput,datar)
    rand2num=randperm(datar,2);
    temp=dinput(:,rand2num(1));
    dinput(:,rand2num(1))=dinput(:,rand2num(2));
    dinput(:,rand2num(2))=temp;
    doutput=dinput;
end 

function S = round2even(x) %cited from mathworks.com 
    if mod(x,2)<1 
        S = fix(x); 
    else 
        S =fix(x) + 1; 
    end
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
