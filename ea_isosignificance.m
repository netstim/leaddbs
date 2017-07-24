function [p,idx]=ea_isosignificance(varargin)
% inputs: values at coordinates, threshold, max dist

vizz=0;
XYZV=varargin{1};
% clear NAN vars:
XYZV(isnan(XYZV(:,4)),:)=[];
XYZV=zscore(XYZV);
try
    thresh=varargin{2}; % standard deviation
catch
    thresh=1.5;
end
try
    mdist=mean(ea_pdist(XYZV(:,1:3)))*varargin{3}; % maximum distance to belong to a cluster
catch
    mdist=(mean(ea_pdist(XYZV(:,1:3)))/2); % maximum distance to belong to a cluster
end

%for pfit=-10:10
% init params:


% setup randomization:

for iter=1:50000
    
    dat=XYZV;
    if iter>1
        dat(:,4)=dat(randperm(length(dat)),4);
    end
    tdat=dat(dat(:,4)>thresh,:);
    distances=squareform(ea_pdist(tdat(:,1:3)));
    tdistances=triu(distances<mdist,1); % number of 1 in this matrix is size of cluster.
    csize(iter)=sum(tdistances(:));
    if iter==1
       idx=tdat; 
    end
end


[counts,centers]=hist(csize(2:end),1000);

[muhat,sigmahat]=normfit(csize(2:end));
pd=makedist('Normal',muhat,sigmahat);
y=pdf(pd,centers);
centers=round(centers);
ix=find(round(centers)==csize(1));
ix=ix(1);
p=trapz(y(ix:end))/trapz(y);

if p<0.05
    disp(['Significance at p=',num2str(p),'.']);
    if vizz
        
        figure
        h=area(linspace(centers(1),centers(end),length(y)),y,'LineStyle',':');
        h(1).FaceColor=[0,0,1];
        hold on
        j=area(linspace(centers(ix),centers(end),length(y(ix:end))),y(ix:end),'LineStyle',':');
        j(1).FaceColor=[1,0,0];
    end
else
    disp(['No significance. (p=',num2str(p),').']);
end



%end



% GM=fitgmdist(XYZV,1);
% P=posterior(GM,XYZV);
% idx=cluster(GM,XYZV);

%mdl = fitglm(XYZV(:,1:end-1),XYZV(:,end))