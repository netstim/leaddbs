function ea_showdiscfibers(M,percent,resultfig)


patlist=M.patient.list(M.ui.listselect);
I=M.clinical.vars{M.ui.clinicallist}(M.ui.listselect);

N=length(patlist);
thresh=0.5; 
tic
reforce=0;
percent=0.2; % how many fibers at least need to be connected for fiber to be taken into account.
if reforce || ~exist([M.ui.groupdir,'correlative_fibertracts_',pointtodash(num2str(percent)),'.mat'],'file')
    for sub=1:length(patlist)
        roilist{sub,1}=[patlist{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_right.nii'];
        roilist{sub,2}=[patlist{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_left.nii'];
    end
    if ~exist([M.ui.groupdir,'connected_fibers.mat'],'file')
        cfile=[ea_getconnectomebase('dMRI'),M.ui.connectomename,filesep,'data.mat'];
    else
        cfile=[M.ui.groupdir,'connected_fibers.mat'];
    end
    [fibsweighted,fibsin]=ea_heatfibertracts(cfile,{roilist},{I},thresh,percent);
    save([M.ui.groupdir,'connected_fibers.mat'],'fibsin');
    save([M.ui.groupdir,'correlative_fibertracts_',pointtodash(num2str(percent)),'.mat'],'fibsweighted','-v7.3');
else
    load([M.ui.groupdir,'correlative_fibertracts_',pointtodash(num2str(percent))]);
end

% visualize:


if reforce || ~exist(['correlative_fibertracts_',pointtodash(num2str(percent)),'_reformatted.mat'],'file')
    fibidx=unique(fibsweighted(:,4));
    fibcell=cell(length(fibidx),1);
    valcell=fibcell;
    vals=zeros(length(fibidx),1);
    cnt=1;
    ea_dispercent(0,'reformatting fibers');
    for fib=fibidx';
        fibcell{cnt}=fibsweighted(fibsweighted(:,4)==fib,1:3);
        vals(cnt)=mean(fibsweighted(fibsweighted(:,4)==fib,5));
        cnt=cnt+1;
        ea_dispercent(cnt/length(fibidx));
    end
    ea_dispercent(1,'end');
    save(['correlative_fibertracts_',pointtodash(num2str(percent)),'_reformatted.mat'],'fibcell','vals','-v7.3');
else
    load(['correlative_fibertracts_',pointtodash(num2str(percent)),'_reformatted.mat']);
end



set(0,'CurrentFigure',resultfig);

%ea_mnifigure;
colormap gray
rb=ea_redblue;
cvals=vals;
cvals(isnan(cvals))=0;
cvals=cvals./max([abs(nanmin(cvals)),abs(nanmax(cvals))]);
alphas=abs(cvals);

% retain only 5% positive and 5% of negative predictive fibers:
posits=cvals(cvals>0);
negits=cvals(cvals<0);
posits=sort(posits,'descend');
negits=sort(negits,'ascend');
% posthresh=posits(round(length(posits)*0.05));
% negthresh=negits(round(length(negits)*0.05));
cutoff=500;
if (length(posits) < cutoff) || (length(negits) < cutoff)
    cutoff=min([length(posits),length(negits)]);
end
posthresh=posits(cutoff);
negthresh=negits(cutoff);
disp(['Fiber colors: Positive (T = ',num2str(posthresh),' - ',num2str(max(posits)),'); Negative (T = ',num2str(negthresh),' - ',num2str(min(negits)),').']);

remove=logical(logical(cvals<posthresh) .* logical(cvals>negthresh));
cvals(remove)=[];
fibcell(remove)=[];
cvals(cvals>0)=ea_minmax(cvals(cvals>0));
cvals(cvals<0)=-ea_minmax(-cvals(cvals<0));
alphas=abs(cvals);

cvals=cvals*31.5;
cvals=cvals+32.5;

cvals=rb(round(cvals),:);
% cvals=vals-min(vals);
% cvals=cvals/max(vals);
% cvals=cvals*63;
% cvals=cvals+1;
% cvals=rb(round(cvals),:);
h=streamtube(fibcell,0.2);
cv=mat2cell(cvals,ones(size(cvals,1),1));

% transform alphas to a logistic curve to highlight more predictive and
% suppress less predictive fibers:
alphas=1./...
    (1+...
    exp((-10)*...
    (alphas-0.5)));

alphas=alphas./nanmax(alphas);
%alphas(alphas<0.5)=alphas(alphas<0.5).^2;
%alphas=double(alphas>0.5);
calph=mat2cell(alphas,ones(size(cvals,1),1));

[h.FaceColor]=cv{:};
[h.FaceAlpha]=calph{:};

nones=repmat({'none'},size(cvals,1),1);
[h.EdgeColor]=nones{:};

function str=pointtodash(str)
str=strrep(str,'.','-');
