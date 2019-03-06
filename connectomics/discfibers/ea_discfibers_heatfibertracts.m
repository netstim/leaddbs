function [fibsweighted,fibsin,fibsval,iaix]=ea_discfibers_heatfibertracts(cfile,roilist,patselection,vals,connthreshold)
% function extracts fibers from a connectome connected to ROIs in the
% roilist and assigns them correlative values based on vals. Vals needs to be of
% same length as roilist, assigning a value for each ROI.

disp('ROI fiber analysis');

fibers=load(cfile);
fn=fieldnames(fibers);
try
    fibers=fibers.fibers;
catch
    fibers=fibers.(fn{1});
end

if ~iscell(roilist)
    roilist={roilist};
    vals={vals};
end

if ~exist('connthreshold','var') % minimum of percentage fibers that need to be connected to a VTA.
    connthreshold=0.2;
end

[fibsin,XYZmm,nii]=ea_discfibers_genroilist_connfibers(fibers, roilist, patselection);
fibsval=zeros(size(fibsin,1),length(patselection)); % 5th column will add up values, 6th will take note how many entries were summed.

% now color fibsin based on predictive value of improvement
ea_dispt('');

ea_dispercent(0,'Iterating ROI');
for roi=1:length(patselection)
    [~,D]=knnsearch(XYZmm{roi}(:,1:3),fibsin(:,1:3),'Distance','chebychev');
    in=D<mean(nii{end}.voxsize);
    fibsval(ismember(fibsin(:,4),unique(fibsin(in,4))),roi)=1;
    ea_dispercent(roi/length(patselection));
end
ea_dispercent(1,'end');

cnt=1;
ea_dispt('Correlating fibers with values');
for group=1:length(roilist)
    thisgroupidx=cnt:(cnt+length(patselection))-1;
    cnt=cnt+length(patselection);
    repvals=repmat(vals{group}',size(fibsval,1),1);
    try
        nfibsval=fibsval(:,thisgroupidx);
        nfibsval(~logical(nfibsval))=nan;
        posvals=repvals.*nfibsval;
        nfibsval=double(~fibsval(:,thisgroupidx));
        nfibsval(~logical(nfibsval))=nan;
        negvals=repvals.*nfibsval;
    catch
        keyboard
    end

    [~,fibidx,iaix]=unique(fibsin(:,4));
    fibsval=fibsval(fibidx,:);
    nfibsval=nfibsval(fibidx,:);

    sumfibsval=sum(fibsval,2);
    % discard fibers with less than connthreshold connections.
    exclude=sumfibsval<size(fibsval,2)*connthreshold;
    % discard fibers with more than 1-connthreshold connections.
    exclude=logical(exclude+(sumfibsval>size(fibsval,2)*(1-connthreshold)));
    fibsval(exclude,:)=[];
    fibsweighted=fibsin;
    fibsweighted((exclude(iaix)),:)=[];

    [~,fibidx,iaix]=unique(fibsweighted(:,4));

    allvals=repmat(vals{1}',size(fibsval,1),1);
    fibsimpval=allvals;
    fibsimpval(~logical(fibsval))=nan;
    nfibsimpval=allvals;
    nfibsimpval(logical(fibsval))=nan;
    [h,p,ci,stats]=ttest2(fibsimpval',nfibsimpval');
    fibsweighted=[fibsweighted,stats.tstat(iaix)'];
end
