function [fibsweighted,fibsin,fibsval,iaix]=ea_discfibers_heatfibertracts_corr(cfile,roilist,patselection,vals,efieldthresh)
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

[fibsin,XYZmm,nii,valsmm]=ea_discfibers_genroilist_connfibers(fibers, roilist, patselection, efieldthresh);
fibsval=zeros(size(fibsin,1),length(patselection)); % 5th column will add up values, 6th will take note how many entries were summed.

% now color fibsin based on predictive value of improvement
ea_dispt('');
ea_dispercent(0,'Iterating ROI');
fibunique=unique(fibsin(:,4))';
for roi=1:length(patselection)
    [IX,D]=knnsearch(XYZmm{roi}(:,1:3),fibsin(:,1:3),'Distance','chebychev');
    in=D<mean(nii{end}.voxsize);
    for fib=fibunique
        fibsel=fibsin(:,4)==fib;
        fibsval(fibsel,roi)=sum(valsmm{roi}(IX(and(in,fibsel))));
    end
    ea_dispercent(roi/length(patselection));
end
ea_dispercent(1,'end');
ea_dispt('Correlating fibers with values');
cnt=1;
for group=1:length(roilist) % groups currently not implemented, should always be one within Lead-DBS
    % thisgroupidx=cnt:(cnt+length(patselection))-1;
    cnt=cnt+length(patselection);

    % reduce to one entry per fiber:
    [~,fibidx,iaix]=unique(fibsin(:,4));
    fibsval=fibsval(fibidx,:);
    %repvals=repmat(vals{group}',size(fibsval,1),1);

    [R]=corr(vals{group},fibsval','rows','pairwise','type','Spearman');

    fibsweighted=fibsin;

    fibsweighted=[fibsweighted,R(iaix)'];
end
