function [fibcell,fibsval,XYZmm,nii,valsmm]=ea_discfibers_heatfibertracts_corr(cfile,roilist,patselection,vals,efieldthresh)
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

% now color fibsin based on predictive value of improvement
ea_dispt('');


% Reformat to cell:
[~,fibiaxfirst]=unique(fibsin(:,4),'first');
[~,fibiaxlast]=unique(fibsin(:,4),'last');
fiblen = fibiaxlast - fibiaxfirst + 1;
fibcell = mat2cell(fibsin(:,1:3),fiblen);

% repair fibsin to be incrementing from 1 to x:
for f=1:length(fibcell)
    fibsin(fibiaxfirst(f):fibiaxlast(f),4)=f;
end
    

dthresh=2*mean(nii{end}.voxsize);
for side=1:2
    fibsval{side}=zeros(size(fibcell,1),length(patselection)); % 5th column will add up values, 6th will take note how many entries were summed.
    
    ea_dispercent(0,['Iterating ROI, side ',num2str(side)]);   
    for roi=1:length(patselection)
        tree=KDTreeSearcher(XYZmm{roi,side}(:,1:3));
        for fib=1:length(fibcell)
            [IX,D]=knnsearch(tree,fibcell{fib},'Distance','chebychev');
            in=D<dthresh;
            if any(in)
                fibsval{side}(fib,roi)=sum(valsmm{roi,side}(IX(in)));
            end
        end
        ea_dispercent(roi/length(patselection));
    end
    ea_dispercent(1,'end');
    
end

%ea_dispt('Correlating fibers with values');
% cnt=1;
% for group=1:length(roilist) % groups currently not implemented, should always be one within Lead-DBS
%     % thisgroupidx=cnt:(cnt+length(patselection))-1;
%     cnt=cnt+length(patselection);
% 
%     % reduce to one entry per fiber:
%     [~,fibidx,iaix]=unique(fibsin(:,4));
%     fibsval=fibsval(fibidx,:);
%     %repvals=repmat(vals{group}',size(fibsval,1),1);
% 
%     [R]=corr(vals{group},fibsval','rows','pairwise','type','Spearman');
% 
%     fibsweighted=fibsin;
% 
%     fibsweighted=[fibsweighted,R(iaix)'];
% end
