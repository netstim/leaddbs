function [fibcell,fibsval,XYZmm,nii]=ea_discfibers_heatfibertracts(cfile,roilist,patselection,vals,connthreshold)
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
    
searchfibsin=round(fibsin(:,1:3).*4)/4; % reduce precision a bit (0.25 mm) to speed up search
[searchfibsin,~,ic]=unique(searchfibsin,'rows');

for side=1:2
    fibsval{side}=zeros(size(fibcell,1),length(patselection)); % 5th column will add up values, 6th will take note how many entries were summed.
    
    % now color fibsin based on predictive value of improvement
    
    ea_dispercent(0,['Iterating ROI, side ',num2str(side)]);
 
    for roi=1:length(patselection)
        [~,D]=knnsearch(XYZmm{roi,side}(:,1:3),searchfibsin,'Distance','chebychev');
        in=D<mean(nii{end}.voxsize);
        in=in(ic); % upscale using unique output
        fibsval{side}(unique(fibsin(in,4)),roi)=1;
        ea_dispercent(roi/length(patselection));
    end
    ea_dispercent(1,'end');
end