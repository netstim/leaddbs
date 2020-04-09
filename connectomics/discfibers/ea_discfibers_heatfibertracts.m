function [fibsval]=ea_discfibers_heatfibertracts(obj,fibcell,fibsin,XYZmm,niivx)
% function extracts fibers from a connectome connected to ROIs in the
% roilist and assigns them correlative values based on vals. Vals needs to be of
% same length as roilist, assigning a value for each ROI.

disp('ROI fiber analysis');

patselection=[obj.patientselection,obj.patientselection+length(obj.allpatients)];


    
searchfibsin=round(fibsin(:,1:3).*4)/4; % reduce precision a bit (0.25 mm) to speed up search
[searchfibsin,~,ic]=unique(searchfibsin,'rows');

for side=1:2
    fibsval{side}=zeros(size(fibcell,1),length(patselection),'logical'); % 5th column will add up values, 6th will take note how many entries were summed.
    
    % now color fibsin based on predictive value of improvement
    
    ea_dispercent(0,['Iterating ROI, side ',num2str(side)]);
 
    for roi=1:length(patselection)
        [~,D]=knnsearch(XYZmm{roi,side}(:,1:3),searchfibsin,'Distance','chebychev');
        in=D<mean(niivx);
        in=in(ic); % upscale using unique output
        fibsval{side}(unique(fibsin(in,4)),roi)=1;
        ea_dispercent(roi/length(patselection));
    end
    ea_dispercent(1,'end');
end