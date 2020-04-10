function [fibsval_sum,fibsval_mean,fibsval_peak,fibsval_5peak]=ea_discfibers_heatfibertracts_corr(obj,fibcell,XYZmm,niivx,valsmm)
% function extracts fibers from a connectome connected to ROIs in the
% roilist and assigns them correlative values based on vals. Vals needs to be of
% same length as roilist, assigning a value for each ROI.

disp('ROI fiber analysis');

patselection=[obj.patientselection,obj.patientselection+length(obj.allpatients)];

dthresh=2*mean(niivx);
for side=1:2
    fibsval_sum{side}=zeros(size(fibcell{side},1),length(patselection),'single'); % 5th column will add up values, 6th will take note how many entries were summed.
    fibsval_mean{side}=zeros(size(fibcell{side},1),length(patselection),'single'); % 5th column will add up values, 6th will take note how many entries were summed.
    fibsval_peak{side}=zeros(size(fibcell{side},1),length(patselection),'single'); % 5th column will add up values, 6th will take note how many entries were summed.
    fibsval_5peak{side}=zeros(size(fibcell{side},1),length(patselection),'single'); % 5th column will add up values, 6th will take note how many entries were summed.

    ea_dispercent(0,['Iterating ROI, side ',num2str(side)]);   
    for roi=1:length(patselection)
        tree=KDTreeSearcher(XYZmm{roi,side}(:,1:3));
        for fib=1:length(fibcell{side})
            [IX,D]=knnsearch(tree,fibcell{side}{fib},'Distance','chebychev');
            in=D<dthresh;
            if any(in)
                tv=valsmm{roi,side}(IX(in));
                fibsval_sum{side}(fib,roi)=sum(tv);
                fibsval_mean{side}(fib,roi)=mean(tv);
                fibsval_peak{side}(fib,roi)=max(tv);
 
%                 switch obj.efieldmetric
%                     case 'sum'
%                         fibsval{side}(fib,roi)=sum(tv);
%                     case 'mean'
%                         fibsval{side}(fib,roi)=mean(tv);
%                     case 'peak'
%                         fibsval{side}(fib,roi)=max(tv);
%                 end
            end
        end
        ea_dispercent(roi/length(patselection));
    end
    ea_dispercent(1,'end');
    
end