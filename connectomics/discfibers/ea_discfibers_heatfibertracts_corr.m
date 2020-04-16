function [fibsval_sum,fibsval_mean,fibsval_peak,fibsval_5peak]=ea_discfibers_heatfibertracts_corr(obj,fibcell,XYZmm,niivx,valsmm)
% function extracts fibers from a connectome connected to ROIs in the
% roilist and assigns them correlative values based on vals. Vals needs to be of
% same length as roilist, assigning a value for each ROI.

disp('ROI fiber analysis');

patInd = 1:2*length(obj.allpatients);
vizz=0;
dthresh=2*mean(niivx);
for side=1:2
    fibsval_sum{side}=zeros(size(fibcell{side},1),length(patInd),'single'); % 5th column will add up values, 6th will take note how many entries were summed.
    fibsval_mean{side}=zeros(size(fibcell{side},1),length(patInd),'single'); % 5th column will add up values, 6th will take note how many entries were summed.
    fibsval_peak{side}=zeros(size(fibcell{side},1),length(patInd),'single'); % 5th column will add up values, 6th will take note how many entries were summed.
    fibsval_5peak{side}=zeros(size(fibcell{side},1),length(patInd),'single'); % 5th column will add up values, 6th will take note how many entries were summed.

    ea_dispercent(0,['Iterating ROI, side ',num2str(side)]);   
    for roi=1:length(patInd)
        tree=KDTreeSearcher(XYZmm{roi,side}(1:2:end,1:3)); % light downsample
        if vizz
           figure
           hold on
           plot3(XYZmm{roi,side}(:,1),XYZmm{roi,side}(:,2),XYZmm{roi,side}(:,3),'r.');
           for fib=1:100
           plot3(fibcell{side}{fib}(:,1),fibcell{side}{fib}(:,2),fibcell{side}{fib}(:,3),'r.');
           end
        end
        for fib=1:length(fibcell{side})
            [IX,D]=knnsearch(tree,fibcell{side}{fib}(1:2:end,:),'Distance','chebychev'); % light downsample
            in=D<dthresh;
            if any(in)
                tv=valsmm{roi,side}(IX(in));
                fibsval_sum{side}(fib,roi)=sum(tv);
                fibsval_mean{side}(fib,roi)=mean(tv);
                fibsval_peak{side}(fib,roi)=max(tv);
                tv=sort(tv,'descend');
                fibsval_5peak{side}(fib,roi)=mean(tv(1:ceil(0.05*length(tv))));
            end
        end
        ea_dispercent(roi/length(patInd));
    end
    ea_dispercent(1,'end');
end
