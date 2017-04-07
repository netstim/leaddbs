function fibsin=ea_heatfibertracts(cfile,roilist,vals,thresh)
% function extracts fibers from a connectome connected to ROIs in the
% roilist and assigns them correlative values based on vals. Vals needs to be of
% same length as roilist, assigning a value for each ROI.


disp('ROI fiber analysis');

[fibers,idx,voxmm,mat]=ea_loadfibertracts(cfile);


%tree=KDTreeSearcher(fibers(:,1:3));
% load in all ROI
ea_dispercent(0,'Aggregating ROI');
for roi=1:length(roilist);
    nii{roi}=ea_load_nii(roilist{roi});
    if exist('thresh','var')
        nii{roi}.img=double(nii{roi}.img>thresh);
    end
    [xx,yy,zz]=ind2sub(size(nii{roi}.img),find(nii{roi}.img));
    XYZvx=[xx,yy,zz,ones(length(xx),1)]';
    XYZmm{roi}=nii{roi}.mat*XYZvx;
    XYZmm{roi}=XYZmm{roi}(1:3,:)';
    if ~exist('AllXYZ','var')
        AllXYZ=XYZmm{roi};
    else
        AllXYZ=[AllXYZ;XYZmm{roi}];
    end
    ea_dispercent(roi/length(roilist));
end
ea_dispercent(1,'end');
ea_dispt('Combining ROI');
% isolate fibers that are connected to any ROI:
AllXYZ=unique(round(AllXYZ),'rows');

[~,D]=knnsearch(AllXYZ(:,1:3),fibers(:,1:3),'Distance','chebychev');

in=D<(mean(nii{end}.voxsize)*2); % larger threshold here since seed vals have been rounded above. Also, this is just to reduce size of the connectome for speed improvements, so we can be more liberal here.

fibsin=fibers(ismember(fibers(:,4),unique(fibers(in,4))),:);
fibsval=[zeros(size(fibsin,1),length(roilist))]; % 5th column will add up values, 6th will take note how many entries were summed.

% now color fibsin based on predictive value of improvement
ea_dispt('');

ea_dispercent(0,'Iterating ROI');

for roi=1:length(roilist);
    [~,D]=knnsearch(XYZmm{roi}(:,1:3),fibsin(:,1:3),'Distance','chebychev');
    in=D<mean(nii{end}.voxsize);
    fibsval(ismember(fibsin(:,4),unique(fibsin(in,4))),roi)=1;
    ea_dispercent(roi/length(roilist));
    
end
ea_dispercent(1,'end');
ea_dispt('Correlating fibers with values');

%R=corr(fibsval',vals,'type','spearman');
repvals=repmat(vals',size(fibsval,1),1);
nfibsval=fibsval; nfibsval(nfibsval==0)=nan; posvals=repvals.*nfibsval;
nfibsval=double(~fibsval); nfibsval(nfibsval==0)=nan; negvals=repvals.*nfibsval;

[h,p,ci,stats]=ttest2(posvals',negvals');

ea_dispt('Exporting fiberset');
fibsin=[fibsin,stats.tstat'];
fibsin(isnan(fibsin(:,5)),:)=[];

ea_dispt('');



