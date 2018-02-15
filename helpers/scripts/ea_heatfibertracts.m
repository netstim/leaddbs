function fibsin=ea_heatfibertracts(cfile,roilist,vals,thresh)
% function extracts fibers from a connectome connected to ROIs in the
% roilist and assigns them correlative values based on vals. Vals needs to be of
% same length as roilist, assigning a value for each ROI.


disp('ROI fiber analysis');

[fibers,idx,voxmm,mat]=ea_loadfibertracts(cfile);

if ~iscell(roilist)
    roilist={roilist};
    vals={vals};
end


allroilist=cat(2,roilist{:});
%tree=KDTreeSearcher(fibers(:,1:3));
    % load in all ROI
    ea_dispercent(0,'Aggregating ROI');
    for roi=1:length(allroilist);
        nii{roi}=ea_load_nii(allroilist{roi});
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
        ea_dispercent(roi/length(allroilist));
    end
    ea_dispercent(1,'end');
    ea_dispt('Selecting connected fibers');
    % isolate fibers that are connected to any ROI:
    AllXYZ=unique(round(AllXYZ),'rows');
    
    [~,D]=knnsearch(AllXYZ(:,1:3),fibers(:,1:3),'Distance','chebychev');
    
    in=D<(mean(nii{end}.voxsize)*2); % larger threshold here since seed vals have been rounded above. Also, this is just to reduce size of the connectome for speed improvements, so we can be more liberal here.
    
    fibsin=fibers(ismember(fibers(:,4),unique(fibers(in,4))),:);
    fibsval=[zeros(size(fibsin,1),length(allroilist))]; % 5th column will add up values, 6th will take note how many entries were summed.
    %for group=1:length(roilist)
        
        % now color fibsin based on predictive value of improvement
        ea_dispt('');
        
        ea_dispercent(0,'Iterating ROI');
        
        for roi=1:length(allroilist);
            [~,D]=knnsearch(XYZmm{roi}(:,1:3),fibsin(:,1:3),'Distance','chebychev');
            in=D<mean(nii{end}.voxsize);
            fibsval(ismember(fibsin(:,4),unique(fibsin(in,4))),roi)=1;
            ea_dispercent(roi/length(allroilist));
            
        end
        ea_dispercent(1,'end');
        ea_dispt('Correlating fibers with values');
        cnt=1;
        for group=1:length(roilist)
            thisgroupidx=cnt:(cnt+length(roilist{group}))-1;
            cnt=cnt+length(roilist{group});
            %R=corr(fibsval',vals,'type','spearman');
            repvals=repmat(vals{group}',size(fibsval,1),1);
            try
            nfibsval=fibsval(:,thisgroupidx); nfibsval(nfibsval==0)=nan; posvals=repvals.*nfibsval;
            nfibsval=double(~fibsval(:,thisgroupidx)); nfibsval(nfibsval==0)=nan; negvals=repvals.*nfibsval;
            catch
                keyboard
            end
            [h,p,ci,stats]=ttest2(posvals',negvals');
            fibsin=[fibsin,stats.tstat'];
        end
            fibsin(all(isnan(fibsin(:,5:end)),2),:)=[];
ea_dispt('Exporting fiberset');

ea_dispt('');
