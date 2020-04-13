function [fibsin,XYZmm,nii,valsmm]=ea_discfibers_genroilist_connfibers(fibers, roilist, D)

prefs = ea_prefs;
thresh = prefs.machine.vatsettings.horn_ethresh;
allroilist=cat(2,roilist{:});
% tree=KDTreeSearcher(fibers(:,1:3));
% load in all ROI

    vizz=0;

if isempty(D)
    cnt=1;
    XYZmm = cell(size(allroilist,1), 1);
    nii = cell(size(allroilist,1), 2);
    valsmm = cell(size(allroilist,1), 1);
    ea_dispercent(0,'Aggregating ROI');
    for roi=1:size(allroilist,1)
        nii{roi,1}=ea_load_nii(allroilist{roi,1});
        ixs=find(nii{roi,1}.img(:)>(thresh*100));
        [xx,yy,zz]=ind2sub(size(nii{roi,1}.img),ixs);
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XY=nii{roi,1}.mat*XYZvx;
        XYZmm{roi,1}=XY(1:3,:)';
        valsmm{roi,1}=nii{roi,1}.img(ixs);
        
        
        nii{roi,2}=ea_load_nii(allroilist{roi,2});
        ixs=find(nii{roi,2}.img(:)>(thresh*100));
        [xx,yy,zz]=ind2sub(size(nii{roi,2}.img),ixs);
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XY=nii{roi,2}.mat*XYZvx;
        XYZmm{roi,2}=XY(1:3,:)';
        valsmm{roi,2}=nii{roi,2}.img(ixs);
        
        if vizz
            figure
            plot3(XYZmm{roi,1}(:,1),XYZmm{roi,1}(:,2),XYZmm{roi,1}(:,3),'r.')
        end

        
        ea_dispercent(roi/size(allroilist,1));
    end
    ea_dispercent(1,'end');
    AllXYZ{1}=cell2mat(XYZmm(:,1));
    AllXYZ{2}=cell2mat(XYZmm(:,1)); 
    
else % aggregate AllXYZ from prior data   
    AllXYZ{1}=cell2mat(D.XYZmm(:,1));
    AllXYZ{2}=cell2mat(D.XYZmm(:,1));    
end
dthreshs=[3,0.7,mean(nii{end}.voxsize)];
dsfactor=[0.5,2,inf];
for side=1:2
    ea_dispt(['Selecting connected fibers, side ',num2str(side)]);
    fibsin{side}=fibers(:,1:4);
    cnt=1;
    for downsample=dsfactor
        disp(['Pass #',num2str(cnt)]);
        % isolate fibers that are connected to any ROI:
        if isinf(downsample) % last pass, no downsampling
            ea_dispt('Final pass, no downsampling');
            dsAllXYZ{side}=AllXYZ{side};
        else
            ea_dispt('Downsampling spatially');
            dsAllXYZ{side}=unique(round(AllXYZ{side}*downsample)/downsample,'rows');
        end
        if vizz
            figure
            plot3(dsAllXYZ{side}(:,1),dsAllXYZ{side}(:,2),dsAllXYZ{side}(:,3),'r.');
        end    
        ea_dispt('KNN Search');
        [~,D]=knnsearch(dsAllXYZ{side}(1:4:end,1:3),fibsin{side}(:,1:3),'Distance','chebychev');
        in=D<(dthreshs(cnt)); % larger threshold here since seed vals have been rounded above. Also, this is just to reduce size of the connectome for speed improvements, so we can be more liberal here.
        fibsin{side}=fibsin{side}(ismember(fibsin{side}(:,4),unique(fibsin{side}(in,4))),:);
        if vizz
            figure
            hold on
            plot3(AllXYZ{side}(:,1),AllXYZ{side}(:,2),AllXYZ{side}(:,3),'r.');
            plot3(fibsin{side}(:,1),fibsin{side}(:,2),fibsin{side}(:,3),'b.');
        end
        cnt=cnt+1;
    end
    ea_dispt('');
end
