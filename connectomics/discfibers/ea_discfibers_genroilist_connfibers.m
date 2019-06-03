function [fibsin, XYZmm, nii, valsmm]=ea_discfibers_genroilist_connfibers(fibers, roilist, patselection, thresh)

if ~exist('thresh', 'var')
    thresh = 0;
    mode = 'binary';
else
    mode = 'efield';
end

allroilist=cat(2,roilist{:});
% tree=KDTreeSearcher(fibers(:,1:3));
% load in all ROI

cnt=1;
XYZmm = cell(length(allroilist), 1);
nii = cell(length(allroilist), 2);
if strcmp(mode, 'binary')
    valsmm = [];
else
    valsmm = cell(length(allroilist), 1);
end
ea_dispercent(0,'Aggregating ROI');
for roi=1:length(allroilist)
    if size(allroilist,2)==2 % left and right entered separately, combine.
        nii{roi,1}=ea_load_nii(allroilist{roi,1});

        [xx,yy,zz]=ind2sub(size(nii{roi,1}.img),find(nii{roi,1}.img(:)>thresh));
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XY=nii{roi,1}.mat*XYZvx;
        XYZmm{roi}=XY(1:3,:)';
        if strcmp(mode, 'efield')
            valsmm{roi}=nii{roi,1}.img((nii{roi,1}.img(:))>thresh);
        end

        nii{roi,2}=ea_load_nii(allroilist{roi,2});
        [xx,yy,zz]=ind2sub(size(nii{roi,2}.img),find(nii{roi,2}.img(:)>thresh));
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XY=nii{roi,2}.mat*XYZvx;
        XYZmm{roi}=[XYZmm{roi};XY(1:3,:)'];
        if strcmp(mode, 'efield')
            valsmm{roi}=[valsmm{roi};nii{roi,2}.img((nii{roi,2}.img(:))>thresh)];
        end
    else
        %if strcmp(mode, 'efield')
            keyboard % this is not maintained and currently not functioning.
        %end

        nii{roi}=ea_load_nii(allroilist{roi});
        nii{roi}.img=double(nii{roi}.img>thresh);

        keyboard
        [xx,yy,zz]=ind2sub(size(nii{roi}.img),find(nii{roi}.img));
        XYZvx=[xx,yy,zz,ones(length(xx),1)]';
        XYZmm{roi}=nii{roi}.mat*XYZvx;
        XYZmm{roi}=XYZmm{roi}(1:3,:)';
    end

    if ~exist('AllXYZ','var')
        AllXYZ=XYZmm{roi};
    else
        AllXYZ=[AllXYZ;XYZmm{roi}];
    end

    if ismember(roi,patselection)
        XYZmmSel{cnt}=XYZmm{roi};
        niiSel{cnt,:}=nii{roi,:};

%         if ~exist('AllXYZSel','var')
%             AllXYZSel=XYZmmSel{cnt};
%         else
%             AllXYZSel=[AllXYZSel;XYZmm{cnt}];
%         end

        cnt=cnt+1;
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
if size(fibsin,2)>4
   fibsin(:,5:end)=[];
end

% AllXYZ=AllXYZSel; % now that connected fibers were selected, replace with selected VTAs only to avoid confusion.
XYZmm = XYZmmSel;
nii = niiSel;
if ~isempty(valsmm)
    valsmm = valsmm(patselection);
end
