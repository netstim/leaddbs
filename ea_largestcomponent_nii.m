function ea_largestcomponent_nii(fname,cnum)
% helper function that reduces the nifti components to the largest
% component.

nii=ea_load_nii(fname);
nii.img(isnan(nii.img))=0;
nii.img=logical(nii.img);
C=bwconncomp(nii.img);


ls=cellfun(@length,C.PixelIdxList,'Uniformoutput',0);
ls=cell2mat(ls);
[~,ix]=sort(ls);

[pth,fn,ext]=fileparts(fname);

assignhemispheres=0;
if ischar(cnum)
    if strcmp(cnum,'hemispheres')
        cnum=2;
        assignhemispheres=1;
    end
end

for cc=0:cnum-1
    exp{cc+1}=nii;
    exp{cc+1}.img(:)=0;
    exp{cc+1}.img(C.PixelIdxList{ix(end-cc)})=1;
    if assignhemispheres
        
        [xx,yy,zz]=ind2sub(size(exp{cc+1}.img),C.PixelIdxList{ix(end-cc)});
        XYZ=[xx,yy,zz,ones(size(xx,1),1)]';
        XYZ=exp{cc+1}.mat*XYZ;
        xx=XYZ(1,:);
        mxx(cc+1)=mean(xx);
        % do not export as of yet (see section below).
    else
        exp{cc+1}.fname=fullfile(pth,[fn,'_c',num2str(cc),ext]);
        spm_write_vol(exp{cc+1},exp{cc+1}.img);
    end
end

% export if in hemispheres mode
if assignhemispheres
    if mxx(1)<mxx(2)
        append={'lh','rh'};
    else
        append={'rh','lh'};
    end
    for cc=1:2
        exp{cc}.fname=fullfile(pth,[fn,'_',append{cc},ext]);
        spm_write_vol(exp{cc},exp{cc}.img);
    end
end
