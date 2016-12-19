function crop_tool(fi,region)

dtd = dtdstruct_read([fi '_DTD.mat']);

dtd.eigenVec_struc.dataAy = dtd.eigenVec_struc.dataAy(region(1,1):region(1,2),region(2,1):region(2,2),region(3,1):region(3,2),:,:);
dtd.eigenVal_struc.dataAy = dtd.eigenVal_struc.dataAy(region(1,1):region(1,2),region(2,1):region(2,2),region(3,1):region(3,2),:,:);
dtd.b0_image_struc.dataAy = dtd.b0_image_struc.dataAy(region(1,1):region(1,2),region(2,1):region(2,2),region(3,1):region(3,2),:,:);
if isfield(dtd,'meanDWI_image_struc'),
   dtd.meanDWI_image_struc.dataAy = dtd.meanDWI_image_struc.dataAy(region(1,1):region(1,2),region(2,1):region(2,2),region(3,1):region(3,2),:,:);
end;

dtdstruct_write(dtd,[fi '_crop_DTD.mat']);

if exist([fi '_HARDI.mat'],'file')
    mr = mrstruct_read([fi '_HARDI.mat']);
    mr.dataAy = mr.dataAy(region(1,1):region(1,2),region(2,1):region(2,2),region(3,1):region(3,2),:,:);
    mrstruct_write(mr,[fi '_crop_HARDI.mat']);
end
sizeAy = size(dtd.b0_image_struc.dataAy);

if exist([fi '_slice1.mat'],'file')
    for m = region(3,1) : region(3,2)
        mr = mrstruct_read([fi '_slice',num2str(m),'.mat']);
        mr.dataAy = mr.dataAy(region(1,1):region(1,2),region(2,1):region(2,2),:);
        mrstruct_write(mr,[fi '_crop_slice',num2str(m),'.mat']);
    end
end

if exist([fi '_WM.mat'],'file')
    mr = mrstruct_read([fi '_WM.mat']);
    mr.dataAy = mr.dataAy(region(1,1):region(1,2),region(2,1):region(2,2),region(3,1):region(3,2));
    mrstruct_write(mr,[fi '_crop_WM.mat']);
end

[file path] = uigetfile('*.mat','Select the ROI to crop!');
if isempty(file)
    return
else
    mask = maskstruct_read(fullfile(path,file));
    no_masks = size(mask.maskCell,1);
    for m = 1 : no_masks
        mask.maskCell{m,1} = mask.maskCell{m,1}(region(1,1):region(1,2),region(2,1):region(2,2),region(3,1):region(3,2));
    end
    mask.sizeAy = sizeAy;
    mask.mrsProp.user.size = sizeAy;
    mask.mrsProp.user.matrix = [sizeAy(1) sizeAy(2)];
    mask.mrsProp.user.SliceNo = sizeAy(3);
    maskstruct_write(mask,fullfile(path,[file(1:end-4),'_crop.mat']));
end
