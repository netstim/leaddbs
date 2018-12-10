function ea_combinenii_generic(fis,mask,ofname)
usenormal=1;
denoise=1;
if exist('mask','var')
    if ~isempty(mask)
        msk=ea_load_nii(mask);
        maskidx=(msk.img(:)>0);
    end
end

[pth,fname,ext]=fileparts(fis{1});
for fi=1:length(fis)
    if denoise
        uuid=ea_generate_uuid;
        if ~exist(fullfile(pth,'tmp'),'dir')
            mkdir(fullfile(pth,'tmp'));
        end
        ea_delete(fullfile(pth,'tmp',['tmp',uuid,'.nii']));
        ea_denoise_mri(fis{fi},fullfile(pth,'tmp',['tmp',uuid,'.nii']));
        try
            nii=ea_load_nii(fullfile(pth,'tmp',['tmp',uuid,'.nii']));
            ea_delete(fullfile(pth,'tmp',['tmp',uuid,'.nii']));
        catch
            nii=ea_load_nii(fis{fi});
        end
    else
        nii=ea_load_nii(fis{fi});
    end
    if ~exist('maskidx','var')
        maskidx=[1:length(nii.img(:))]';
    end
    nii.img(~maskidx)=0;
    if usenormal
    nii.img(maskidx)=ea_normal(nii.img(maskidx),1,'TRUE');
    else
        nii.img(maskidx)=zscore(nii.img(maskidx));
    end
    if ~exist('cnii','var')
        cnii=nii;
    else
        cnii.img=cnii.img.*nii.img;
    end
    
end


if exist('ofname','var')
    cnii.fname=ofname;
else
    cnii.fname=fullfile(pth,[fname,'_combined','.nii']);
end
if usenormal
    cnii.img(maskidx)=ea_normal(cnii.img(maskidx),1,'TRUE');
else
    cnii.img(maskidx)=zscore(cnii.img(maskidx));
end
cnii.img=cnii.img.*100;
cnii.img(cnii.img>255)=255;
cnii.img(cnii.img<-255)=-255;
ea_write_nii(cnii);



