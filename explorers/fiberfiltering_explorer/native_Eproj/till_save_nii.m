function till_save_nii(vol,img)
file = vol.fname;
if contains(file,'.nii.gz')
    for i=1:numel(vol)
        vol(i).fname = file(1:end-3);
    end
    if ndims(img)<4
        spm_write_vol(vol,img);
    elseif ndims(img) == 4
        nimg = size(img,4);
        if ~isequal(vol(:).mat)
            warning('4D nifti has different .mats!')
            keyboard
        end
        myfiles={};
        for i = 1:nimg
            voltmp = vol(i);
            voltmp.n = [1,1];
            imgtmp = squeeze(img(:,:,:,i));
            voltmp.fname = ['./temp' num2str(i) '.nii'];
            spm_write_vol(voltmp,imgtmp);
            myfiles = vertcat(myfiles,{voltmp.fname});
        end
        spm_file_merge(myfiles,vol(1).fname);
        for i = 1:nimg
            delete(myfiles{i})
        end
    end
    gzip(file(1:end-3))
    delete(file(1:end-3))
elseif ~contains(file,'.nii.gz') && contains(file,'.nii')
    if ndims(img)<4
        spm_write_vol(vol,img);
    elseif ndims(img) == 4
        nimg = size(img,4);
        if ~isequal(vol(:).mat)
            warning('4D nifti has different .mats!')
            keyboard
        end
        myfiles={};
        for i = 1:nimg
            voltmp = vol(i);
            voltmp.n = [1,1];
            imgtmp = squeeze(img(:,:,:,i));
            voltmp.fname = ['./temp' num2str(i) '.nii'];
            spm_write_vol(voltmp,imgtmp);
            myfiles = vertcat(myfiles,{voltmp.fname});
        end
        spm_file_merge(myfiles',vol(1).fname);
        for i = 1:nimg
            delete(myfiles{i})
        end
    end
else
    error('File seems to be neither .nii nor .nii.gz!')
end
end