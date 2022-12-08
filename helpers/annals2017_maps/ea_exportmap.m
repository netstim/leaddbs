function ea_exportmap(n,M,fis,regressor,output,mask,sk,corrtype)

n.img(:)=M;
n.fname=output;

n.dt(1) = 16;
ea_write_nii(n);
if ~exist('mask','var')
    mask=nan;
end

if exist('sk','var')
    switch sk
        case 'k'
            dok=1; dos=0;
        case 's'
            dos=1; dok=0;
        case 'sk'
            dos=1; dok=0;
        otherwise
            dos=0; dok=0;
    end
else
    dos=0;
    dok=0;
end

if isnan(mask)
    mask=1:numel(n.img);
    if dok
        warning('If using k option should apply a mask');
    end
end

if dok
    nanix=(mask);
    nanix(isnan(n.img(nanix)))=0;
    n.img(nanix)=ea_normal(n.img(nanix));
    ea_write_nii(n);
end

if dos
    matlabbatch{1}.spm.spatial.smooth.data = {output};
    matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run',{matlabbatch});
    clear matlabbatch
end
