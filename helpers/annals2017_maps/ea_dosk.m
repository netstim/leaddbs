function fina=ea_dosk(fina,mask)
[pth,fn,ext]=fileparts(fina);
nii=ea_load_nii(fina);
nii.img(~mask)=nan;
nii.img(mask)=ea_normal(nii.img(mask));
nii.fname=fullfile(pth,['k',fn,ext]);
ea_write_nii(nii);
matlabbatch{1}.spm.spatial.smooth.data = {[fullfile(pth,['k',fn]),'.nii,1']};
matlabbatch{1}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});
clear matlabbatch
delete(fullfile(pth,['k',fn,ext]));
if ~isBIDSFileName(fina)
    fina=fullfile(pth,['sk',fn,ext]);
else
    parsed = parseBIDSFilePath(fina);
    if ~isfield(parsed, 'desc')
        parsed.desc = '';
    end
    fina = setBIDSEntity(fina, 'desc', [parsed.desc,'NormSmooth']);
    movefile(fullfile(pth,['sk',fn,ext]), fina);
end
