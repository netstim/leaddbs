function ea_addshoot

root=[ea_getearoot,'templates',filesep,'dartel',filesep];
for dt=1:6
    nii=ea_load_nii([root,'dartelmni_',num2str(dt),'.nii']);
    if size(nii.img,4)>3
        disp(['File: ',[root,'dartelmni_',num2str(dt),'.nii'],' already has 4 entries in 4th dim.']);
        return
    else
    matlabbatch{1}.spm.util.split.vol = {[root,'dartelmni_',num2str(dt),'.nii,1']};
    matlabbatch{1}.spm.util.split.outdir = {root};
    spm_jobman('run',{matlabbatch});
    clear matlabbatch
    

        X=-sum(nii.img,4);
        X=X+1;
        
        nii=ea_load_nii([root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',1),'.nii']);
        nii.img=X;
        nii.fname=[root,'/dartelmni_',num2str(dt),'_',sprintf('%05.0f',4),'.nii'];
        ea_write_nii(nii);
        
        for i=1:4
            matlabbatch{1}.spm.util.cat.vols{i} = [root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',i),'.nii'];
        end
        matlabbatch{1}.spm.util.cat.vols=matlabbatch{1}.spm.util.cat.vols';
        matlabbatch{1}.spm.util.cat.name = ['dartelmni_',num2str(dt),'.nii'];
        matlabbatch{1}.spm.util.cat.dtype = 0;
        spm_jobman('run',{matlabbatch});
        clear matlabbatch
    end
    for i=1:4 % cleanup
        delete([root,'dartelmni_',num2str(dt),'_',sprintf('%05.0f',i),'.nii']);
    end
    delete([root,'dartelmni_',num2str(dt),'.mat']);
end