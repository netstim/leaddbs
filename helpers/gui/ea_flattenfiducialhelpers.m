function ea_flattenfiducialhelpers(~,~,handles)

uipatdir=getappdata(handles.leadfigure,'uipatdir');
cnt=1;
for pt=1:length(uipatdir)
    fids=dir([uipatdir{pt},filesep,'fiducials',filesep,'*.nii.gz']);
    
    for fi=1:length(fids)
        % pt folder
        gunzip([uipatdir{pt},filesep,'fiducials',filesep,fids(fi).name]);
        delete([uipatdir{pt},filesep,'fiducials',filesep,fids(fi).name]);
        fidcell{fi}=[uipatdir{pt},filesep,'fiducials',filesep,fids(fi).name(1:end-3)];
        
        % templates
        if exist([ea_space,'fiducials',filesep,fids(fi).name],'file')
            gunzip([ea_space,'fiducials',filesep,fids(fi).name]);
            delete([ea_space,'fiducials',filesep,fids(fi).name]);
            tfidcell{fi}=[ea_space,'fiducials',filesep,fids(fi).name(1:end-3)];
        end
    end
    
    % flatten ROI:
    
    fguid=ea_generate_uuid;
    if exist(tfidcell{1},'file') % in case multiple pts selected could be that template fiducials have already been flattened (if they were the same ones)
        matlabbatch{1}.spm.util.imcalc.input = tfidcell';
        matlabbatch{1}.spm.util.imcalc.output = [fguid,'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[ea_space,'fiducials']};
        matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',{matlabbatch});
        smoothgzip([ea_space,'fiducials'],[fguid,'.nii']);
        ea_delete(tfidcell);
    end
    for pt=1:length(uipatdir)
        matlabbatch{1}.spm.util.imcalc.input = fidcell';
        matlabbatch{1}.spm.util.imcalc.output = [fguid,'.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {[uipatdir{pt},filesep,'fiducials']};
        matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',{matlabbatch});
        smoothgzip([uipatdir{pt},filesep,'fiducials'],[fguid,'.nii']);
        ea_delete(fidcell);
    end
end


function smoothgzip(pathn,filen)
matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pathn,filen)};
matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{1}.spm.spatial.smooth.dtype = 8;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});
movefile(fullfile(pathn,['s',filen]),fullfile(pathn,filen));
gzip(fullfile(pathn,filen));
delete(fullfile(pathn,filen));


