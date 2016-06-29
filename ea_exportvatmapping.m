function ea_exportvatmapping(M,options,handles)

selectedregressor=M.clinical.vars{get(handles.clinicallist,'Value')};
zselectedregressor=zscore(selectedregressor);
    mkdir([options.root,options.patientname,filesep,'statvat_results']);
allV{1}=[options.earoot,'templates',filesep,'bb.nii'];
zallV{1}=[options.earoot,'templates',filesep,'bb.nii'];
tempzallV{1}=[options.earoot,'templates',filesep,'bb.nii'];

if size(selectedregressor,2)==1;
    bihemispheric=0;
elseif size(selectedregressor,2)==2
    bihemispheric=1;
else
    ea_error('Please select a regressor with entries for each hemisphere or each patient to perform this action.');
end
cnt=2;
for pt=1:length(M.patient.list)
    pdate=0;
    stims=dir([M.patient.list{pt},filesep,'stimulations']);
    for stim=1:length(stims)
        if ~strcmp(stims(stim).name(1),'.')
            if stims(stim).datenum>pdate


                pdate=stims(stim).datenum;

                mostrecentstim=stim;
            end
        end
    end

    try
    Vvatr=ea_load_nii([M.patient.list{pt},filesep,'stimulations',filesep,stims(mostrecentstim).name,filesep,'vat_right.nii']);
    Vvatl=ea_load_nii([M.patient.list{pt},filesep,'stimulations',filesep,stims(mostrecentstim).name,filesep,'vat_left.nii']);
    catch

       keyboard

    end
    Vvatr.img(Vvatr.img==0)=nan;     Vvatl.img(Vvatl.img==0)=nan;
    zVvatr=Vvatr; zVvatl=Vvatl;
    if bihemispheric
        Vvatr.img=Vvatr.img*selectedregressor(pt,1);
        Vvatl.img=Vvatl.img*selectedregressor(pt,2);
        zVvatr.img=zVvatr.img*zselectedregressor(pt,1);
        zVvatl.img=zVvatl.img*zselectedregressor(pt,2);
    else
        Vvatr.img=Vvatr.img*selectedregressor(pt,1);
        Vvatl.img=Vvatl.img*selectedregressor(pt,1);
        zVvatr.img=zVvatr.img*zselectedregressor(pt,1);
        zVvatl.img=zVvatl.img*zselectedregressor(pt,1);
    end

    % writeout

    Vvatr.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'];
    Vvatl.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','lh','.nii'];
    zVvatr.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','rh','.nii'];
    zVvatl.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','lh','.nii'];

    Vvatr.dt=[16,2]; Vvatl.dt=[16,2];      zVvatr.dt=[16,2]; zVvatl.dt=[16,2];
    spm_write_vol(Vvatr,Vvatr.img);     spm_write_vol(Vvatl,Vvatl.img);
    spm_write_vol(zVvatr,zVvatr.img);     spm_write_vol(zVvatl,zVvatl.img);
    Vright{pt}=Vvatr.fname;
    Vleft{pt}=Vvatl.fname;
    zVright{pt}=zVvatr.fname;
    zVleft{pt}=zVvatl.fname;
    ea_flip_lr([options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'],[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh_flipped','.nii']);
    ea_flip_lr([options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','rh','.nii'],[options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','rh_flipped','.nii']);
    allV{cnt}=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','lh','.nii'];
    allV{cnt+1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh_flipped','.nii'];
    zallV{cnt}=[options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','lh','.nii'];
    zallV{cnt+1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','rh_flipped','.nii'];
    tempzallV{cnt}=[options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','lh','_temp.nii'];
    tempzallV{cnt+1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'zs',num2str(pt),'_','rh_flipped','_temp.nii'];
    cnt=cnt+2;
end


% export mean
for z=0:1
    switch z
        case 0
            matlabbatch{1}.spm.util.imcalc.input = allV;
            matlabbatch{1}.spm.util.imcalc.output = 'statvat_mean.nii';
        case 1
            matlabbatch{1}.spm.util.imcalc.input = zallV;
            matlabbatch{1}.spm.util.imcalc.output = 'zstatvat_mean.nii';
    end
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname,filesep,'statvat_results',filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = 'ea_nanmean(X)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 16;

    spm_jobman('run',{matlabbatch});
            ea_crop_nii([options.root,options.patientname,filesep,'statvat_results',filesep,matlabbatch{1}.spm.util.imcalc.output]);

    clear matlabbatch



end


allV{1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_mean.nii'];
zallV{1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_mean.nii'];
tempzallV{1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_mean.nii'];


for pt=1:length(tempzallV)-1

    matlabbatch{1}.spm.util.imcalc.input = zallV;
    matlabbatch{1}.spm.util.imcalc.output = tempzallV{pt+1};
    matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname,filesep,'statvat_results',filesep]};
    matlabbatch{1}.spm.util.imcalc.expression = ['i',num2str(pt+1)];
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

end

dostats=0;
if dostats


% spm factorial:
try
    rmdir([options.root,options.patientname,filesep,'statvat_results',filesep,'SPM'],'s');
end
mkdir([options.root,options.patientname,filesep,'statvat_results',filesep,'SPM']);
matlabbatch{1}.spm.stats.factorial_design.dir = {'/PA/Neuro/_projects/LEAD_pilot/lg/statvat_results/spm'};
%%
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = tempzallV(2:end)';
%%
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    matlabbatch{1}.spm.stats.fmri_est.spmmat = {[options.root,options.patientname,filesep,'statvat_results',filesep,'SPM',filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    matlabbatch{1}.spm.stats.con.spmmat = {[options.root,options.patientname,filesep,'statvat_results',filesep,'SPM',filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'mainFX';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    matlabbatch{1}.spm.stats.results.spmmat = {[options.root,options.patientname,filesep,'statvat_results',filesep,'SPM',filesep,'SPM.mat']};
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = false;
matlabbatch{1}.spm.stats.results.write.tspm.basename = 'spm_result.nii';
    spm_jobman('run',{matlabbatch});
    clear matlabbatch

    % cleanup
    delete([options.root,options.patientname,filesep,'statvat_results',filesep,'*_temp.nii']);
end

