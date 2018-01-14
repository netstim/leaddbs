function ea_exportvatmapping(M,options,handles)

selectedregressor=M.clinical.vars{M.ui.clinicallist};
selectedregressor=selectedregressor(M.ui.listselect,:);
if size(selectedregressor,2)==1
    bh='lateralized';
elseif size(selectedregressor,2)==2
    bh='global';
else
    ea_error('Please select a regressor with entries for each hemisphere or each patient to perform this action.');
end

if ~exist([options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_bb.nii'],'file')
    ea_exportstatvatfiles(M,options,handles)
end

% generate unique hash-ID based on patient selection
hshid=ea_datahash(M.ui.listselect);

if ~exist([options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_',M.clinical.labels{M.ui.clinicallist},'_mean_',hshid,'.nii'],'file') % file already generated
    if ~exist([options.root,options.patientname,filesep,'statvat_results',filesep,'X_',bh,'.mat'],'file')
        [X,XR]=ea_getXvat(M,options);
        keyboard
        save([options.root,options.patientname,filesep,'statvat_results',filesep,'X_',bh,'.mat'],'X','XR');
    else
        load([options.root,options.patientname,filesep,'statvat_results',filesep,'X_',bh,'.mat']);
    end
    
    [T,Mn,Medn,N,nii]=ea_getMvat(M,{X,XR});
    % write out results
    ea_writeMvat(M,nii,options,{T,Mn,Medn,N},{'T','mean','median','N'}); 
end








function ea_exportstatvatfiles(M,options,handles)

disp('Need to export VTAs in proper format, this may take a while');

mkdir([options.root,options.patientname,filesep,'statvat_results']);
copyfile([ea_space(options),'bb.nii'],[options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii']);
nii=ea_load_nii([options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii']);
nii.dt=[16,1];
nii.img(:)=nan;
ea_write_nii(nii);
allV{1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii'];


cnt=2;
for pt=1:length(M.patient.list)
    if M.ui.detached % process locally in lead group directory
        Vvatr=ea_load_nii([options.root,options.patientname,filesep,M.patient.list{pt},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_right.nii']);
        Vvatl=ea_load_nii([options.root,options.patientname,filesep,M.patient.list{pt},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_left.nii']);
    else
        Vvatr=ea_load_nii([M.patient.list{pt},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_right.nii']);
        Vvatl=ea_load_nii([M.patient.list{pt},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat_left.nii']);
    end
    Vvatr.img(Vvatr.img==0)=nan;     Vvatl.img(Vvatl.img==0)=nan;

    % writeout
    Vvatr.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'];
    Vvatl.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','lh','.nii'];
    
    Vvatr.dt=[16,2]; Vvatl.dt=[16,2];  
    spm_write_vol(Vvatr,Vvatr.img);     spm_write_vol(Vvatl,Vvatl.img);
    Vright{pt}=Vvatr.fname;
    Vleft{pt}=Vvatl.fname;
    ea_flip_lr_nonlinear([options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'],[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh_flipped','.nii'],0);
    allV{cnt}=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','lh','.nii'];
    allV{cnt+1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'];
    cnt=cnt+2;
end


% export mean to get bounding box
matlabbatch{1}.spm.util.imcalc.input = allV';
matlabbatch{1}.spm.util.imcalc.output = 'statvat_bb.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {[options.root,options.patientname,filesep,'statvat_results',filesep]};
matlabbatch{1}.spm.util.imcalc.expression = 'ea_nanmean(X)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 16;

spm_jobman('run',{matlabbatch});
ea_crop_nii([options.root,options.patientname,filesep,'statvat_results',filesep,matlabbatch{1}.spm.util.imcalc.output],'','nn');
clear matlabbatch

% no conform each VTA to bb
nii=ea_load_nii([options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_bb.nii']);
nii.dt=[2,0];
delete([options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_bb.nii']);
ea_write_nii(nii);
for pt=1:length(M.patient.list)
    fns={[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_lh.nii'],...
        [options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_rh.nii'],...
        [options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_rh_flipped.nii']};
    for f=1:length(fns)
        fn=fns{f};
        ea_conformspaceto([options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_bb.nii'],...
            fn,0);
        nii=ea_load_nii(fn); delete(fn); nii.dt=[2,0]; ea_write_nii(nii);
    end
end
delete([options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii']);