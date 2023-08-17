function ea_exportvatmapping(M,options,handles)

selectedregressor=M.clinical.vars{M.ui.clinicallist};
selectedregressor=selectedregressor(M.ui.listselect,:);
if size(selectedregressor,2)==1
    bh='global';
elseif size(selectedregressor,2)==2
    bh='lateralized';
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
        [X,XR,nii]=ea_getXvat(M,options);
        
        save([options.root,options.patientname,filesep,'statvat_results',filesep,'X_',bh,'.mat'],'X','XR','nii','-v7.3');
    else
        load([options.root,options.patientname,filesep,'statvat_results',filesep,'X_',bh,'.mat']);
    end
    
    [T,Mn,Medn,N]=ea_getMvat(M,{X,XR});
    
    % write out results
    ea_writeMvat(M,nii,options,{T,Mn,Medn,N},{'T','mean','median','N'}); 
end


function ea_exportstatvatfiles(M,options,handles)

disp('Need to export VTAs in proper format, this may take a while');

ea_mkdir([options.root,options.patientname,filesep,'statvat_results']);
copyfile([ea_space(options),'bb.nii'],[options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii']);
nii=ea_load_nii([options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii']);
nii.dt(1) = 16;
nii.img(:)=nan;
ea_write_nii(nii);
allV{1}=[options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii'];

cnt=2;

modelLabel = ea_simModel2Label(M.vatmodel);

for pt=1:length(M.patient.list)
    %process left and right separately, check if file exists
    
    is_left_present=false;
    is_right_present=false;

    patStimDir = [M.patient.list{pt},filesep,'stimulations',filesep,ea_nt(options),'gs_',M.guid];

    fname_l = ea_regexpdir(patStimDir, ['_sim-binary_model-',modelLabel,'_hemi-L\.nii$']);
    fname_r = ea_regexpdir(patStimDir, ['_sim-binary_model-',modelLabel,'_hemi-R\.nii$']);

    if ~isempty(fname_l)
        %if file exist, process the right side
        is_left_present=true;
        fname_l = fname_l{1};
    end
    if ~isempty(fname_r)
        %if file exist, process the right side
        is_right_present=true;
        fname_r = fname_r{1};
    end
    
    fns_isprocessed=[false false false];%lh, rh, rh_flipped
    
    %process left side
    if is_left_present
        Vvatl=ea_load_nii(fname_l);
        Vvatl.img(Vvatl.img==0)=nan;
        % writeout
        Vvatl.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','lh','.nii'];
        Vvatl.dt=[16,2];
        spm_write_vol(Vvatl,Vvatl.img);
        Vleft{pt}=Vvatl.fname;
    
        allV{cnt}=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','lh','.nii'];
        
        cnt=cnt+1;
        
        fns_isprocessed(1)=true;
    end

    % process right side
    if is_right_present
        Vvatr=ea_load_nii(fname_r);
        Vvatr.img(Vvatr.img==0)=nan;
        % writeout
        Vvatr.fname=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'];
        Vvatr.dt=[16,2];
        spm_write_vol(Vvatr,Vvatr.img);
        Vright{pt}=Vvatr.fname;
        
        ea_flip_lr_nonlinear( ...
            [options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'], ...
            [options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh_flipped','.nii'], ...
            0);
        
        allV{cnt}=[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_','rh','.nii'];
        
        cnt=cnt+1;
        
        fns_isprocessed(2)=true;% rh
        fns_isprocessed(3)=true;% rh_flipped
    end
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
nii.dt(1) = 2;
delete([options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_bb.nii']);
ea_write_nii(nii);
for pt=1:length(M.patient.list)
    fns={[options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_lh.nii'],...
        [options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_rh.nii'],...
        [options.root,options.patientname,filesep,'statvat_results',filesep,'s',num2str(pt),'_rh_flipped.nii']};

    %prune based on which side was actually processed
    fns=fns(fns_isprocessed);
    
    for f=1:length(fns)
        fn=fns{f};
        ea_conformspaceto([options.root,options.patientname,filesep,'statvat_results',filesep,'statvat_bb.nii'],...
            fn,0);
        nii=ea_load_nii(fn); delete(fn); nii.dt(1) = 2; ea_write_nii(nii);
    end
end
delete([options.root,options.patientname,filesep,'statvat_results',filesep,'bb_nan.nii']);
