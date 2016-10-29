function ea_lcm(options)
% main wrapper for lead connectome mapper runs
if strcmp(options.lcm.seeddef,'vats')
    options.lcm.seeds=ea_resolvevatseeds(options);
    if isempty(options.lcm.odir)
        options.lcm.odir=[fileparts(options.lcm.seeds{1}),filesep];
    end
end

% the order may be crucial since VAT-seeds are resliced in fMRI
if options.lcm.struc.do
    ea_lcm_struc(options);
end
if options.lcm.func.do
    ea_lcm_func(options);
end



function seeds=ea_resolvevatseeds(options)
disp('Preparing VATs as seedfiles...');
vatdir=[options.root,options.patientname,filesep,'stimulations',filesep,options.lcm.seeds,filesep];
switch options.prefs.lcm.vatseed
    case 'binary'
        addstr='';
    case 'efield_gauss'
        addstr='_efield_gauss';
    case 'efield'
        addstr='_efield';
end
cnt=1;
for side=1:2
    switch side
        case 1
            sidec='right';
        case 2
            sidec='left';
    end
    
    if exist([vatdir,'vat',addstr,'_',sidec,'.nii'],'file')
        copyfile([vatdir,'vat',addstr,'_',sidec,'.nii'],[vatdir,'tmp_',sidec,'.nii']);
        ea_conformspaceto([ea_getearoot,'templates',filesep,'mni_hires_bb.nii'],[vatdir,'tmp_',sidec,'.nii'],1);
        nii(cnt)=ea_load_nii([vatdir,'tmp_',sidec,'.nii']);
        cnt=cnt+1;
    end
end
Cnii=nii(1);
for n=2:length(nii)
    Cnii.img=Cnii.img+nii(n).img;
end
Cnii.fname=[vatdir,'vat_seed_compound.nii'];
ea_write_nii(Cnii);
ea_crop_nii(Cnii.fname);
delete([vatdir,'tmp_*']);
seeds{1}=[vatdir,'vat_seed_compound.nii'];
disp('Done.');
