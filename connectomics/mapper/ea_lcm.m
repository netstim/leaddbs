function ea_lcm(options)


% the order may be crucial since VAT-seeds are resliced in fMRI
if options.lcm.struc.do
    
    % main wrapper for lead connectome mapper runs
    if strcmp(options.lcm.seeddef,'vats')
        options.lcm.seeds=ea_resolvevatseeds(options,'dMRI');
        if isempty(options.lcm.odir)
            options.lcm.odir=[fileparts(options.lcm.seeds{1}),filesep];
        end
    end
    
    ea_lcm_struc(options);
end
if options.lcm.func.do
    
    if strcmp(options.lcm.seeddef,'vats')
        options.lcm.seeds=ea_resolvevatseeds(options,'fMRI');
        if isempty(options.lcm.odir)
            options.lcm.odir=[fileparts(options.lcm.seeds{1}),filesep];
        end
    end
    
    ea_lcm_func(options);
end



function seeds=ea_resolvevatseeds(options,modality)
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

% prepare for dMRI
switch modality
    case 'dMRI'
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
        Cnii.fname=[vatdir,'vat_seed_compound_dMRI.nii'];
        ea_write_nii(Cnii);
        ea_crop_nii(Cnii.fname);
        delete([vatdir,'tmp_*']);
        seeds{1}=[vatdir,'vat_seed_compound_dMRI.nii'];
        disp('Done.');
        
    case 'fMRI'
        % prepare for fMRI
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
                cname=options.lcm.func.connectome;
                if ismember('>',cname)
                    delim=strfind(cname,'>');
                    subset=cname(delim+2:end);
                    cname=cname(1:delim-2);
                end
                d=load([ea_getconnectomebase,'fMRI',filesep,cname,filesep,'dataset_info.mat']);
                d.dataset.vol.space.fname=[vatdir,'tmp_space.nii'];
                ea_write_nii(d.dataset.vol.space);
                ea_conformspaceto(d.dataset.vol.space.fname,[vatdir,'tmp_',sidec,'.nii'],1);
                nii(cnt)=ea_load_nii([vatdir,'tmp_',sidec,'.nii']);
                cnt=cnt+1;
            end
        end
        Cnii=nii(1);
        for n=2:length(nii)
            Cnii.img=Cnii.img+nii(n).img;
        end
        Cnii.fname=[vatdir,'vat_seed_compound_fMRI.nii'];
        ea_write_nii(Cnii);
        delete([vatdir,'tmp_*']);
        seeds{1}=[vatdir,'vat_seed_compound_fMRI.nii'];
        disp('Done.');
end
