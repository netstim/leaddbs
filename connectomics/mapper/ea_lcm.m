function ea_lcm(options)


% the order may be crucial since VAT-seeds are resliced in fMRI
if options.lcm.struc.do
    
    % main wrapper for lead connectome mapper runs
    if strcmp(options.lcm.seeddef,'vats')
        originalseeds=options.lcm.seeds;
        options.lcm.seeds=ea_resolvevatseeds(options,'dMRI');
        if isempty(options.lcm.odir)
            options.lcm.odir=[fileparts(options.lcm.seeds{1}),filesep];
        end
    end
    
    ea_lcm_struc(options);
end
if options.lcm.func.do
    
    if strcmp(options.lcm.seeddef,'vats')
        if exist('originalseeds','var')
            options.lcm.seeds=originalseeds;
        end
        options.lcm.seeds=ea_resolvevatseeds(options,'fMRI');
        if isempty(options.lcm.odir)
            options.lcm.odir=[fileparts(options.lcm.seeds{1}),filesep];
        end
    end
    
    ea_lcm_func(options);
end

% convert VTA seeds also if neither func or struc conn is chosen.
if (~options.lcm.func.do) && (~options.lcm.struc.do)
    if strcmp(options.lcm.seeddef,'vats')
        ea_resolvevatseeds(options,'dMRI');
        ea_resolvevatseeds(options,'fMRI');
    end
end



function seeds=ea_resolvevatseeds(options,modality)
disp('Preparing VATs as seedfiles...');
vatdir=[options.root,options.patientname,filesep,'stimulations',filesep,options.lcm.seeds,filesep];

suffices={'binary','efield','efield_gauss'};

[~,dowhich]=ismember(options.prefs.lcm.vatseed,suffices);
for suffix=dowhich
    
    switch suffices{suffix}
        case 'binary'
            addstr='';
            dinterp=0;
        case 'efield_gauss'
            addstr='_efield_gauss';
            dinterp=1;
        case 'efield'
            addstr='_efield';
            dinterp=1;
    end
    if strcmp(options.prefs.lcm.vatseed,suffices{suffix})
        keepthisone=1;
    else
        keepthisone=0;
    end
    
    % prepare for dMRI
    switch modality
        case 'dMRI'
            %if ~exist([vatdir,'vat_seed_compound_dMRI',addstr,'.nii'],'file');
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
                        ea_conformspaceto([ea_space,'bb.nii'],[vatdir,'tmp_',sidec,'.nii'],dinterp);
                        nii(cnt)=ea_load_nii([vatdir,'tmp_',sidec,'.nii']);
                        cnt=cnt+1;
                    end
                end
                Cnii=nii(1);
                for n=2:length(nii)
                    Cnii.img=Cnii.img+nii(n).img;
                end
                Cnii.fname=[vatdir,'vat_seed_compound_dMRI',addstr,'.nii'];
                ea_write_nii(Cnii);
                ea_crop_nii(Cnii.fname);
                delete([vatdir,'tmp_*']);
                
                ea_split_nii_lr(Cnii.fname);
                disp('Done.');
            %end
            if keepthisone
                seeds{1}=[vatdir,'vat_seed_compound_dMRI',addstr,'.nii'];
            end
        case 'fMRI'
            % prepare for fMRI
            %if ~exist([vatdir,'vat_seed_compound_fMRI',addstr,'.nii'],'file');
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
                        tnii=ea_load_nii([vatdir,'tmp_',sidec,'.nii']);
                        tnii.dt=[16,0];
                        ea_write_nii(tnii);
                        cname=options.lcm.func.connectome;
                        if ismember('>',cname)
                            delim=strfind(cname,'>');
                            subset=cname(delim+2:end);
                            cname=cname(1:delim-2);
                        end
                        d=load([ea_getconnectomebase,'fMRI',filesep,cname,filesep,'dataset_info.mat']);
                        d.dataset.vol.space.fname=[vatdir,'tmp_space.nii'];
                        d.dataset.vol.space.dt=[16,0];
                        ea_write_nii(d.dataset.vol.space);
                        ea_conformspaceto(d.dataset.vol.space.fname,[vatdir,'tmp_',sidec,'.nii'],1);
                        
                        nii(cnt)=ea_load_nii([vatdir,'tmp_',sidec,'.nii']);
                        nii(cnt).img(isnan(nii(cnt).img))=0;
                        if ~any(nii(cnt).img(:))
                            msgbox(['Created empty VTA for ',options.patientname,', ',sidec,' hemisphere.']);
                        end
                        cnt=cnt+1;
                    end
                    
                end
                Cnii=nii(1);
                for n=2:length(nii)
                    Cnii.img=Cnii.img+nii(n).img;
                end
                Cnii.fname=[vatdir,'vat_seed_compound_fMRI',addstr,'.nii'];
                ea_write_nii(Cnii);
                delete([vatdir,'tmp_*']);
             
                ea_split_nii_lr(Cnii.fname);
                disp('Done.');
            %end
            if keepthisone
                seeds{1}=[vatdir,'vat_seed_compound_fMRI',addstr,'.nii'];
            end
    end
end
