function redo=ea_prepare_dti(options)
% general DTI preprocessing goes here

directory=[options.root,options.patientname,filesep];

if exist([directory,'ea_ftupsampling_applied.mat'],'file')
    ft_upsampling_applied=load([directory,'ea_ftupsampling_applied.mat']);
else % legacy
    ft_upsampling_applied.usfactor=1;
end

if ~exist([options.root,options.patientname,filesep,options.prefs.b0],'file')
    disp('Building DTI files...');
    try %unring
        dti=ea_load_untouch_nii([options.root,options.patientname,filesep,options.prefs.dti]);
        dti.img=ea_unring(dti.img);
        ea_save_untouch_nii(dti,[options.root,options.patientname,filesep,options.prefs.dti]);
    end

    % export B0
    if ~exist([options.root,options.patientname,filesep,options.prefs.b0],'file')
        ea_exportb0(options);
    end

    % export FA
    if ~exist([options.root,options.patientname,filesep,options.prefs.fa],'file')
        ea_isolate_fa(options);
    end
else
    disp('B0 found, no need to rebuild.');
end

usfactor=resolve_usfactor(options.lc.struc.ft.upsample);
if ~isequal(usfactor,ft_upsampling_applied.usfactor)
    redo=1;
else
    redo=0;
end
ft_upsampling_applied.usfactor=usfactor;
save([directory,'ea_ftupsampling_applied.mat'],'-struct','ft_upsampling_applied');


if options.lc.struc.ft.upsample>1
    % stash dti / b0 files:
    if exist([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],'file')
        movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],...
            [options.root,options.patientname,filesep,options.prefs.dti]);
        redo=1; % apparently prior run crashed ? redo to be safe.
    end
    if exist([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],'file')
        movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],...
            [options.root,options.patientname,filesep,options.prefs.b0]);
        redo=1; % apparently prior run crashed ? redo to be safe.
    end
    copyfile([options.root,options.patientname,filesep,options.prefs.dti],...
        [options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)]);
    copyfile([options.root,options.patientname,filesep,options.prefs.b0],...
        [options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)]);
    % upsample:
    V=ea_open_vol([options.root,options.patientname,filesep,options.prefs.dti]);
    
    ea_reslice_nii([options.root,options.patientname,filesep,options.prefs.dti],...
        [options.root,options.patientname,filesep,options.prefs.dti],V.voxsize./usfactor,0,0,1,[],[],1);
    ea_reslice_nii([options.root,options.patientname,filesep,options.prefs.b0],...
        [options.root,options.patientname,filesep,options.prefs.b0],V.voxsize./usfactor,0,0,1,[],[],1);

    %% add methods dump:
    cits={
        'Dyrby, T. B., Lundell, H., Burke, M. W., Reislev, N. L., Paulson, O. B., Ptito, M., & Siebner, H. R. (2013). Interpolation of diffusion weighted imaging datasets. NeuroImage, 103(C), 1?12. http://doi.org/10.1016/j.neuroimage.2014.09.005'
        };
    ea_methods(options,['Raw diffusion data was upsampled using bspline-interpolation with a factor of ',num2str(usfactor),' following the concept described in Dyrby et al. 2014.'],cits);

end


function factor=resolve_usfactor(popupmenuentry)
switch popupmenuentry
    case 1
        factor=1;
    case 2
        factor=1.5;
    case 3
        factor=2;
    case 4
        factor=3;
end