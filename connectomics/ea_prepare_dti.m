function redo=ea_prepare_dti(options)
% general DTI preprocessing goes here

directory=[options.root,options.patientname,filesep];

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

% restore dti
if exist([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],'file')
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],...
        [options.root,options.patientname,filesep,options.prefs.dti]);
    redo=1; % apparently prior run crashed - redo to be safe.
end
% restore b0
if exist([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],'file')
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],...
        [options.root,options.patientname,filesep,options.prefs.b0]);
    redo=1; % apparently prior run crashed - redo to be safe.
end
