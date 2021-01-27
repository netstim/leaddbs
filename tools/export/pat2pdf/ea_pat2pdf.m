function  ea_pat2pdf(uipatdir,handles)

options=ea_handles2options(handles);
options.native=0; % MNI space only

[options.root,options.patientname]=fileparts(uipatdir);
options.root=[options.root,filesep];
options=ea_resolve_elspec(options);
options.prefs=ea_prefs(options.patientname);
options.d2 = ea_tdhandles2options([], options.d2);
options.d3.expdf=1;
options.d3.verbose='off';

if ~ea_checkslicespresent(options)
    cuts=ea_writeplanes(options);
end

ea_elvis(options);

fname = [uipatdir, filesep, 'Lead-DBS_Electrode_Localization'];
movefile([fname, '.u3d'], [uipatdir, filesep, 'export', filesep, 'pdf']);
movefile([fname, '.pdf'], [uipatdir, filesep, 'export', filesep, 'pdf']);


function present=ea_checkslicespresent(options)
present=1;
directory=[options.root,options.patientname,filesep];
for el=1:length(options.elspec.contactnames)
    for tracor=1:3
        switch tracor
            case 1
                ppend='_axial.png';
            case 2
                ppend='_coronal.png';
            case 3
                ppend='_sagittal.png';
        end
        if ~exist([directory,options.elspec.contactnames{el},ppend],'file');
            present=0;
        end
    end
end
