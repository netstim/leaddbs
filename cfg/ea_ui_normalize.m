function ea_ui_normalize(P)
% Function to receive cfg_util jobs.
options=ea_defaultoptions;
[options.root,patientname]=fileparts(fileparts(P.foldername{1}));
options.root=[options.root,filesep];
options.prefs=ea_prefs(patientname);
% load reconstruction results


switch P.method
    case 1 % Schönecker, linear threestep
        ea_normalize_schoenecker(options);
    case 2 % Schönecker, linear threestep, incl. preop data
        ea_normalize_schoenecker_pre(options);
    case 3 % Witt, nonlinear.
        ea_normalize_witt(options);
end