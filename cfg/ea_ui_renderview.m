function ea_ui_renderview(P)
% Function to receive cfg_util jobs.
options=ea_defaultoptions;
[options.root,patientname]=fileparts(fileparts(P.foldername{1}));
options.root=[options.root,filesep];
options.prefs=ea_prefs(patientname);
% load reconstruction results


options.elmodel=P.elspec;
options=ea_resolve_elspec(options);
options.showatlases=P.showatlases;

options.showfibers=P.showfibers;

options.showconnectivities=P.showconnectivities;
options.prolongelectrode=P.prolongelectrode;
options.dostimulation=0;
% Prior Results are loaded here inside the function (this way, function
% can be called just by giving the patient directory.

[coords_mm,resultfig]=ea_render_view(patientname,options);

% save scene as matlab figure
saveas(resultfig,[options.root,patientname,filesep,'eAuto_scene.fig']);
%figure2xhtml([options.root,patientname,filesep,'eAuto_scene'],resultfig);