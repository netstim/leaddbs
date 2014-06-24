function ea_ui_manualheight(P)
% Function to receive cfg_util jobs.
options=ea_defaultoptions;
[options.root,patientname]=fileparts(fileparts(P.foldername{1}));
options.root=[options.root,filesep];
options.prefs=ea_prefs(patientname);
% load reconstruction results
try
    load([P.foldername{1},'ea_reconstruction']);
catch
    error('No reconstruction information found. Please run reconstruction first.');
end
mcfig=figure('name',[patientname,': Manual Height Correction'],'numbertitle','off');

[mcfig,coords_mm,trajectory]=ea_manualheightcorrection(mcfig,coords_mm,trajectory,patientname,options);
close(mcfig)

% save results
ea_write_fiducials(coords_mm,fullfile(options.root,patientname,'ea_coords.fcsv'));
save([options.root,patientname,filesep,'ea_reconstruction'],'trajectory','coords_mm');
