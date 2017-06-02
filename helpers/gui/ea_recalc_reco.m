function [coords_mm,trajectory,markers]=ea_recalc_reco(~,~,handles)

if ischar(handles)
    uipatdirs={handles};
else
    uipatdirs=getappdata(handles.leadfigure,'uipatdir');
end
options.native=1;
options.hybridsave=1;
options.sides=1:2;

for pt=1:length(uipatdirs)
    options=ea_getptopts(uipatdirs{pt},options);
    disp(['Re-propagating reconstruction from native (postop) -> native (preop) -> template space: ', options.patientname]); 
    [~,~,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
    options.elmodel=elmodel;
    options=ea_resolve_elspec(options);
    [coords_mm,trajectory,markers]=ea_resolvecoords(markers,options);
    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options);
    if nargout % export coords in MNI space.
        options.native=0;
        [coords_mm,trajectory,markers]=ea_load_reconstruction(options);
    end
end