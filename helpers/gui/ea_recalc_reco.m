function ea_recalc_reco(~,~,handles)


uipatdirs=getappdata(handles.leadfigure,'uipatdir');
options.native=1;
options.hybridsave=1;
options.sides=1:2;

for pt=1:length(uipatdirs)
    options=ea_getptopts(uipatdirs{pt},options);
    disp(['Re-propagating reconstruction from native (postop) -> native (preop) -> template space: ', options.patientname]); 
    [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options);
end