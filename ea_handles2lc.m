function lc=ea_handles2lc(handles)

% General settings
parcellations = get(handles.parcellation, 'String');
lc.general.parcellation = parcellations{get(handles.parcellation,'Value')};


% Graph options:
lc.graph.struc_func_sim=get(handles.struc_func_sim,'Value');
lc.graph.nodal_efficiency=get(handles.nodal_efficiency,'Value');
lc.graph.eigenvector_centrality=get(handles.eigenvector_centrality,'Value');
lc.graph.degree_centrality=get(handles.degree_centrality,'Value');
lc.graph.fthresh=str2double(get(handles.fthresh,'String'));
lc.graph.sthresh=str2double(get(handles.sthresh,'String'));


% functional options:
lc.func.compute_CM=get(handles.compute_CM_func,'Value');
lc.func.compute_GM=get(handles.compute_GM_func,'Value');
lc.func.prefs.TR=str2double(get(handles.TR,'String'));


% structural options:
lc.struc.compute_CM=get(handles.compute_CM_struc,'Value');
lc.struc.compute_GM=get(handles.compute_GM_struc,'Value');
ftFunctions = getappdata(handles.leadfigure, 'ftFunctions');
lc.struc.ft.upsample.factor=get(handles.upsamplingfactor,'Value');
lc.struc.ft.upsample.how=get(handles.use_internal_upsampling,'Value');
lc.struc.ft.method = ftFunctions{get(handles.ftmethod,'Value')};
lc.struc.ft.do=get(handles.perf_ft,'Value');
lc.struc.ft.normalize=get(handles.normalize_fibers,'Value');

if strcmp(get(handles.fiber_count,'Visible'), 'on')
    lc.struc.ft.dsistudio.fiber_count = str2double(get(handles.fiber_count,'String'));
end
