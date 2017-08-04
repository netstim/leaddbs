function ea_deselectall_dicom(handles)
if get(handles.dicomcheck,'Value') || get(handles.assignnii,'Value')
    
    % DBS:
    try set(handles.coregct_checkbox,'Enable','off'); end
    try set(handles.coregctcheck,'Enable','off'); end
    try set(handles.normalize_checkbox,'Enable','off'); end
    try set(handles.normcheck,'Enable','off'); end
    try set(handles.doreconstruction_checkbox,'Enable','off'); end
    try set(handles.manualheight_checkbox,'Enable','off'); end
    try set(handles.include_lead_connectome_subroutine,'Enable','off'); end
    try set(handles.writeout2d_checkbox,'Enable','off'); end
    try set(handles.render_checkbox,'Enable','off'); end
    try set(handles.scrf,'Enable','off'); end

    % Connectome:
    try set(handles.perf_ft,'Enable','off'); end
    try set(handles.normalize_fibers,'Enable','off'); end
    try set(handles.compute_CM_struc,'Enable','off'); end
    try set(handles.compute_GM_struc,'Enable','off'); end
    try set(handles.compute_CM_func,'Enable','off'); end
    try set(handles.compute_GM_func,'Enable','off'); end
    try set(handles.degree_centrality,'Enable','off'); end
    try set(handles.struc_func_sim,'Enable','off'); end
    try set(handles.eigenvector_centrality,'Enable','off'); end
    try set(handles.nodal_efficiency,'Enable','off'); end


else
    % DBS:
    try set(handles.coregct_checkbox,'Enable','on'); end
    try set(handles.coregctcheck,'Enable','on'); end
    try set(handles.normalize_checkbox,'Enable','on'); end
    try set(handles.normcheck,'Enable','on'); end
    try set(handles.doreconstruction_checkbox,'Enable','on'); end
    try set(handles.manualheight_checkbox,'Enable','on'); end
    try set(handles.include_lead_connectome_subroutine,'Enable','on'); end
    try set(handles.writeout2d_checkbox,'Enable','on'); end
    try set(handles.render_checkbox,'Enable','on'); end
    try set(handles.scrf,'Enable','on'); end
    try    ea_switchctmr(handles,get(handles.MRCT,'Value')); end

    
    % Connectome:
    try set(handles.perf_ft,'Enable','on'); end
    try set(handles.normalize_fibers,'Enable','on'); end
    try set(handles.compute_CM_struc,'Enable','on'); end
    try set(handles.compute_GM_struc,'Enable','on'); end
    try set(handles.compute_CM_func,'Enable','on'); end
    try set(handles.compute_GM_func,'Enable','on'); end
    try set(handles.degree_centrality,'Enable','on'); end
    try set(handles.struc_func_sim,'Enable','on'); end
    try set(handles.eigenvector_centrality,'Enable','on'); end
    try set(handles.nodal_efficiency,'Enable','on'); end
    
    
end

