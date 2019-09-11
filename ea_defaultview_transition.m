function [] = ea_defaultview_transition(v,togglestates)
% transition between current view and defaultview

H = findall(0,'type','figure');
resultfig = H(~cellfun(@isempty,strfind({H(:).Name},{'Electrode-Scene'})));
resultfig = resultfig(1); % take the first if there are many.

togglestates_init = getappdata(resultfig,'togglestates');
togglestates_init.xyztransparencies(~togglestates_init.xyztoggles) = 0;
togglestates_init.xyztoggles = [1 1 1]; 
togglestates_init.refreshview = 1;

v_init = ea_view();

steps = 50;

xyzmm_diff = (togglestates.xyzmm - togglestates_init.xyzmm) / steps;
xyztransparencies_diff = (togglestates.xyztransparencies - togglestates_init.xyztransparencies) / steps;
v_az_diff = (v.az - v_init.az) / steps;
v_el_diff = (v.el - v_init.el) / steps;
v_camva_diff = (v.camva - v_init.camva) / steps;
v_camup_diff = (v.camup - v_init.camup) / steps;
v_out.camproj = 'orthographic';
v_camtarget_diff = (v.camtarget - v_init.camtarget) / steps;
v_campos_diff = (v.campos - v_init.campos) / steps;


for i = 1:steps
    %togglestates_init.xyzmm = togglestates_init.xyzmm + xyzmm_diff;
    %togglestates_init.xyztransparencies = togglestates_init.xyztransparencies + xyztransparencies_diff;
    %ea_anatomyslices(resultfig,togglestates_init,struct,[]);
    v_out.az = v_init.az + v_az_diff * i;
    v_out.el = v_init.el + v_el_diff * i;
    v_out.camva = v_init.camva + v_camva_diff * i;
    v_out.camup = v_init.camup + v_camup_diff * i;
    v_out.camtarget = v_init.camtarget + v_camtarget_diff * i;
    v_out.campos = v_init.campos + v_campos_diff * i;
    ea_view(v_out);
    drawnow
end
