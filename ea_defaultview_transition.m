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


xyzmm_diff = (togglestates.xyzmm - togglestates_init.xyzmm);
xyztransparencies_diff = (togglestates.xyztransparencies - togglestates_init.xyztransparencies);
v_az_diff = (v.az - v_init.az);
v_el_diff = (v.el - v_init.el);
v_camva_diff = (v.camva - v_init.camva);
v_camup_diff = (v.camup - v_init.camup);
v_out.camproj = 'orthographic';
v_camtarget_diff = (v.camtarget - v_init.camtarget);
v_campos_diff = (v.campos - v_init.campos);

speed_factor = 60;
steps = abs(v_camva_diff) / 40;
steps = steps + max(abs(v_campos_diff)) / 800;
steps = steps + max(abs(v_camtarget_diff)) / 500;

steps = round(steps * speed_factor);

for i = 1:steps
    %togglestates_init.xyzmm = togglestates_init.xyzmm + xyzmm_diff;
    %togglestates_init.xyztransparencies = togglestates_init.xyztransparencies + xyztransparencies_diff;
    %ea_anatomyslices(resultfig,togglestates_init,struct,[]);
    v_out.az = v_init.az + v_az_diff / steps * i;
    v_out.el = v_init.el + v_el_diff / steps * i;
    v_out.camva = v_init.camva + v_camva_diff / steps * i;
    v_out.camup = v_init.camup + v_camup_diff / steps * i;
    v_out.camtarget = v_init.camtarget + v_camtarget_diff / steps * i;
    v_out.campos = v_init.campos + v_campos_diff / steps * i;
    ea_view(v_out);
    drawnow
end
