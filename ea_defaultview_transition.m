function [] = ea_defaultview_transition(v,togglestates)
% transition between current view and defaultview

H = findall(0,'type','figure');
resultfig = H(~cellfun(@isempty,strfind({H(:).Name},{'Electrode-Scene'})));
resultfig = resultfig(1); % take the first if there are many.

togglestates_init = getappdata(resultfig,'togglestates');
togglestates_init.xyztransparencies(~togglestates_init.xyztoggles) = 0;
togglestates_init.xyztoggles(logical(togglestates_init.xyztoggles)) = 1;
togglestates_init.refreshview = 1;

steps = 50;

xyzmm_diff = (togglestates.xyzmm - togglestates_init.xyzmm) / steps;
xyztransparencies_diff = (togglestates.xyztransparencies - togglestates_init.xyztransparencies) / steps;

campos_diff = (v.campos - campos) / steps;
camtarget_diff = (v.camtarget - camtarget) / steps;


for i = 1:steps
    togglestates_init.xyzmm = togglestates_init.xyzmm + xyzmm_diff;
    togglestates_init.xyztransparencies = togglestates_init.xyztransparencies + xyztransparencies_diff;
    ea_anatomyslices(resultfig,togglestates_init,struct,[]);
    camdolly(camtarget_diff(1),camtarget_diff(2),camtarget_diff(3),'movetarget','data')
    camdolly(campos_diff(1)-camtarget_diff(1),campos_diff(2)-camtarget_diff(2),campos_diff(3)-camtarget_diff(3),'fixtarget','data')
    drawnow
end
