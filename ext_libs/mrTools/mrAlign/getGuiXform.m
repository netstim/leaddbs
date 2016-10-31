function xform = getGuiXform(handles)
% function xform = getGuiXform(handles)
%
% Returns 4x4 xform deteremined by translation and rotation text fields
%
% djh 1/5/2007

global ALIGN

% Get trans from GUI
trans = [str2double(get(handles.transX,'String')),...
    str2double(get(handles.transY,'String')),...
    str2double(get(handles.transZ,'String'))];
% str2double returns nan for non-numeric strings, change to zeros
nanIndices = find(isnan(trans));
trans(nanIndices) = zeros(size(nanIndices));
		
% Get rot from GUI
rot = [str2double(get(handles.rotX,'String')),...
    str2double(get(handles.rotY,'String')),...
    str2double(get(handles.rotZ,'String'))];
% str2double returns nan for non-numeric strings, change to zeros
nanIndices = find(isnan(rot));
rot(nanIndices) = zeros(size(nanIndices));
		
% Convert from deg to radians
rot = pi/180 * rot;
		
% Stuff rot and trans into 4x4 matrix
rotXform = eye(4);
cosx = cos(rot(1));		sinx = sin(rot(1));
cosy = cos(rot(2));		siny = sin(rot(2));
cosz = cos(rot(3));		sinz = sin(rot(3));
rotXform(1,1) = cosz*cosy+sinz*sinx*siny;
rotXform(1,2) = sinz*cosy-cosz*sinx*siny;
rotXform(1,3) = cosx*siny;
rotXform(2,1) = -sinz*cosx;
rotXform(2,2) = cosz*cosx;
rotXform(2,3) = sinx;
rotXform(3,1) = sinz*sinx*cosy-cosz*siny;
rotXform(3,2) = -cosz*sinx*cosy-sinz*siny;
rotXform(3,3) = cosx*cosy;

inpCenter = ones(4,1);
inpCenter(1:3) = ALIGN.inplaneSize/2;
volCenter = ALIGN.xform * inpCenter;

prerotXform = eye(4);
prerotXform(:,4) = -volCenter;

transXform = eye(4);
transXform(1,4) = trans(1);
transXform(2,4) = trans(2);
transXform(3,4) = trans(3);

xform = transXform * inv(prerotXform) * rotXform * prerotXform;

end
