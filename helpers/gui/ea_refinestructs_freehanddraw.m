function [lineobj,xs,ys] = ea_freehanddraw(ax,handles,drawvoxz,varargin)
% [LINEOBJ,XS,YS] = FREEHANDDRAW(ax_handle,handles,line_options)
%
% Draw a smooth freehand line object on the current axis (default),
% or on the axis specified by handle in the first input argument.
% Left-click to begin drawing, right-click to terminate, or double-click
% to close contour and terminate.
% 
%
% INPUT ARGUMENTS:  First:      axis handle (optional)
%                  Additional: valid line property/value pairs
%
% OUTPUT ARGUMENTS: 1) Handle to line object
%                  2) x-data
%                  3) y-data
% (Note that output args 2 and 3 can also be extracted from the first output
% argument.)
%
% Ex: [myobj,xs,ys] = freehanddraw(gca,'color','r','linewidth',3);
%     freehanddraw('linestyle','--');
%
% Written by Brett Shoelson, PhD
% shoelson@helix.nih.gov
% 3/29/05
% Modified:
% 10/05/05: Now allows double-click closing of contour
trafun=get(handles.tra,'ButtonDownFcn');
set(handles.tra,'ButtonDownFcn', '');
corfun=get(handles.cor,'ButtonDownFcn');
set(handles.cor,'ButtonDownFcn', '');
sagfun=get(handles.sag,'ButtonDownFcn');
set(handles.sag,'ButtonDownFcn', '');

set(handles.discardfiducial,'Visible','off');
set(handles.approvefiducial,'Visible','off');
ea_csremovedrawings(handles);	
%Get current figure and axis parameters
oldvals = get(gcf);
oldhold = ishold(gca);

hold on;
set(handles.refinestatus,'String','Click on a point to start and begin drawing.');

set(handles.checkstructures,'doublebuffer','on');

%Get the initial point
[xs,ys,button] = ginput(1);
set(handles.refinestatus,'String','Great! Now right-click to terminate drawing.');

%Create and store line object
	lineobj = line(xs,ys,'tag','tmpregsel',varargin{1:end});

setappdata(handles.checkstructures,'lineobj',lineobj);

%Modify wbmf of current figure to update lineobject on mouse motion
set(handles.checkstructures,'windowbuttonmotionfcn',@wbmfcn);
set(handles.checkstructures,'windowbuttondownfcn',@wbdfcn);
%Wait for right-click or double-click
while ~strcmp(get(handles.checkstructures,'SelectionType'),'alt') & ~strcmp(get(handles.checkstructures,'SelectionType'),'open')
	drawnow;
end

% Everything from here on is already "wrap up" code, i.e. when doubleclick or
% rightclick was issued.

%Extract xyz data from line object for return in output variables
%(Also retrievable from first output argument)
	xs = get(getappdata(handles.checkstructures,'lineobj'),'xdata')';

	ys = get(getappdata(handles.checkstructures,'lineobj'),'ydata')';

% get z coordinate:

linefiducial=zeros(length(xs),3);

[planedim,onedim,secdim]=ea_getdims(drawvoxz(3),1);
linefiducial(:,planedim)=drawvoxz(1); % fill z coord
linefiducial(:,onedim)=xs;
linefiducial(:,secdim)=ys;
setappdata(handles.checkstructures,'linefiducial',linefiducial);
setappdata(handles.checkstructures,'fiducialview',drawvoxz(3));
set(handles.refinestatus,'String','');

set(handles.discardfiducial,'Visible','on');
set(handles.approvefiducial,'Visible','on');
%Clear temporary variables from base workspace
evalin('caller','clear tmpx tmpy tmpz done gca lineobj');

%Reset figure parameters
set(handles.checkstructures,'Pointer',oldvals.Pointer,...
	'windowbuttonmotionfcn',oldvals.WindowButtonMotionFcn,...
    'windowbuttondownfcn',oldvals.WindowButtonDownFcn);
if isfield(oldvals, 'DoubleBuffer')
    set(handles.checkstructures,'doublebuffer',oldvals.DoubleBuffer);
end
%Reset hold value of the axis
if ~oldhold, hold off; end 

ea_setatlascline(handles);
setappdata(handles.checkstructures,'bdfuns',{trafun,corfun,sagfun});
% set(handles.tra,'ButtonDownFcn', trafun);
% set(handles.cor,'ButtonDownFcn', corfun);
% set(handles.sag,'ButtonDownFcn', sagfun);

function wbmfcn(varargin)
lineobj = getappdata(gcf,'lineobj');
if strcmp(get(gcf,'selectiontype'),'normal');
    tmpx = get(lineobj,'xdata');
    tmpy = get(lineobj,'ydata');
    a=get(gca,'currentpoint');
    set(lineobj,'xdata',[tmpx,a(1,1)],'ydata',[tmpy,a(1,2)]);
    drawnow;
else
    setappdata(gcf,'lineobj',lineobj);
end

function wbdfcn(varargin)
lineobj = getappdata(gcf,'lineobj');
if strcmp(get(gcf,'selectiontype'),'open')
    tmpx = get(lineobj,'xdata');
    tmpy = get(lineobj,'ydata');
    a=get(gca,'currentpoint');
    set(lineobj,'xdata',[tmpx,tmpx(1)],'ydata',[tmpy,tmpy(1)]);
    setappdata(gcf,'lineobj',lineobj);
    drawnow;
end
return