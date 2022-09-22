function ea_show_light(resultfig,viz)
% if viz = 1; turn on
% if viz = 0; turn off
% default viz=1;
if nargin<2
    viz=1;
end

set(0,'CurrentFigure',resultfig);

prefs = ea_prefs;

CamLight=light('style','infinite','Color',prefs.d3.camlightcolor); % not modifiable, infinite light.
%set(cam_lamp,'Color',[1,1,1]);
camlight(CamLight,'headlight'); % move light object.

CeilingLight=light('Position',[0 0 10],'style','local','Color',prefs.d3.ceilinglightcolor); % not modifiable, infinite light.
%set(ceiling_lamp,'Color',[0.5,0.5,0.5]);

RightLight=light('Position',[-100 0 0],'style','infinite','Color',prefs.d3.rightlightcolor); % not modifiable, infinite light.
%set(right_lamp,'Color',[0.5,0.3,0.5]);

LeftLight=light('Position',[100 0 0],'style','infinite','Color',prefs.d3.leftlightcolor); % not modifiable, infinite light.
%set(left_lamp,'Color',[0.5,0.5,0.3]);

%lightobj=light('Position',[30 30 30],'Style','local');
%lightbulb=plot3(30, 30, 30,'o','MarkerSize',20,'MarkerFaceColor','y','MarkerEdgeColor','k');
%setappdata(resultfig,'lightbulb',lightbulb);
%setappdata(resultfig,'lightobj',lightobj);

if ~viz
    set(CamLight,'Visible','off')
    set(CeilingLight,'Visible','off')
    set(RightLight,'Visible','off')
    set(LeftLight,'Visible','off')
end

setappdata(resultfig,'RightLight',RightLight);
setappdata(resultfig,'LeftLight',LeftLight);
setappdata(resultfig,'CeilingLight',CeilingLight);
setappdata(resultfig,'CamLight',CamLight);

%ea_draggable(lightbulb,@movelight);


function movelight(lightbulb)
xx=get(lightbulb,'xdata');
yy=get(lightbulb,'ydata');
zz=get(lightbulb,'zdata');
lightobj=getappdata(gcf,'lightobj');
set(lightobj,'Position',[xx yy zz]);
