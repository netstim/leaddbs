function ea_show_light(resultfig,viz)
% if viz = 1; turn on
% if viz = 0; turn off
% default viz=1;
if nargin<2
    viz=1;
end
        set(0,'CurrentFigure',resultfig); 
cam_lamp=camlight('headlight'); % not modifiable, infinite light.
%set(cam_lamp,'Color',[1,1,1]);
ceiling_lamp=light('Position',[0 0 100]); % not modifiable, infinite light.

%set(ceiling_lamp,'Color',[0.5,0.5,0.5]);

right_lamp=light('Position',[-100 0 0]); % not modifiable, infinite light.
%set(right_lamp,'Color',[0.5,0.3,0.5]);

left_lamp=light('Position',[100 0 0]); % not modifiable, infinite light.
%set(left_lamp,'Color',[0.5,0.5,0.3]);

%lightobj=light('Position',[30 30 30],'Style','local');
%lightbulb=plot3(30, 30, 30,'o','MarkerSize',20,'MarkerFaceColor','y','MarkerEdgeColor','k');
%setappdata(resultfig,'lightbulb',lightbulb);
%setappdata(resultfig,'lightobj',lightobj);

if ~viz
    set(cam_lamp,'Visible','off')
    set(ceiling_lamp,'Visible','off')
    set(right_lamp,'Visible','off')
    set(left_lamp,'Visible','off')
end

setappdata(resultfig,'right_lamp',right_lamp);
setappdata(resultfig,'left_lamp',left_lamp);
setappdata(resultfig,'ceiling_lamp',ceiling_lamp);
setappdata(resultfig,'cam_lamp',cam_lamp);

%ea_draggable(lightbulb,@movelight);


function movelight(lightbulb)
xx=get(lightbulb,'xdata');
yy=get(lightbulb,'ydata');
zz=get(lightbulb,'zdata');
lightobj=getappdata(gcf,'lightobj');
set(lightobj,'Position',[xx yy zz]);
