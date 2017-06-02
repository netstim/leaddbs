function mertoggles = ea_getsettogglestates(handles,varargin)

% Example type A: 
%   varargin{1}=mertoggles;
%
% Example type B: ea_getsettogglestates(handles,1,mermarkers(Ridx(1)).dat.leaddepth);
%   varargin{1}=side;
%   varargin{2}=dist;
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh, Brain Modulation Lab
%
% Ari Kappel

if size(varargin,2)==1
    mertoggles = varargin{1};
elseif size(varargin,2)==2
    side = varargin{1};
    dist = varargin{2};
end

try
    resultfig=getappdata(handles.mercontrolfig,'resultfig');
catch
    handles = getappdata(handles,'UsedByGUIData_m');
    resultfig=getappdata(handles.mercontrolfig,'resultfig');
end
merhandles = getappdata(resultfig,'merhandles');

% reset states based on gui:
if ~exist('mertoggles','var')

mertoggles.keycontrol=[get(handles.keycentral_left,'Value'),get(handles.keyanterior_left,'Value'),...
    get(handles.keyposterior_left,'Value'),get(handles.keylateral_left,'Value'),get(handles.keymedial_left,'Value');...
    get(handles.keycentral_right,'Value'),get(handles.keyanterior_right,'Value'),get(handles.keyposterior_right,'Value'),...
    get(handles.keylateral_right,'Value'),get(handles.keymedial_right,'Value')];
mertoggles.togglestates=[get(handles.togglecentral_left,'Value'),get(handles.toggleanterior_left,'Value'),...
    get(handles.toggleposterior_left,'Value'),get(handles.togglelateral_left,'Value'),get(handles.togglemedial_left,'Value');...
    get(handles.togglecentral_right,'Value'),get(handles.toggleanterior_right,'Value'),get(handles.toggleposterior_right,'Value'),...
    get(handles.togglelateral_right,'Value'),get(handles.togglemedial_right,'Value')];

else
    
side=1; sidestr='right';
set(handles.(['keycentral_',sidestr]),'Value',mertoggles.keycontrol(side,1))
    set(handles.(['keyanterior_',sidestr]),'Value',mertoggles.keycontrol(side,2))
    set(handles.(['keyposterior_',sidestr]),'Value',mertoggles.keycontrol(side,3))
    set(handles.(['keylateral_',sidestr]),'Value',mertoggles.keycontrol(side,4))
    set(handles.(['keymedial_',sidestr]),'Value',mertoggles.keycontrol(side,5))
set(handles.(['togglecentral_',sidestr]),'Value',mertoggles.togglestates(side,1))
    set(handles.(['toggleanterior_',sidestr]),'Value',mertoggles.togglestates(side,2))
    set(handles.(['toggleposterior_',sidestr]),'Value',mertoggles.togglestates(side,3))
    set(handles.(['togglelateral_',sidestr]),'Value',mertoggles.togglestates(side,4))
    set(handles.(['togglemedial_',sidestr]),'Value',mertoggles.togglestates(side,5))

side=2; sidestr='left';
set(handles.(['keycentral_',sidestr]),'Value',mertoggles.keycontrol(side,1))
    set(handles.(['keyanterior_',sidestr]),'Value',mertoggles.keycontrol(side,2))
    set(handles.(['keyposterior_',sidestr]),'Value',mertoggles.keycontrol(side,3))
    set(handles.(['keylateral_',sidestr]),'Value',mertoggles.keycontrol(side,4))
    set(handles.(['keymedial_',sidestr]),'Value',mertoggles.keycontrol(side,5))
set(handles.(['togglecentral_',sidestr]),'Value',mertoggles.togglestates(side,1))
    set(handles.(['toggleanterior_',sidestr]),'Value',mertoggles.togglestates(side,2))
    set(handles.(['toggleposterior_',sidestr]),'Value',mertoggles.togglestates(side,3))
    set(handles.(['togglelateral_',sidestr]),'Value',mertoggles.togglestates(side,4))
    set(handles.(['togglemedial_',sidestr]),'Value',mertoggles.togglestates(side,5))

end
setappdata(getappdata(handles.mercontrolfig,'resultfig'),'mertoggles',mertoggles); % also store toggle data in resultfig.

% set togglestates 
% toggles right
side=1;
if get(handles.togglecentral_right,'Value') && isvalid(merhandles.central{side})
    set(merhandles.central{side},'Visible','on')
elseif ~get(handles.togglecentral_right,'Value') && isvalid(merhandles.central{side})
    set(merhandles.central{side},'Visible','off')
end
if get(handles.toggleanterior_right,'Value') && isvalid(merhandles.anterior{side})
    set(merhandles.anterior{side},'Visible','on')
elseif ~get(handles.toggleanterior_right,'Value') && isvalid(merhandles.anterior{side})
    set(merhandles.anterior{side},'Visible','off')
end
if get(handles.toggleposterior_right,'Value') && isvalid(merhandles.posterior{side})
    set(merhandles.posterior{side},'Visible','on')
elseif ~get(handles.toggleposterior_right,'Value') && isvalid(merhandles.posterior{side})
    set(merhandles.posterior{side},'Visible','off')
end
if get(handles.togglelateral_right,'Value') && isvalid(merhandles.lateral{side})
    set(merhandles.lateral{side},'Visible','on')
elseif ~get(handles.togglelateral_right,'Value') && isvalid(merhandles.lateral{side})
    set(merhandles.lateral{side},'Visible','off')
end
if get(handles.togglemedial_right,'Value') && isvalid(merhandles.medial{side})
    set(merhandles.medial{side},'Visible','on')
elseif ~get(handles.togglemedial_right,'Value') && isvalid(merhandles.medial{side})
    set(merhandles.medial{side},'Visible','off')
end
% toggles left
side=2;
sidestr='left';
if get(handles.(['togglecentral_',sidestr]),'Value') && isvalid(merhandles.central{side})
    set(merhandles.central{side},'Visible','on')
elseif ~get(handles.(['togglecentral_',sidestr]),'Value') && isvalid(merhandles.central{side})
    set(merhandles.central{side},'Visible','off')
end
if get(handles.(['toggleanterior_',sidestr]),'Value') && isvalid(merhandles.anterior{side})
    set(merhandles.anterior{side},'Visible','on')
elseif ~get(handles.(['toggleanterior_',sidestr]),'Value') && isvalid(merhandles.anterior{side})
    set(merhandles.anterior{side},'Visible','off')
end
if get(handles.(['toggleposterior_',sidestr]),'Value') && isvalid(merhandles.posterior{side})
    set(merhandles.posterior{side},'Visible','on')
elseif ~get(handles.(['toggleposterior_',sidestr]),'Value') && isvalid(merhandles.posterior{side})
    set(merhandles.posterior{side},'Visible','off')
end
if get(handles.(['togglelateral_',sidestr]),'Value') && isvalid(merhandles.lateral{side})
    set(merhandles.lateral{side},'Visible','on')
elseif ~get(handles.(['togglelateral_',sidestr]),'Value') && isvalid(merhandles.lateral{side})
    set(merhandles.lateral{side},'Visible','off')
end
if get(handles.(['togglemedial_',sidestr]),'Value') && isvalid(merhandles.medial{side})
    set(merhandles.medial{side},'Visible','on')
elseif ~get(handles.(['togglemedial_',sidestr]),'Value') && isvalid(merhandles.medial{side})
    set(merhandles.medial{side},'Visible','off')
end

% set pos handles
if exist('dist','var')
    if side==2
    sidestr='left';
    set(handles.(['poscentral_',sidestr]),'String',num2str(dist));
    set(handles.(['posanterior_',sidestr]),'String',num2str(dist));
    set(handles.(['posposterior_',sidestr]),'String',num2str(dist));
    set(handles.(['poslateral_',sidestr]),'String',num2str(dist));
    set(handles.(['posmedial_',sidestr]),'String',num2str(dist));
    elseif side==1
    sidestr='right';
    set(handles.(['poscentral_',sidestr]),'String',num2str(dist));
    set(handles.(['posanterior_',sidestr]),'String',num2str(dist));
    set(handles.(['posposterior_',sidestr]),'String',num2str(dist));
    set(handles.(['poslateral_',sidestr]),'String',num2str(dist));
    set(handles.(['posmedial_',sidestr]),'String',num2str(dist));
    end
end
setappdata(resultfig,'merhandles',merhandles)