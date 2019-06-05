function [] = ea_defaultview(varargin)
% saves and sets default view preferences
% there must be Electrode-Scene figure
%
% ea_defaultview() saves current view as default view
%
% ea_defaultview(v,togglestates) sets view and togglesates

H = findall(0,'type','figure');
resultfig = H(~cellfun(@isempty,strfind({H(:).Name},{'Electrode-Scene'})));
resultfig = resultfig(1); % take the first if there are many.
togglestates = getappdata(resultfig,'togglestates');
set(0,'CurrentFigure',resultfig);

if nargin == 0
    % save current view and togglesates
    ea_setprefs('view',ea_view)
    ea_setprefs('togglestates',getappdata(resultfig,'togglestates'))

elseif nargin == 2
    % set preferences specified in vararg in
    % togglestates
    togglestates.xyzmm = varargin{2}.xyzmm;
    togglestates.xyztoggles = varargin{2}.xyztoggles;
    togglestates.xyztransparencies = varargin{2}.xyztransparencies;
    togglestates.refreshview = 1;
    ea_anatomyslices(resultfig,togglestates,struct,[]);
    % camera view
    ea_view(varargin{1});
    % update togglestates
    setappdata(resultfig,'togglestates',togglestates);
    % update anatomy control
    close(getappdata(resultfig,'awin'));
    options = getappdata(resultfig,'options');
    awin = ea_anatomycontrol(resultfig,options);
    setappdata(resultfig,'awin',awin);
end

end
