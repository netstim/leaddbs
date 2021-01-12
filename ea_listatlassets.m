function ea_listatlassets(options,handles,mninative,oldatlas)
% mninative==1: MNI space, ==2: native space

if ~exist('oldatlas','var')
    oldatlas='';
end

% dir 'atlases' folder
atlases=dir(ea_space(options,'atlases'));
atlases = {atlases(cell2mat({atlases.isdir})).name};    % only keep folders
atlases = atlases(cellfun(@(x) ~strcmp(x(1),'.'), atlases));  % also remove '.', '..' and '.*' folders from dir results

if options.prefs.env.dev
    atlases{end+1}='Segment patient anatomy';
end

if ~(isfield(handles, 'output') && strcmp(get(handles.output,'Name'),'FEM-based VAT model setting'))
    atlases{end+1}='Use none';
end

natlases=cell(0);
if mninative==2
    if ~strcmp(get(handles.patdir_choosebox,'String'),'Choose Patient Directory')
        % sweep pt dir for atlases
        natlases=dir([get(handles.patdir_choosebox,'String'),filesep,'atlases',filesep]);
        natlases = {natlases(cell2mat({natlases.isdir})).name};
        natlases = natlases(cellfun(@(x) ~strcmp(x(1),'.'), natlases));
        natlases = cellfun(@(x) {['Local atlas: ', x]}, natlases);
    end
end

try
    set(handles.atlassetpopup,'String',[atlases,natlases]);
catch
    keyboard
end

[~,defix]=ismember(options.prefs.machine.atlases.default,[atlases,natlases]);
if defix
    set(handles.atlassetpopup,'Value',defix);
end

if ~isempty(oldatlas)
    [~,oldix]=ismember(oldatlas,[atlases,natlases]);
    if oldix
        set(handles.atlassetpopup,'Value',oldix);
    end
end
