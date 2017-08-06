function ea_listatlassets(options,handles,mninative,oldatlas)
% mninative==1: MNI space, ==2: native space

if ~exist('oldatlas','var')
    oldatlas='';
end

atlases=struct2cell(dir(ea_space(options,'atlases')));  % dir 'atlases' folder
atlases=atlases(1,cell2mat(atlases(5,:)));  % only keep folders, the 5th row of the struct2cell result is 'isdir'
atlases=atlases(cellfun(@(x) ~strcmp(x(1),'.'),atlases)); % remove '.', '..' and '.*' folders from dir results

if options.prefs.env.dev
    atlases{end+1}='Segment patient anatomy';
end

atlases{end+1}='Use none';

natlases=cell(0);
if mninative==2
    if ~strcmp(get(handles.patdir_choosebox,'String'),'Choose Patient Directory')
        % sweep pt dir for atlases
        natlases=struct2cell(dir([get(handles.patdir_choosebox,'String'),filesep,'atlases',filesep]));
        natlases=natlases(1,cell2mat(natlases(5,:)));
        natlases=natlases(cellfun(@(x) ~strcmp(x(1),'.'),natlases));
        natlases=cellfun(@(x) ['Local atlas: ', x],natlases,'Uniform',0);
    end
end

try
    set(handles.atlassetpopup,'String',[atlases,natlases]);
catch
    keyboard
end

[~,defix]=ismember(options.prefs.atlases.default,[atlases,natlases]);
if defix
    set(handles.atlassetpopup,'Value',defix);
end

if ~isempty(oldatlas)
    [~,oldix]=ismember(oldatlas,[atlases,natlases]);
    if oldix
        set(handles.atlassetpopup,'Value',oldix);
    end
end
