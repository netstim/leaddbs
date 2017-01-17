function ea_listatlassets(handles,mninative)
% mninative==1: MNI space, ==2: native space

options.prefs=ea_prefs('');

as=dir([ea_space(options,'atlases')]);

asc=cell(0);
cnt=1;
for i=1:length(as)
    if as(i).isdir && ~strcmp(as(i).name(1),'.')
        asc{cnt}=as(i).name;
        cnt=cnt+1;
    end
end
options.prefs=ea_prefs('');

if options.prefs.env.dev
    asc{end+1}='Segment patient anatomy';
end
asc{end+1}='Use none';
nasc=cell(0);
if mninative==2
    
    if ~strcmp(get(handles.patdir_choosebox,'String'),'Choose Patient Directory')
       
        % sweep pt dir for atlases
        nas=dir([get(handles.patdir_choosebox,'String'),filesep,'atlases',filesep]);
        cnt=1;
        for i=1:length(nas)
            if nas(i).isdir && ~strcmp(nas(i).name(1),'.')
                nasc{cnt}=['Local atlas: ',nas(i).name];
                cnt=cnt+1;
            end
        end

    end
    
end


set(handles.atlassetpopup,'String',[asc,nasc]);

[~,defix]=ismember(options.prefs.atlases.default,asc);
if defix
    set(handles.atlassetpopup,'Value',defix);
end

