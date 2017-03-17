function ea_elvis_surfice(~,~,handles,native)


oexp=ea_exportatlas([],[],'PLY',handles);
if oexp

script=['BEGIN;',...
    ' RESETDEFAULTS;',...
    ' MESHLOAD(''',oexp,''');',...
    'END.'];
ea_surfice(script,0); % no hold

end
