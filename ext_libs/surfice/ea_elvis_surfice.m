function ea_elvis_surfice(~,~,handles,native)


uipatdir=getappdata(handles.leadfigure,'uipatdir');
 ea_exportpat([],[],'PLY',handles)
for pt=1:length(uipatdir)
   thispt=uipatdir{pt};
   script=['BEGIN;',...
       ' RESETDEFAULTS;',...
       ' MESHLOAD(''',[thispt,filesep,'export',filesep,'ply',filesep,'combined_scene.ply'],''');',...
       'END.'];
   ea_surfice(script,0); % no hold
end


