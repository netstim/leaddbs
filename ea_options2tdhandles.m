function d2=ea_options2tdhandles(tdhandles,d2)

try   
    allbackdrops=get(tdhandles.tdbackdrop,'String');
    [~,bdix]=ismember(d2.backdrop,allbackdrops);
    set(tdhandles.tdbackdrop,'Value',bdix);

    if get(tdhandles.tdbackdrop,'Value')>length(allbackdrops) || get(tdhandles.tdbackdrop,'Value')<1
        set(tdhandles.tdbackdrop,'Value',1);
    end

    set(tdhandles.tdcolorscheck,'Value',d2.col_overlay);
    set(tdhandles.tdcontourcheck,'Value',d2.con_overlay);
    setappdata(tdhandles.tdcontourcolor,'color',d2.con_color);
    set(tdhandles.tdlabelcheck,'Value',d2.lab_overlay);
    set(tdhandles.bbsize,'String',num2str(d2.bbsize));

    try % additional values when called from lead_anatomy
        set(tdhandles.tracor,'Value',d2.tracor);
        set(tdhandles.depth,'String',num2str(d2.depth));
    end

    try % fiducials check
        set(tdhandles.tdfidcheck,'Value',d2.fid_overlay);
    end
end
