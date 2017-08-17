function d2 = ea_tdhandles2options(tdhandles,fd2)

% defaults:
d2.col_overlay=1;
d2.con_overlay=1;
d2.con_color=[1,1,1];
d2.lab_overlay=1;
d2.bbsize=50;
d2.backdrop='ICBM 152 2009b NLIN Asym T2';

if isempty(tdhandles) % call from somewhere else, return defaults or loaded file
    try
        p = ea_prefs('');
        d2 = p.machine.d2;
    end
end

if nargin>1 % fuse with prior d2 structure
    fieldNames = fieldnames(fd2);
    for i = 1:size(fieldNames,1)
        d2.(fieldNames{i}) = fd2.(fieldNames{i});
    end
end

if isempty(tdhandles)
    return
end

try d2.backdrop=get(tdhandles.tdbackdrop,'String'); end
try d2.backdrop=d2.backdrop{get(tdhandles.tdbackdrop,'Value')}; end
try d2.col_overlay=get(tdhandles.tdcolorscheck,'Value'); end
try	d2.con_overlay=get(tdhandles.tdcontourcheck,'Value'); end
try d2.con_color=getappdata(tdhandles.tdcontourcolor,'color'); end
if isempty(d2.con_color)
    d2.con_color=[1,1,1]; % white
end
try d2.lab_overlay=get(tdhandles.tdlabelcheck,'Value'); end
try d2.bbsize=str2double(get(tdhandles.bbsize,'String')); end
try	d2.fid_overlay=get(tdhandles.tdfidcheck,'Value'); end
