function txt=ea_sub2space(txt) % replaces subscores with spaces
if ~iscell(txt)
    txt(txt=='_')=' ';
else
    for c=1:length(txt)
        txt{c}(txt{c}=='_')=' ';
    end
end