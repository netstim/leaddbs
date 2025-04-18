function conname=ea_conn2connid(conname)

if ~isempty(conname)
    conname = regexprep(conname, '\W', '');
end
