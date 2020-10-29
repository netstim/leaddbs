function str=ea_sub2space(str) % replaces subscores with spaces
if ~iscell(str)
    str(str=='_')=' ';
else
    for c=1:length(str)
    	str{c}(str{c}=='_')=' ';
    end
end
