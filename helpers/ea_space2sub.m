function str=ea_space2sub(str) % replaces subscores with spaces
if ~iscell(str)
str(str==' ')='_';
str(str==')')='';
str(str=='(')='';
else
    for c=1:length(str)
        
      str{c}(str{c}==' ')='_';
      str{c}(str{c}=='(')='';
      str{c}(str{c}==')')='';
      
    end
end