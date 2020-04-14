function str=ea_space2sub(str) % replaces subscores with spaces
if ~iscell(str)
str(str==' ')='_';
else
    for c=1:length(str)
        
      str{c}(str{c}==' ')='_';
  
    end
end