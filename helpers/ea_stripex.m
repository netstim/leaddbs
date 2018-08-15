function fn=ea_stripex(fn)
if iscell(fn)
   for f=1:length(fn)
       
       [~,fn{f}]=fileparts(fn{f});
       [~,fn{f}]=fileparts(fn{f}); % .gz support

   end
else
[~,fn]=fileparts(fn);
[~,fn]=fileparts(fn); % .gz support
end