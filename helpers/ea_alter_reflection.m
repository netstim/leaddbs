function ea_alter_reflection(resultfig,factor,reset)
% assuming all patches have same properties, slightly alter them

if ~exist('resultfig','var')
   resultfig=gcf; 
end

if ~exist('factor','var')
    factor=1;
end
allpatches = findall(resultfig,'Type','patch');

if reset
    origpatchlighting=getappdata(resultfig,'origpatchlighting');
    if isempty(origpatchlighting)
        return
    end
    
    [allpatches.AmbientStrength]=deal(origpatchlighting(1));
    [allpatches.DiffuseStrength]=deal(origpatchlighting(2));
    [allpatches.SpecularStrength]=deal(origpatchlighting(3));
    [allpatches.SpecularExponent]=deal(origpatchlighting(4));
    return
end




origpatchlighting=getappdata(resultfig,'origpatchlighting');
if isempty(origpatchlighting) % function has not yet been used, store original values in case of reset
    setappdata(resultfig,'origpatchlighting',...
        [allpatches(1).AmbientStrength,...
        allpatches(1).DiffuseStrength,...
        allpatches(1).SpecularStrength,...
        allpatches(1).SpecularExponent]);
end

% add jitter to settings
[allpatches.AmbientStrength]=deal(maxminvals(allpatches(1).AmbientStrength+(randn(1)*0.001*factor)));
[allpatches.DiffuseStrength]=deal(maxminvals(allpatches(1).DiffuseStrength+(randn(1)*0.001*factor)));
[allpatches.SpecularStrength]=deal(maxminvals(allpatches(1).SpecularStrength+(randn(1)*0.001*factor)));
try
[allpatches.SpecularExponent]=deal(maxvals(allpatches(1).SpecularExponent+(randn(1)*0.1*factor)));
end
function v=maxminvals(v)


v(v>1)=1;
v(v<0)=0;


function v=maxvals(v)
v(v<0)=0;



