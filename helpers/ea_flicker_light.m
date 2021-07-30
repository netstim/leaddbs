function ea_flicker_light(resultfig,factor,reset)

if ~exist('factor','var')
    factor=1;
end
if ~exist('reset','var')
    reset=0;
end
lights=findall(resultfig,'Type','light');

if reset
    origlightpos=getappdata(resultfig,'origlightpos');
    if isempty(origlightpos)
        return
    end
    
    
    origlightcol=getappdata(resultfig,'origlightcol');
    for l=1:length(lights)
        lights(l).Position=origlightpos{l};
        lights(l).Color=origlightcol{l};
    end
    return
end


origlightpos=getappdata(resultfig,'origlightpos');
if isempty(origlightpos) % function has not yet been used
    setappdata(resultfig,'origlightpos',{lights.Position});
    setappdata(resultfig,'origlightcol',{lights.Color});
end

% now jitter settings a bit

for l=1:length(lights)
    lights(l).Position=lights(l).Position+lights(l).Position.*(randn(1,3)*0.1*factor);
    lights(l).Color=maxvals(lights(l).Color+(randn(1,3)*0.001*factor));
end



function v=maxvals(v)


v(v>1)=1;
v(v<0)=0;