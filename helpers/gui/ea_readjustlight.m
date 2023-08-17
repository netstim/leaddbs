function ea_readjustlight(hfig)

% get light handles
RightLight = getappdata(hfig, 'RightLight');
LeftLight = getappdata(hfig, 'LeftLight');
CeilingLight = getappdata(hfig, 'CeilingLight');
CamLight = getappdata(hfig, 'CamLight');

prefs = ea_prefs;

try
    camlight(CamLight, 'headlight'); % move light object.
end

try
    set(CeilingLight, 'Position', [0 0 10], 'style', 'local', 'Color', prefs.d3.ceilinglightcolor); % not modifiable, infinite light.
end

try
    set(RightLight, 'Position', [-100 0 0], 'style', 'infinite', 'Color', prefs.d3.rightlightcolor); % not modifiable, infinite light.
end

try
    set(LeftLight, 'Position', [100 0 0], 'style', 'infinite', 'Color', prefs.d3.leftlightcolor); % not modifiable, infinite light.
end
